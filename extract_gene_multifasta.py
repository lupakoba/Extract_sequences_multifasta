#!/usr/bin/env python3
"""
extract_gene_multifasta.py

Extrae secuencias homologas a un gen de referencia (ej: cirA de E. coli)
desde multiples genomas ensamblados (FASTA), usando BLASTN local, y genera
un multifasta con la referencia como PRIMERA secuencia, seguida de la
secuencia extraida de cada genoma.

DISEÑO (por que funciona asi):
- Se usa blastn con -task blastn (no megablast) y word_size bajo, que es
  mas sensible a divergencia/indels/SNPs que el modo por defecto. Esto es
  clave para no perder el hit cuando hay un stop prematuro o un frameshift,
  ya que blastn segmenta el alineamiento en HSPs (fragmentos) en vez de
  fallar en encontrar el gen completo.
- NO se traduce ni se filtra por marco de lectura abierto: la idea es traer
  la secuencia nucleotidica cruda, stops/frameshifts incluidos, para que tu
  script de analisis de cambios aminoacidicos la procese despues.
- Si el gen aparece fragmentado en varios HSPs (por un indel grande, una
  region muy divergente, o porque el ensamblaje corta el gen entre dos
  contigs), el script ORDENA los fragmentos segun su posicion en la
  secuencia de referencia (query) -no en el genoma- y los concatena. Asi
  reconstruye el gen "en el orden correcto" sin importar si viene de 1 o
  varios contigs.
- Cada hit puede quedar en hebra + o -; se hace reverse-complement segun
  corresponda antes de concatenar.
- Se eliminan hits redundantes/anidados (mismo tramo del query cubierto por
  mas de un hit) quedandose con el de mayor bitscore.
- Se genera un reporte TSV indicando, por genoma: si el gen se encontro,
  en cuantos fragmentos, en que contig(s), y que % del query esta cubierto
  -para que puedas revisar manualmente los casos raros (multi-contig,
  cobertura baja, etc).

LIMITACIONES A TENER EN CUENTA:
- Esto NO corrige frameshifts ni rellena gaps: entrega la secuencia cruda
  extraida del genoma, que es justo lo que pediste para tu script aguas
  abajo.
- Si el gen esta duplicado (parálogos) en el mismo genoma, este script se
  queda con la mejor cadena colineal respecto al query; los hits
  redundantes/alternativos no se usan. Si sospechas duplicaciones reales,
  revisa el TSV crudo de blast (se guarda temporalmente; puedes activar
  --keep_blast_tsv para conservarlo).
- Usa blastn -subject (sin crear una BD indexada) por simplicidad; para
  cientos/miles de genomas es mas eficiente crear una BD con makeblastdb
  por genoma y paralelizar. Se deja indicado en el codigo donde cambiarlo.

Requiere:
    - BLAST+ instalado y en el PATH (blastn). Ej: conda install -c bioconda blast
    - biopython (pip install biopython)

Uso basico (crea resultados/cirA/cirA_multifasta.fasta, etc.):
    python extract_gene_multifasta.py \
        --ref cirA_reference.fasta \
        --genomes_dir genomas/

Corriendo varios genes sobre los MISMOS genomas (cada uno en su propia
subcarpeta dentro de --results_dir, sin pisarse entre si):
    python extract_gene_multifasta.py --ref cirA_reference.fasta --genomes_dir genomas/
    python extract_gene_multifasta.py --ref fepA_reference.fasta --genomes_dir genomas/
    python extract_gene_multifasta.py --ref iroN_reference.fasta --genomes_dir genomas/

    -> resultados/cirA/cirA_multifasta.fasta, resultados/cirA/extraction_report.tsv
    -> resultados/fepA/fepA_multifasta.fasta, resultados/fepA/extraction_report.tsv
    -> resultados/iroN/iroN_multifasta.fasta, resultados/iroN/extraction_report.tsv

Uso con parametros ajustados y carpeta raiz / nombre de gen personalizados:
    python extract_gene_multifasta.py \
        --ref cirA_reference.fasta \
        --genomes_dir genomas/ \
        --results_dir mis_resultados \
        --gene_name cirA \
        --min_identity 65 \
        --min_cov 15 \
        --keep_blast_tsv
"""

import argparse
import os
import subprocess
import sys
import tempfile
import shutil
import logging
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)


def check_blast_installed():
    if shutil.which("blastn") is None:
        log.error("No se encontro 'blastn' en el PATH. Instala BLAST+ "
                   "(ej: conda install -c bioconda blast) e intenta de nuevo.")
        sys.exit(1)


def run_blastn(query_fasta, subject_fasta, out_tsv, evalue=1e-10, word_size=11):
    """
    Corre blastn: query = gen de referencia, subject = genoma completo.
    -task blastn (no megablast) + word_size bajo = mayor sensibilidad a
    secuencias divergentes/con indels, a costa de un poco mas de tiempo.
    -dust no: no enmascara regiones de baja complejidad (podrian ser parte
    real del gen).
    """
    cmd = [
        "blastn",
        "-query", str(query_fasta),
        "-subject", str(subject_fasta),
        "-out", str(out_tsv),
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen "
                   "qstart qend sstart send evalue bitscore qlen slen",
        "-evalue", str(evalue),
        "-word_size", str(word_size),
        "-dust", "no",
        "-max_target_seqs", "50",
        "-task", "blastn",
    ]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)


def parse_hits(tsv_path, min_identity, min_len_frac):
    """Filtra hits por identidad minima y longitud minima relativa al query."""
    hits = []
    if not os.path.exists(tsv_path) or not os.path.getsize(tsv_path):
        return hits
    with open(tsv_path) as fh:
        for line in fh:
            f = line.strip().split("\t")
            (qseqid, sseqid, pident, length, mismatch, gapopen,
             qstart, qend, sstart, send, evalue, bitscore, qlen, slen) = f
            pident = float(pident)
            length = int(length)
            qlen = int(qlen)
            if pident < min_identity:
                continue
            if (length / qlen) * 100 < min_len_frac:
                continue
            hits.append({
                "qseqid": qseqid, "sseqid": sseqid, "pident": pident,
                "length": length, "qstart": int(qstart), "qend": int(qend),
                "sstart": int(sstart), "send": int(send),
                "evalue": float(evalue), "bitscore": float(bitscore),
                "qlen": qlen, "slen": int(slen),
            })
    return hits


def remove_redundant_hits(hits, overlap_thresh=0.6):
    """
    Descarta hits cuyo tramo de QUERY se solapa fuertemente con uno de
    mayor bitscore (evita duplicar tramos del gen por hits alternativos/
    parálogos/repeticiones).
    """
    hits_sorted = sorted(hits, key=lambda h: -h["bitscore"])
    kept = []
    for h in hits_sorted:
        qs, qe = sorted([h["qstart"], h["qend"]])
        redundant = False
        for k in kept:
            ks, ke = sorted([k["qstart"], k["qend"]])
            inter = max(0, min(qe, ke) - max(qs, ks))
            shorter = min(qe - qs, ke - ks) or 1
            if inter / shorter > overlap_thresh:
                redundant = True
                break
        if not redundant:
            kept.append(h)
    return kept


def order_hits_by_query(hits):
    """Ordena los hits (ya no redundantes) segun su posicion en la referencia,
    sin importar de que contig del genoma provienen. Esto es lo que permite
    reconstruir un gen fragmentado entre 2+ contigs, o partido por un indel
    grande dentro del mismo contig."""
    return sorted(hits, key=lambda h: min(h["qstart"], h["qend"]))


def extract_fragment(genome_records, sseqid, sstart, send):
    rec = genome_records[sseqid]
    lo, hi = sorted([sstart, send])
    sub = str(rec.seq[lo - 1:hi])  # coords blast son 1-based inclusivas
    if sstart > send:
        sub = str(Seq(sub).reverse_complement())
    return sub


def process_genome(ref_fasta, genome_fasta, workdir, min_identity, min_len_frac,
                    keep_blast_tsv, output_tsv_dir):
    genome_name = Path(genome_fasta).stem
    out_tsv = workdir / f"{genome_name}.blast.tsv"
    run_blastn(ref_fasta, genome_fasta, out_tsv)

    if keep_blast_tsv:
        dest = Path(output_tsv_dir) / f"{genome_name}.blast.tsv"
        shutil.copy(out_tsv, dest)

    hits = parse_hits(out_tsv, min_identity, min_len_frac)
    if not hits:
        log.warning(f"[{genome_name}] Sin hits significativos para el gen de referencia.")
        return None, {"genome": genome_name, "status": "SIN_HIT", "n_fragmentos": 0,
                       "contigs": [], "cobertura_query_%": 0}

    hits = remove_redundant_hits(hits)
    ordered = order_hits_by_query(hits)

    genome_records = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))

    fragments = []
    contigs_used = []
    total_query_bp = 0
    for h in ordered:
        seq_frag = extract_fragment(genome_records, h["sseqid"], h["sstart"], h["send"])
        fragments.append(seq_frag)
        contigs_used.append(h["sseqid"])
        total_query_bp += abs(h["qend"] - h["qstart"]) + 1

    combined_seq = "".join(fragments)
    qlen = ordered[0]["qlen"]
    coverage_pct = round(100 * total_query_bp / qlen, 1)

    n_contigs = len(set(contigs_used))
    if n_contigs > 1:
        status = "FRAGMENTADO_ENTRE_CONTIGS"
    elif len(ordered) > 1:
        status = "FRAGMENTADO_MISMO_CONTIG"
    else:
        status = "OK"

    contigs_unicos = list(dict.fromkeys(contigs_used))  # preserva orden, sin duplicados
    info = {
        "genome": genome_name,
        "status": status,
        "n_fragmentos": len(ordered),
        "contigs": contigs_unicos,
        "cobertura_query_%": coverage_pct,
    }

    header = (f"{genome_name} | status={status} | fragmentos={len(ordered)} "
              f"| contigs={','.join(contigs_unicos)} | cobertura_query={coverage_pct}%")

    record = SeqRecord(Seq(combined_seq), id=genome_name, description=header)
    return record, info


def sanitize_name(name):
    """Limpia el nombre del gen/ref para usarlo como nombre de carpeta/archivo."""
    keep = "-_."
    return "".join(c if c.isalnum() or c in keep else "_" for c in name).strip("_")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ref", required=True,
                     help="FASTA con la secuencia de referencia del gen (una sola secuencia)")
    ap.add_argument("--genomes_dir", required=True,
                     help="Directorio con genomas .fasta/.fa/.fna (uno por archivo)")
    ap.add_argument("--results_dir", default="resultados",
                     help="Carpeta raiz donde se guardan los resultados de TODAS las corridas "
                          "(default: 'resultados'). Dentro se crea automaticamente una "
                          "subcarpeta con el nombre del gen (de --ref o --gene_name), asi puedes "
                          "correr distintos --ref sobre los mismos genomas sin pisar resultados.")
    ap.add_argument("--gene_name", default=None,
                     help="Nombre del gen, usado para la subcarpeta y archivos de salida "
                          "(ej: cirA). Si no se indica, se usa el id de la secuencia en --ref.")
    ap.add_argument("--min_identity", type=float, default=70.0,
                     help="Identidad minima %% por hit (default 70)")
    ap.add_argument("--min_cov", type=float, default=20.0,
                     help="Longitud minima de cada hit relativa al query, en %% (default 20). "
                          "Bajalo si esperas hits muy fragmentados por indels grandes.")
    ap.add_argument("--keep_blast_tsv", action="store_true",
                     help="Conserva el TSV crudo de blast por genoma (para revision manual "
                          "de casos raros / posibles parálogos)")
    args = ap.parse_args()

    check_blast_installed()

    genomes_dir = Path(args.genomes_dir)
    genome_files = sorted([p for p in genomes_dir.iterdir()
                            if p.suffix.lower() in (".fasta", ".fa", ".fna")])
    if not genome_files:
        log.error(f"No se encontraron genomas .fasta/.fa/.fna en {genomes_dir}")
        sys.exit(1)

    ref_records = list(SeqIO.parse(args.ref, "fasta"))
    if len(ref_records) != 1:
        log.error("El archivo --ref debe contener EXACTAMENTE una secuencia (la referencia del gen).")
        sys.exit(1)
    ref_record = ref_records[0]

    # nombre del gen -> nombre de la subcarpeta de resultados
    gene_name = sanitize_name(args.gene_name if args.gene_name else ref_record.id)

    ref_record.id = "REFERENCIA_" + ref_record.id
    ref_record.description = "secuencia de referencia"

    # Estructura de carpetas:
    #   <results_dir>/<gen>/<gen>_multifasta.fasta
    #   <results_dir>/<gen>/extraction_report.tsv
    #   <results_dir>/<gen>/blast_raw_hits/*.tsv   (solo si --keep_blast_tsv)
    gene_outdir = Path(args.results_dir) / gene_name
    gene_outdir.mkdir(parents=True, exist_ok=True)

    output_fasta = gene_outdir / f"{gene_name}_multifasta.fasta"
    report_path = gene_outdir / "extraction_report.tsv"
    blast_tsv_dir = gene_outdir / "blast_raw_hits"
    if args.keep_blast_tsv:
        blast_tsv_dir.mkdir(exist_ok=True)

    log.info(f"Resultados de este gen se guardaran en: {gene_outdir}/")

    out_records = [ref_record]
    report_rows = []

    with tempfile.TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        for gf in genome_files:
            log.info(f"Procesando {gf.name} ...")
            rec, info = process_genome(args.ref, gf, workdir,
                                        args.min_identity, args.min_cov,
                                        args.keep_blast_tsv, blast_tsv_dir)
            if rec is not None:
                out_records.append(rec)
            report_rows.append(info)

    SeqIO.write(out_records, output_fasta, "fasta")
    log.info(f"Multifasta generado: {output_fasta} ({len(out_records)} secuencias)")

    with open(report_path, "w") as fh:
        fh.write("genoma\tstatus\tn_fragmentos\tcontigs\tcobertura_query_%\n")
        for r in report_rows:
            fh.write(f"{r['genome']}\t{r['status']}\t{r['n_fragmentos']}\t"
                      f"{';'.join(r['contigs'])}\t{r['cobertura_query_%']}\n")
    log.info(f"Reporte de extraccion: {report_path}")

    n_frag = sum(1 for r in report_rows if r["status"] != "OK" and r["status"] != "SIN_HIT")
    n_missing = sum(1 for r in report_rows if r["status"] == "SIN_HIT")
    if n_frag:
        log.warning(f"{n_frag} genoma(s) con el gen fragmentado (revisa '{report_path}').")
    if n_missing:
        log.warning(f"{n_missing} genoma(s) sin hit detectado para el gen.")


if __name__ == "__main__":
    main()
