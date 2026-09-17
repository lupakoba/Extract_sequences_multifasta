# extract_gene_multifasta.py — Changelog & Usage

## 1. What the script does

Extracts, from multiple assembled genomes (FASTA), the region homologous to a
reference gene (e.g. *E. coli* `cirA`), using **local BLASTN**, and generates a
**multifasta** containing:

1. The reference sequence (always first).
2. One sequence per genome, extracted from that genome.

The extracted sequence is delivered **raw** (not translated, stops and
frameshifts not corrected), so that a downstream script can analyze amino acid
changes, premature stops, or frameshifts.

## 2. Changes made in this revision

**Problem being solved:** running the script with different `--ref` files
(different genes) against the same genomes caused each run's results to
overwrite the previous one, because the output path was set manually
(`--output`).

**Solution implemented:**

- Removed the `--output` and `--report` arguments (manual file paths).
- Added `--results_dir` (root results folder, default `resultados/`).
- Added `--gene_name` (optional; if not given, the sequence `id` from `--ref`
  is used instead).
- The script now automatically creates a **subfolder per gene** inside
  `--results_dir`, and all files from that run (multifasta, TSV report, and
  the raw blast TSVs if `--keep_blast_tsv` is used) are placed inside that
  subfolder.

**Resulting folder structure:**

```
resultados/                      <- --results_dir (default: "resultados")
├── cirA/
│   ├── cirA_multifasta.fasta
│   ├── extraction_report.tsv
│   └── blast_raw_hits/          <- only if --keep_blast_tsv is used
│       ├── genoma1.blast.tsv
│       └── genoma2.blast.tsv
├── fepA/
│   ├── fepA_multifasta.fasta
│   └── extraction_report.tsv
└── iroN/
    ├── iroN_multifasta.fasta
    └── extraction_report.tsv
```

With this, running the script multiple times with different genes on the same
genomes no longer overwrites previous results; each gene ends up neatly
organized in its own folder.

### Some additional features

- `blastn -task blastn` with a low `word_size` and `-dust no`: higher
  sensitivity to divergence, indels, and low-complexity regions, so the hit
  isn't lost when there's a premature stop or a frameshift.
- Hits are ordered by their position in the **reference** (not in the genome),
  which allows reconstructing genes fragmented **within the same contig**
  (due to a large indel) or **across two or more contigs** (when the assembly
  breaks the gene apart).
- Automatic reverse-complement of fragments that fall on the `-` strand.
- Removal of redundant/nested hits (the one with the highest bitscore is
  kept).
- Per-genome TSV report with status (`OK`, `FRAGMENTADO_MISMO_CONTIG`,
  `FRAGMENTADO_ENTRE_CONTIGS`, `SIN_HIT`), number of fragments, contigs
  involved, and % coverage relative to the reference.

## 3. Requirements

- **BLAST+** installed and available in the `PATH` (the `blastn` command).
  Typical install: `conda install -c bioconda blast` or
  `apt install ncbi-blast+`.
- **Biopython**: `pip install biopython`.

## 4. Usage

### Basic usage

```bash
python extract_gene_multifasta.py \
    --ref cirA_reference.fasta \
    --genomes_dir genomas/
```

This generates:

```
resultados/cirA/cirA_multifasta.fasta
resultados/cirA/extraction_report.tsv
```

### Multiple genes against the same genomes

```bash
python extract_gene_multifasta.py --ref cirA_reference.fasta --genomes_dir genomas/
python extract_gene_multifasta.py --ref fepA_reference.fasta --genomes_dir genomas/
python extract_gene_multifasta.py --ref iroN_reference.fasta --genomes_dir genomas/
```

Each run is placed in its own subfolder (`resultados/cirA/`,
`resultados/fepA/`, `resultados/iroN/`), without overwriting the previous
ones.

### Custom parameters

```bash
python extract_gene_multifasta.py \
    --ref cirA_reference.fasta \
    --genomes_dir genomas/ \
    --results_dir mis_resultados \
    --gene_name cirA \
    --min_identity 65 \
    --min_cov 15 \
    --keep_blast_tsv
```

### Available arguments

| Argument | Required | Default | Description |
|---|---|---|---|
| `--ref` | Yes | — | FASTA with the reference gene sequence (a single sequence). |
| `--genomes_dir` | Yes | — | Directory with genomes `.fasta` / `.fa` / `.fna` (one per file). |
| `--results_dir` | No | `resultados` | Root folder where per-gene subfolders are created. |
| `--gene_name` | No | id of the sequence in `--ref` | Name used for the subfolder and output files. |
| `--min_identity` | No | `70.0` | Minimum identity (%) per hit to be considered valid. |
| `--min_cov` | No | `20.0` | Minimum length of each hit relative to the query (%). Lower it if you expect hits to be heavily fragmented by large indels. |
| `--keep_blast_tsv` | No | disabled | Keeps the raw blast TSV per genome (useful for manually reviewing edge cases or possible paralogs). |

## 5. Interpreting the report (`extraction_report.tsv`)

| Column | Meaning |
|---|---|
| `genoma` | Genome file name (without extension). |
| `status` | `OK` (single fragment), `FRAGMENTADO_MISMO_CONTIG`, `FRAGMENTADO_ENTRE_CONTIGS`, or `SIN_HIT`. |
| `n_fragmentos` | Number of fragments (HSPs) used to reconstruct the sequence. |
| `contigs` | Source contig(s), in the order they were concatenated. |
| `cobertura_query_%` | Percentage of the reference sequence covered by the extracted fragments. |

It's recommended to manually review genomes with a `status` other than `OK`,
especially those with low `cobertura_query_%` or `n_contigs > 1`, since these
could reflect either genuinely fragmented genes or assembly artifacts.

## 6. Known limitations

- Does not correct frameshifts or fill gaps: it delivers the raw sequence as
  found in the genome (which is exactly the intended input for the downstream
  amino-acid-change analysis script).
- If the gene is duplicated (paralogs) within a single genome, the colinear
  chain with the best bitscore relative to the reference is kept; alternative
  hits are discarded. Use `--keep_blast_tsv` to review these cases manually.
- Uses `blastn -subject` (no indexed database) for simplicity; for hundreds
  or thousands of genomes, it's more efficient to build a database with
  `makeblastdb` per genome and parallelize the runs.

