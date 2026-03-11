# ScaffMatch-UD Genome Scaffolding

ScaffMatch-UD is a Python implementation of a lightweight genome assembly/scaffolding workflow for unbalanced sequencing datasets. The project builds a de Bruijn graph from FASTQ reads, applies simple graph-cleaning heuristics, extracts contigs, and reports assembly quality metrics.

## What the pipeline does

The main pipeline (`src/pipeline.py`) performs the following steps:

1. Loads FASTQ reads and skips any read containing `N`.
2. Counts encoded k-mers (`k=21` by default).
3. Builds a weighted directed de Bruijn graph.
4. Cleans graph artifacts using:
   - tip removal
   - bubble removal
5. Extracts contigs by traversing non-branching graph paths.
6. Computes assembly metrics:
   - total contigs
   - average contig length
   - max contig length
   - N50
7. Saves a contig length histogram to `output/contig_length_distribution.png`.

## Repository structure

```text
src/
  pipeline.py            # End-to-end assembly pipeline
  read_loader.py         # FASTQ reader (filters reads with N)
  kmer_counter.py        # 2-bit k-mer encoding and counting
  graph_builder.py       # de Bruijn graph construction
  graph_cleaner.py       # tip and bubble cleanup
  contig_extractor.py    # contig path extraction
  metrics.py             # contig stats and N50
  plotter.py             # contig length histogram
  experiment_plots.py    # static experiment/result figure generation
requirements.txt
README.md
```

## Requirements

- Python 3.10+
- Dependencies in `requirements.txt`

Install dependencies:

```bash
pip install -r requirements.txt
```

## How to run

From the repository root:

```bash
python src/pipeline.py
```

By default, the pipeline expects input at:

- `dataset/sample.fastq`

If your FASTQ file lives elsewhere, update `file_path` in `src/pipeline.py`.

## Outputs

The pipeline prints progress and summary stats to stdout and writes plots to `output/`:

- `contig_length_distribution.png`
- `debruijn_graph_sample.png` (only when optional graph visualization is enabled)

The experimental plotting script can also generate additional figures:

```bash
python src/experiment_plots.py
```

This writes multiple `fig_*.png` files under `output/`.

## Notes

- Graph visualization is disabled by default (`ENABLE_GRAPH_VISUALIZATION = False`) to keep runtime manageable.
- Low-support edges are filtered during graph construction (`count < 2` is skipped).
- Progress bars are shown for read processing, graph building, and graph traversal via `tqdm`.
