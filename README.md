# miRBooking

Implementation of the miRBooking algorithm and metrics in C

 - fast and memory efficient
 - usable from Python, JavaScript and Vala via GObject introspection
 - usable from Java via JNA
 - memory-mapped score tables, target and mirnas (for parallel execution)
 - memory-mapped and zero-copy for input files containing sequences (i.e. FASTA)
 - binary with support for static linking for more portability
 - stdin/stdout for piping from and into other tools

## Usage

```bash
mirbooking --mirnas mirbase.fa --targets human-genome.fa --score-table scores [--output output.tsv] [--threshold] [--log-base] [quantities]
```

## Parallelization (using GNU parallel)

```bash
parallel mirbooking --mirnas mature.fa --targets hg38.fa --score-table scores ::: wildtype.tsv over-expression.tsv
```

The targets, mirnas and score table files will reuse the same physical memory
across all processes, yielding a minimal memory usage.

