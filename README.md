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

The command line program expects a number of inputs:

 - a FASTA containing mRNA transcripts where the identifier is the accession
   (i.e. NM_002710.3)
 - a FASTA containing mature miRNAs where the first token in the comment is the
   accession (i.e. MIMAT0004792)
 - a two columns TSV document mapping target accessions to coding regions the
   'a..b' format where a and b are inclusive 1-based indexes
 - a score table with hybridization probabilities for seeds
 - a quantity file mapping target and mirna accessions to expressed quantity in
   any comparable units

The only requirements for the quantities is to be in the same magnitude order.

The output is a TSV with the following columns:

 - Target Accession
 - MiRNA Accession
 - Position
 - Location
 - Probability
 - Occupancy
 - Silencing
 
## Installation

You'll need [Meson](http://mesonbuild.com/) and [Ninja](http://ninja-build.org/) 
as well as GLib development files installed on your system.

```bash
mkdir build
meson --buildtype=release
ninja
ninja install
```

You can perform a local installation using `meson --prefix=$HOME/.local`, but you'll need `LD_LIBRARY_PATH` set accordingly since the `mirbooking` program uses a shared library. Otherwise, a static linkage can be done by calling `meson --default-library=static`.

## Parallelization (using GNU parallel)

```bash
parallel mirbooking --mirnas mature.fa --targets hg38.fa --score-table scores ::: wildtype.tsv over-expression.tsv
```

The targets, mirnas and score table files will reuse the same physical memory
across all processes, yielding a minimal memory usage.

