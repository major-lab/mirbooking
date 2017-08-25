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
mirbooking --mirnas mirbase.fa
           --targets human-genome.fa
           --score-table scores
           [--cds-regions cds-regions.tsv]
           [--threshold 0.0179]
           [--log-base 512]
           [--5prime-footprint 26]
           [--3prime-footprint 19]
           [--5prime-utr-multiplier 0.1]
           [--cds-multiplier 0.1]
           [--3prime-utr-multiplier 1.0]
           [--quantities quantities.tsv]
           [--output output.tsv]
```

To obtain detailed usage and options, launch `mirbookng --help`.

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

| Column      | Description                                   |
| ----------- | --------------------------------------------- |
| target      | Target accession with version                 |
| mirna       | miRNA accession                               |
| position    | Site position on the target                   |
| location    | Location (e.g. 5', CDS, 3' or N/A if unknown) |
| probability | Hybridization probability                     |
| occupancy   | Number of occupied sites at this position     |
| vacancy     | Number of vacant sites at this position       |
| silencing   | Induced silencing                             |

## Installation

You'll need [Meson](http://mesonbuild.com/) and [Ninja](http://ninja-build.org/)
as well as GLib development files installed on your system.

```bash
mkdir build
meson --buildtype=release
ninja
ninja install
```

To generate fast code, configure with `CFLAGS='-Ofast -march=native meson'`.

You can perform a local installation using `meson --prefix=$HOME/.local`, but
you'll need `LD_LIBRARY_PATH` set accordingly since the `mirbooking` program
uses a shared library. Otherwise, a static linkage can be done by calling
`meson --default-library=static`.

To generate introspection metadata, use `meson -Dwith_introspection=true`. To
generate Vala bindings, use `meson -Dwith_vapi=true`.

## Other tools

In addition to the `mirbooking` binary, this package ship a number of utilities
to process quantity files and compute summaries based on Pandas.

The `mirbooking-calibrate` tool is expecting a transcript and miRNA
quantification (e.g. two-column TSV document mapping accession to quantity) and
process it such that it contains approximately the same amount of each kind by
rescaling the smallest one toward the biggest one. It emits a calibrated output
suitable for `mirbooking`.

```bash
mirbooking-calibrate targets.tsv mirnas.tsv | mirbooking [...]
```

The `mirbooking-aggregate` tool perform an sum aggregation on both the
`occupancy` and `silencing` columns. This can be used for various scenarios:

Telling how many duplex have been formed per target and their effective
silencing:

```bash
mirbooking [...] | mirbooking-aggregate target
```

Telling how many duplex have been formed per miRNA and their induced silencing:

```bash
mirbooking [...] | mirbooking-aggregate mirna
```

Telling how effectively miRNAs are occupying their targets:

```bash
mirbooking [...] | mirbooking-aggregate target mirna
```

Note that the output contains the aggregated columns, the occupancy and the
silencing. Both metrics are summed among the group.

The `mirbooking-extract-cds-regions` tool extract the coding regions from
a GenBank file.

```bash
mirbooking-extract-cds-regions < GCR.gbff > cds-regions.tsv
```

## C API

The API is conform to the GLib style and enable a wide range of use. It is fairly easy to use and a typical
experimentation session is:

 1. create a broker via `mirbooking_broker_new`
 2. create some sequence objects with `mirbooking_target_new` and `mirbooking_mirna_new`
 3. setup quantities via `mirbooking_broker_set_sequence_quantity`
 4. call `mirbooking_broker_run` to perform a full hybridization
 5. retrieve and inspect the microtargetome with `mirbooking_broker_get_target_sites`

For a more detailed usage and code example, the [program source](https://github.com/major-lab/mirbooking/blob/master/bin/mirbooking.c)
is very explicit as it perform a full session and fully output the target sites.

## Parallelization (using GNU parallel)

```bash
parallel mirbooking --mirnas mature.fa --targets hg38.fa --score-table scores ::: wildtype.tsv over-expression.tsv
```

The targets, mirnas and score table files will reuse the same physical memory
across all processes, yielding a minimal memory usage.

