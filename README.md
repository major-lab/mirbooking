# miRBooking

Implementation of the miRBooking algorithm and metrics in C

 - fast and memory efficient
 - usable from Python, JavaScript and Vala via GObject introspection
 - usable from Java via JNA
 - memory-mapped score tables, target and mirnas (for parallel execution)
 - gzip-compressed score tables
 - memory-mapped and zero-copy for input files containing sequences (i.e. FASTA)
 - binary with support for static linking for more portability
 - stdin/stdout for piping from and into other tools

## Usage

```bash
mirbooking --mirnas mature.fa
           --targets GCF_000001405.37_GRCh38.p11_rna.fna
           --score-table scores|scores.gz
           [--cds-regions cds-regions.tsv]
           [--kappa 6.4135955621023e-5]
           [--step-size 5e-6]
           [--tolerance 1e-3]
           [--max-iterations 1e8]
           [--5prime-footprint 26]
           [--3prime-footprint 19]
           [--5prime-utr-multiplier 0.1]
           [--cds-multiplier 0.1]
           [--3prime-utr-multiplier 1.0]
           [--quantities quantities.tsv]
           [--output output.tsv]
```

To obtain detailed usage and options, launch `mirbooking --help`.

The command line program expects a number of inputs:

 - `--targets`, a FASTA containing mRNA transcripts where the identifier is the
   accession (i.e. NM_002710.3)
 - `--mirnas`, a FASTA containing mature miRNAs where the first token in the
   comment is the accession (i.e. MIMAT0004792)
 - `--score-table`, a score table with free energy for seeds
 - `--quantities`, a quantity file mapping target and mirna accessions to
   expressed quantity in FPKM/RPKM/RPM

The `--kappa` parameter indicates how many nM (1e-9 Molar) of concentration
a FPKM represents in the quantification. It can be estimated using spike-ins if
their prior concentration is known. The default value has been calculated from
ENCODE's HeLa S3 reference epigenome[^hela-s3-encode] using the following
procedure:

 1. Obtain spike-ins concentration from Thermofisher and convert attomoles/µL
    to nM
 2. Take ~2% of that total concentration since this is the mixture
 3. Compute the total FPKM for the spike-ins for both replicate
 4. Divide the concentration in nM by the mean FPKM from both replicates to
    obtain $\kappa$

[^hela-s3-encode]: https://www.encodeproject.org/reference-epigenomes/ENCSR068MRQ/
[^spike-ins-thermofisher]: https://www.thermofisher.com/order/catalog/product/4456740

For control conditions (i.e. wildtype or empty vector), long RNA-Seq (in RPKM
or FPKM) with small RNA-Seq (RPM) quantifications should be combined as-is. If
you have spike-ins, use them to adjust the fragment or read counts.

For simulated over-expression conditions, an arbitrary value (i.e 1e5) can be
added to endogenous or as a synthetic microRNA. No renormalization should be
performed: the resulting library size will just be bigger, as expected for such
an experiment.

For knock-out conditions, zero the miRNA or mRNA expression.

To compute the silencing, `--cds-regions` is a two columns TSV document mapping
target accessions to coding regions the 'a..b' format where a and b are
inclusive 1-based indexes

The output is a TSV with the following columns:

| Column           | Description                                            |
| ---------------- | ------------------------------------------------------ |
| target_accession | Target accession with version                          |
| target_name      | Name of the target or N/A if unknown                   |
| target_quantity  | Number of targets                                      |
| target_silencing | Proportion targets silenced by the miRNAs              |
| position         | Site position on the target                            |
| location         | Location (e.g. 5', CDS, 3' or N/A if unknown)          |
| occupancy        | Proportion of occupied sites at this position          |
| mirna_accession  | miRNA accession                                        |
| mirna_name       | Name of the miRNA or N/A if unknown                    |
| mirna_quantity   | Number of miRNAs                                       |
| score            | Molar Gibbs free energy of the miRNA::MRE duplex       |
| quantity         | Number of occupants miRNAs at a this target position   |

## Installation

You'll need [Meson](http://mesonbuild.com/) and [Ninja](http://ninja-build.org/)
as well as GLib development files installed on your system.

```bash
mkdir build && cd build
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

FFTW can be optionally used to compute more accurate silencing by specifying
`meson -Dwith_fftw3=true`. If you redistribute miRBooking source code, be
careful not to enable this as a default because of the GPL license covering
this dependency.

OpenMP can be optionally used to parallelize the evaluation of partial
derivatives by specifying `-Dwith_openmp=true`.

## Other tools

In addition to the `mirbooking` binary, this package ship a number of utilities
to process quantity files and compute summaries based on Pandas.

Te `mirbooking-generate-score-table` compute a hybridization energy table for
the given seed size and upper bound on the number of mismatches.
[MC-Flashfold](https://major.iric.ca/mc-tools/) is required.

```bash
mirbooking-generate-score-table [--mcff=mcff]
                                [--mcff-args]
                                [--seed-length=7]
                                [--max-mismatches=1]
                                [--max-workers=1]
                                --output scores
```

The `mirbooking-calibrate` tool is expecting a transcript and miRNA
quantification (e.g. two-column TSV document mapping accession to quantity) and
process it such that it contains approximately the same amount of each kind by
rescaling the smallest one toward the biggest one. It emits a calibrated output
suitable for `mirbooking`.

```bash
mirbooking-calibrate targets.tsv mirnas.tsv | mirbooking [...]
```

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

