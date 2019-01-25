# miRBooking

Implementation of the miRBooking algorithm and metrics in C

 - fast and memory efficient
 - usable from Python, JavaScript and Vala via GObject introspection
 - memory-mapped score tables, target and mirnas FASTA for  low memory
   footprint in parallel execution
 - binary with support for static linking for more portability
 - stdin/stdout for piping from and into other tools

## Usage

```bash
mirbooking --targets GCF_000001405.37_GRCh38.p11_rna.fna
           --mirnas mature.fa
           --seed-scores scores-7mer-3mismatch-ending
           [--accessibility-scores accessibility-scores[.gz]]
           [--supplementary-scores scores-3mer]
           [--input stdin]
           [--output stdout]
           [--output-format tsv|gff3]
           [--sparse-solver superlu]
           [--tolerance 1e-8]
           [--max-iterations 100]
           [--5prime-footprint 9]
           [--3prime-footprint 7]
           [--cutoff 100]
           [--verbose]
```

To obtain detailed usage and options, launch `mirbooking --help`.

The command line program expects a number of inputs:

 - `--targets`, a FASTA containing RNA transcripts where the identifier is the
   accession (i.e. NM_002710.3)
 - `--mirnas`, a FASTA containing mature miRNAs where the first token in the
   comment is the accession (i.e. MIMAT0004792)
 - `--seed-scores`, a sparse score table of seed free energies which can be
   generated using `mirbooking-generate-score-table` program described below
 - `--accessibilitiy-scores` contains entries with position-wise free energy
   contribution (or penalty) on the targets
 - `--supplementary-scores` contains either 4mer or 3mer
 - `--input`, a quantity file mapping target and mirna accessions to
   expressed quantity in pM

The `--cutoff` parameter can exploit a known upper bound on the complex
concentration to adjust the granularity of the model. Only interaction that can
ideally reach the specified picomolar concentration will be modeled.

The output is a TSV with the following columns:

| Column           | Description                                            |
| ---------------- | ------------------------------------------------------ |
| target_accession | Target accession with version                          |
| target_name      | Name of the target or N/A if unknown                   |
| target_quantity  | Total target concentration in picomolars               |
| position         | Site position on the target                            |
| mirna_accession  | miRNA accession                                        |
| mirna_name       | Name of the miRNA or N/A if unknown                    |
| mirna_quantity   | Total miRNA concentration in picomolars                |
| score            | Michalis-Menten constant of the miRNA::MRE duplex      |
| quantity         | miRNA::MRE duplex concentration this target position   |

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

MPI can be optionally used to distribute the computation across multiple
machine on supported solvers (i.e. `mkl-cluster`) by specifying `-Dwith_mpi=true`.

## Numerical integration

In addition to determine the steady state, miRBooking can also perform
numerical integration of the microtargetome using the programming API.

## Other tools

In addition to the `mirbooking` binary, this package ship a number of
utilities.

Te `mirbooking-generate-score-table` compute a hybridization energy table for
the given seed size and upper bound on the number of mismatches.
[MC-Flashfold](https://major.iric.ca/mc-tools/) is required.

The number of workers can be tuned by setting `OMP_NUM_THREADS` environment
variable.

```bash
mirbooking-generate-score-table [--mcff=mcff]
                                [--mask=....xxx]
                                --output scores
```

To use `RNAcofold` instead from ViennaRNA package, a `mcff-ViennaRNA` script is
provided to perform the conversion:

```bash
mirbooking-generate-score-table --mcff=scripts/mcff-ViennaRNA
                                [--mask=....xxx]
                                --output scores
```

## C API

The API is conform to the GLib style and enable a wide range of use. It is fairly easy to use and a typical
experimentation session is:

 1. create a broker via `mirbooking_broker_new`
 2. create some sequence objects with `mirbooking_target_new` and `mirbooking_mirna_new`
 3. setup quantities via `mirbooking_broker_set_sequence_quantity`
 4. call `mirbooking_broker_evaluate` and `mirbooking_broker_step` repeatedly
    to perform a full hybridization
 5. retrieve and inspect the microtargetome with `mirbooking_broker_get_target_sites`

For a more detailed usage and code example, the [program source](https://github.com/major-lab/mirbooking/blob/master/bin/mirbooking.c)
is very explicit as it perform a full session and fully output the target sites.
