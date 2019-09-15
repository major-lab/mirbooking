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
           [--supplementary-model none]
           [--supplementary-scores scores-3mer]
           [--input stdin]
           [--output stdout]
           [--output-format tsv]
           [--sparse-solver best-available]
           [--max-iterations 100]
           [--5prime-footprint 9]
           [--3prime-footprint 7]
           [--cutoff 100]
           [--relative-cutoff 0]
```

To obtain detailed usage and options, launch `mirbooking --help`.

The command line program expects a number of inputs:

 - `--targets`, a FASTA containing RNA transcripts where the identifier is the
   accession (i.e. NM_002710.3)
 - `--mirnas`, a FASTA containing mature miRNAs where the first token in the
   comment is the accession (i.e. MIMAT0004792)
 - `--seed-scores`, a sparse score table of seed free energies which can be
   generated using `generate-score-table` program described below
 - `--accessibilitiy-scores` contains entries with position-wise free energy
   contribution (or penalty) on the targets
 - `--supplementary-scores` contains either 4mer or 3mer
 - `--input`, a quantity file mapping target and mirna accessions to
   expressed quantity in picomolars

Tables for seed and supplementary scores are provided in the `data` folder.
These were computed with RNAcofold binding energy from ViennaRNA package.

Note that Yan et al. (`--supplementary-model=yan-et-al-2018`) model requires a
3mer table whereas Zamore et al. (`--supplementary-model=zamore-et-al-2012`)
require a 4mer table.

Tables for seed and supplementary bindings are automatically located (new in 2.3).

The `--cutoff` parameter can exploit a known upper bound on the complex
concentration to adjust the granularity of the model. Only interaction that can
ideally reach the specified picomolar concentration will be modeled.

The `--relative-cutoff` parameter is similar, but instead filter based on the
ideal substrate bound fraction.

The output is a TSV with the following columns:

| Column           | Description                                            |
| ---------------- | ------------------------------------------------------ |
| gene_accession   | Gene accession with version (new in 2.3)               |
| gene_name        | Name of the gene or N/A if unknown (new in 2.3)        |
| target_accession | Target accession with version                          |
| target_name      | Name of the target or N/A if unknown                   |
| target_quantity  | Total target concentration in picomolars               |
| position         | Site position on the target                            |
| mirna_accession  | miRNA accession                                        |
| mirna_name       | Name of the miRNA or N/A if unknown                    |
| mirna_quantity   | Total miRNA concentration in picomolars                |
| score            | Michalis-Menten constant of the miRNA::MRE duplex      |
| quantity         | miRNA::MRE duplex concentration this target position   |

The detailed TSV output which expands the score structure in its constituents
can be used with `--output-format=tsv-detailed` (new in 2.3). In this mode, the
`score` column is replaced by `kf`, `kr`, `kcleave`, `krelease`, `kcat`,
`kother`, `kd` and `km`.

The GFF3 output can be used with `--output-format=gff3`. The score will
indicate the bound fraction of the position.

Wiggle output can also be produced with `--output-format=wig`. The score will
be the position-wise bound fraction of substrate which properly account for
overlapping microRNA.

## Installation

You'll need [Meson](http://mesonbuild.com/) and [Ninja](http://ninja-build.org/)
as well as GLib development files installed on your system.

```bash
mkdir build && cd build
meson --buildtype=release
ninja
ninja install
```

To generate fast code, configure with `meson -Doptimization=3`.

You can perform a local installation using `meson --prefix=$HOME/.local`, but
you'll need `LD_LIBRARY_PATH` set accordingly since the `mirbooking` program
uses a shared library. Otherwise, a static linkage can be done by calling
`meson --default-library=static`.

To generate introspection metadata, use `meson -Dwith_introspection=true`. To
generate Vala bindings, use `meson -Dwith_vapi=true`.

FFTW can be optionally used to compute more accurate silencing by specifying
`meson -Dwith_fftw3=true`. If you redistribute miRBooking source code, be
careful not to enable this as a default because of the GPL license covering
this dependency. If you have access to Intel MKL, you can alternatively use its
FFTW3 implementation with `-Dwith_mkl_fftw3=true`.

OpenMP can be optionally used to parallelize the evaluation of partial
derivatives and some supported solvers by specifying `-Dwith_openmp=true`.

MPI can be optionally used to distribute the computation across multiple
machine on supported solvers (i.e. `mkl-cluster`) by specifying `-Dwith_mpi=true`.

| Solver      | Build Options                                                                      |
| ----------- | ---------------------------------------------------------------------------------- |
| SuperLU     | `-Dwith_superlu=true`                                                              |
| SuperLU MT  | `-Dwith_superlu_mt=true`                                                           |
| UMFPACK     | `-Dwith_umfpack=true`                                                              |
| cuSOLVER    | `-Dwith_cuda=<cuda_toolkit_api_version> -Dwith_cusolver=true`                      |
| MKL DSS     | `-Dwith_mkl=true -Dmkl_root=<path to mkl> -Dwith_mkl_dss=true`                     |
| MKL Cluster | `-Dwith_mpi=true -Dwith_mkl=true -Dmkl_root=<path to mkl> -Dwith_mkl_cluster=true` |
| MKL LAPACK  | `-Dwith_mkl=true -Dmkl_root=<path to mkl> -Dwith_mkl_lapack=true`                  |
| PARDISO     | `-Dwith_pardiso=true`                                                              |

cuSOLVER require CUDA toolkit whose API version is to be specified with
`-Dwith_cuda=<cuda_toolkit_api_version>`.

MKL DSS and MKL Cluster can benefit from [TBB](https://www.threadingbuildingblocks.org/)
instead of OpenMP, which can be enabled with `-Dwith_mkl_tbb=true`.

MKL LAPACK is not a sparse linear solver and thus will not handle typical
workload very well, but it will perform orders of magnitude faster on dense
jacobians.

MKL DSS and MKL Cluster can be used with the 64 bit interface, allowing much
larger systems to be solved with `-Dwith_mkl_ilp64=true`. However, this will
break other solvers as it will load a 64 bit BLAS.

PARDISO cannot be used along with MKL DSS because they define common symbols.

By default, the best sparse solver available among the following will be used
(new in 2.3):

 1. MKL-DSS
 2. PARDISO
 2. UMFPACK
 3. SuperLU (even if not available)

## Numerical integration

In addition to determine the steady state, miRBooking can also perform
numerical integration of the microtargetome using the programming API.

## Other tools

In addition to the `mirbooking` binary, this package ship a number of
utilities.

Te `generate-score-table` compute a hybridization energy table for a given seed
mask. Either [ViennaRNA](https://www.tbi.univie.ac.at/RNA/) or
[mcff](https://major.iric.ca/mc-tools) is required to compute energies.

```bash
generate-score-table [--method=RNAcofold]
                     [--temperature=310.5]
                     [--mask=||||...]
                     [--hard-mask=||||...]
                      --output scores
```

The seed mask defines folding constraints on the target with `|` for
a canonical match, `x` for a canonical mismatch and `.` for no constraint. It
also determines the seed length. If a hard mask is provided, unsatisfying
interactions are filtered out (new in 2.3).

It's also possible to ajust the folding temperature (new in 2.3).

The number of workers can be tuned by setting `OMP_NUM_THREADS` environment
variable.

## C API

The API is conform to the GLib style and enable a wide range of use. It is fairly easy to use and a typical
experimentation session is:

 1. create a broker via `mirbooking_broker_new`
 2. create some sequence objects with `mirbooking_target_new` and `mirbooking_mirna_new`
 3. setup quantities via `mirbooking_broker_set_sequence_quantity`
 4. call `mirbooking_broker_evaluate` and `mirbooking_broker_step` repeatedly
    to perform a full hybridization or numerical integration
 5. retrieve and inspect the microtargetome with `mirbooking_broker_get_target_sites`

For a more detailed usage and code example, the main program source in
`bin/mirbooking.c` is very explicit as it perform a full session and fully
output the target sites.
