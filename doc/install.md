# Installation

## Requirements

- A Fortran compiler: `gfortran` or Intel `ifx`. The classic `ifort` works too.
- CMake 3.16 or newer (3.20 or newer with `ifx`). Without CMake, `make` and
  `ar` are enough for the makefiles.

Nothing else is needed. The few LAPACK routines the code uses are included in
`Src/lib/lapack.f`.

## Compiling with CMake

From the root of the repository:

```sh
cmake -B build
cmake --build build
```

The executables are written to `Src/bin/`, where the examples and this
documentation expect them. The `build/` directory only holds intermediate
files and can be deleted afterwards.

### Choosing the compiler

Left to itself, the build takes `gfortran` if it is installed, otherwise
`ifx`, otherwise `ifort`. The first step prints its choice:

```
-- Insert_Sias: Fortran compiler GNU 16.2.0 (/usr/bin/gfortran)
```

To choose another compiler, name it and use a build directory of its own:

```sh
cmake -B build-ifx -DCMAKE_Fortran_COMPILER=ifx
cmake --build build-ifx
```

Setting the environment variable `FC` has the same effect
(`FC=ifx cmake -B build-ifx`). A build directory keeps the compiler it was
created with, hence the separate directory. Both builds write to `Src/bin/`,
so the last one built is the one you run.

### Options

| Option | Effect |
|---|---|
| `-DCMAKE_Fortran_COMPILER=ifx` | Compiler to use |
| `-DCMAKE_BUILD_TYPE=Debug` | No optimisation, array-bounds checking and a traceback on errors. The default is `Release` |
| `-DINSERT_SIAS_BIN_DIR=/some/dir` | Directory that receives the executables instead of `Src/bin/` |

`cmake --install build --prefix /some/prefix` copies the executables to
`/some/prefix/bin`.

The CMake build was checked with gfortran 14, 15 and 16 and with ifx 2025.
The four compilers give the same structures, digit for digit. A few auxiliary
files (`geom_new.data`, `geom_org.data`, `unit_cell`, `fort.12`) differ
between gfortran and ifx in the layout of some numbers.

## Compiling with the makefiles

The makefiles in `Src/` do the same job without CMake.

```sh
cd Src
cp Makefile_gfortran Makefile    # GNU compiler
make
```

For the Intel compiler use `Makefile_ifx` instead. It calls `ifx`; for the
classic compiler set `FC = ifort` in it, or type `make FC=ifort`.

The executables are written to `Src/bin/`, the module files to `Src/Modules/`
and the object files to `Src/obj/`. The first two directories are part of the
repository, so they exist after cloning; `make` creates the third.

Both makefiles were checked, with gfortran 16 and ifx 2025. They give the same
results as the CMake build.

## The executables

| Executable | Source | Purpose |
|---|---|---|
| `cube.x` | `Src/programs/cube.f` | Builds the perfect bcc or fcc reference crystal |
| `put_vct.x` | `Src/programs/put_vct.f90` and `Src/lib/` | Inserts the defects |
| `pre_put_vct.x` | `Src/programs/pre_put_vct.f90` | Generates the list of sites of a dislocation loop |
| `read_write_n00.x` | `Src/programs/read_write_n00.f90` | Translates a list of sites to another place in the box |
| `cube_shell.x` | `Src/programs/cube_shell.f90` | Builds a reference crystal with the atoms of an inner block listed first |

Most work needs only `cube.x` and `put_vct.x`, plus `pre_put_vct.x` for loops.
See [programs.md](programs.md) for what each one asks, and `Src/README` for a
map of the source files.

There is no install step. Call the executables by their path, as the examples
do (`../../Src/bin/put_vct.x`), or add `Src/bin` to your `PATH`.

## Recompiling after a change

Run the build command again, `cmake --build build` or `make`. Both recompile
what the change affects.

`make clean` removes what `make` created: the object files, the module files,
the library and the executables.

## Compiled parameters

Several physical parameters are constants in the source. To change one, edit
the file and recompile.

| Parameter | Where | Default | Meaning |
|---|---|---|---|
| `a0` | `Src/programs/put_vct.f90`, top of the program | 2.8553 Å | Lattice parameter of the reference crystal |
| `as110`, `as111`, `as100` | the same line | 2.05, 1.90, 1.70 Å | Distance between the two atoms of a `<110>`, `<111>`, `<100>` dumbbell |
| `a0`, `as110`, `as111`, `as100` | `Src/programs/pre_put_vct.f90`, first lines of code | 2.8553, 1.6, 1.6, 1.6 Å | The same for `pre_put_vct.x` |
| `deca` | `Src/programs/pre_put_vct.f90`, bi-plane block | (1, 0, 1) | Shift between the two planes of a bi-plane loop, in units of a0/2 |
| `pdist` | `Src/lib/kazuto_style.f90`, routine `write_out_kazuto` | 25 | Size of the cut for defect type `67` |

The dumbbell lengths in `put_vct.f90` suit bcc. The line just below them in the
source, commented out and marked `!fcc`, holds the fcc values (1.55 Å for all
three). Swap the two lines before building loops in fcc: with the bcc values,
neighbouring `<110>` dumbbells of an fcc loop end up less than 1 Å apart.

Changing `a0` is rarely needed. The usual route is to keep 2.8553 Å and rescale
the result, as described in [concepts.md](concepts.md#the-compiled-lattice-parameter).
