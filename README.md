# Insert_Sias

Fortran programs that build starting configurations of point defects and defect
clusters in bcc and fcc crystals: self-interstitial (SIA) dumbbells, dislocation
loops, C15 clusters, icosahedral clusters, di-interstitials and vacancies.

Originally written by M.-C. Marinica and R. Alexander.

## How it works

Every calculation follows the same three steps.

1. **Reference crystal.** `cube.x` builds a perfect bcc or fcc supercell and
   writes it to `formd.str`.
2. **Defect centres.** `no00.inp` lists the lattice sites where the defects
   sit, as integer coordinates in units of half the lattice parameter. You
   write it by hand for small clusters; `pre_put_vct.x` generates it for
   dislocation loops.
3. **Insertion.** `put_vct.x` asks for the defect type, reads `formd.str` and
   `no00.inp`, and writes the defective crystal in several formats: its own
   text formats (`formd_new.str`, `formd_new_ndm_scl.str`), a LAMMPS data
   file (`formd_new.lmp`) and `xyz` files for visualisation.

The output is a geometric starting point. It is meant to be relaxed afterwards
with an interatomic potential or DFT.

## Quick start

```sh
# build with gfortran or ifx; the executables go to Src/bin
cmake -B build
cmake --build build

# one <110> dumbbell in a 4x4x4 bcc cell
mkdir run && cd run
printf '2\n2.8553\n4 4 4\n' | ../Src/bin/cube.x   # bcc, a0 = 2.8553 A, 4x4x4 cells
printf '1\n4 4 4\n' > no00.inp                    # one defect, at the site (4,4,4)
echo 1 | ../Src/bin/put_vct.x                     # 1 = <110> dumbbell
```

`formd.str` holds the 128-atom reference crystal and `formd_new.str` the
129-atom crystal with the dumbbell. The same two structures are written as
LAMMPS data files, `formd.lmp` and `formd_new.lmp`.

Without CMake, use the makefiles in `Src/`; see [doc/install.md](doc/install.md).

> **Always build the reference crystal with a0 = 2.8553 Å**, for bcc and fcc
> alike and whatever the element. This value is compiled into `put_vct.x`, which
> uses it to locate the lattice sites. Rescale the result to your own lattice
> parameter afterwards; see [doc/concepts.md](doc/concepts.md#the-compiled-lattice-parameter).

## Defect types

The number in the first column is what `put_vct.x` asks for.

| Code | Defect | Atoms per entry of `no00.inp` |
|---|---|---|
| `1` | `<110>` dumbbell | +1 |
| `2` | `<111>` dumbbell | +1 |
| `21` | `<111>` dumbbell with the neighbours along the row spread out (crowdion-like) | +1 |
| `210` | `<111>` dumbbell with one atom shifted by a0/6 [2 -1 -1], for Frank loops in fcc | +1 |
| `3` | `<100>` dumbbell | +1 |
| `31` | `<100>` dumbbell with the neighbours along the row spread out | +1 |
| `22`, `221` | initial and final states of a `<111>` dumbbell jump | +1 |
| `4` | Gao di-interstitial, from a six-atom template | +2 |
| `5` | vacancy | -1 |
| `6` | C15 unit (bcc) | +2, less where units share atoms |
| `67` | C15 cluster, also cut out in a rotated cell | as `6` |
| `7` | icosahedral unit (fcc) | +7, less where units share atoms |
| `77` | icosahedral units plus extra dumbbells and single atoms (fcc) | see the catalogue |

A loop or a cluster is simply a list of many entries of the same type.

## Documentation

| Page | Content |
|---|---|
| [doc/install.md](doc/install.md) | Compiling with CMake or the makefiles, the executables, changing the compiled parameters |
| [doc/concepts.md](doc/concepts.md) | Integer site coordinates, the compiled lattice parameter, using the result in LAMMPS, rescaling and converting it |
| [doc/programs.md](doc/programs.md) | Each executable: what it asks, what it reads and writes |
| [doc/defects.md](doc/defects.md) | Every defect type: geometry, format of `no00.inp`, restrictions |
| [doc/files.md](doc/files.md) | Layout of all input and output files |
| [doc/examples.md](doc/examples.md) | The example directories and copy-paste recipes |
| [doc/troubleshooting.md](doc/troubleshooting.md) | Error messages and known limitations |

## Repository layout

| Path | Content |
|---|---|
| `CMakeLists.txt` | CMake build of the five executables |
| `Src/programs/` | One source file per executable |
| `Src/lib/` | The routines of `put_vct.x`, one file per task |
| `Src/bin/` | The executables, once built |
| `Examples/` | Worked examples and a library of C15 cluster inputs |
| `doc/` | This documentation |

## Background reading

- C15 clusters in bcc metals: M.-C. Marinica, F. Willaime and J.-P. Crocombette,
  Phys. Rev. Lett. **108**, 025501 (2012).
- Loops and C15 clusters across sizes: R. Alexander *et al.*,
  Phys. Rev. B **94**, 024103 (2016).
- Icosahedral (A15) clusters in fcc metals: A. M. Goryaeva *et al.*,
  Nat. Commun. **14**, 3003 (2023).
