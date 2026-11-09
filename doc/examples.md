# Examples

The first part of this page gives short recipes that can be pasted into a
terminal. The second part describes the directories under `Examples/`.

The recipes assume an empty working directory two levels below the repository
root, like the example directories, so that the executables are at
`../../Src/bin/`. Each one was run to produce the atom counts quoted.

## Recipes

### A single dumbbell in bcc

```sh
printf '2\n2.8553\n4 4 4\n' | ../../Src/bin/cube.x
printf '1\n4 4 4\n' > no00.inp
echo 1 | ../../Src/bin/put_vct.x
```

128 atoms become 129. Replace the last `1` by `2` or `3` for a `<111>` or
`<100>` dumbbell.

As in every recipe below, the result is in `formd_new.str` and, for LAMMPS, in
`formd_new.lmp`; the perfect crystal is in `formd.str` and `formd.lmp`.

### A dislocation loop in bcc

A circular loop of `<111>` dumbbells in a (111) plane, in a 16 x 16 x 16 cell.
The answers to `pre_put_vct.x` are, in order: single plane, (111) habit plane,
circle, radius 12 Å, and the dumbbell type for its own output.

```sh
printf '2\n2.8553\n16 16 16\n' | ../../Src/bin/cube.x
printf '1\n2\n2\n12.0\n2\n'    | ../../Src/bin/pre_put_vct.x
cp no00_new.inp no00.inp
echo 2 | ../../Src/bin/put_vct.x
```

The loop holds 31 dumbbells: 8192 atoms become 8223.

### A Frank loop in fcc

The same procedure in fcc, finished with type `210`.

```sh
printf '1\n2.8553\n12 12 12\n' | ../../Src/bin/cube.x
printf '1\n2\n2\n8.0\n2\n'     | ../../Src/bin/pre_put_vct.x
cp no00_new.inp no00.inp
echo 210 | ../../Src/bin/put_vct.x
```

55 sites: 6912 atoms become 6967. Finishing with `1` instead of `210` gives a
loop of `<110>` dumbbells on the same sites; compile the fcc dumbbell lengths
first ([install.md](install.md#compiled-parameters)).

### A vacancy cluster

A di-vacancy on first neighbours in a 6 x 6 x 6 bcc cell.

```sh
printf '2\n2.8553\n6 6 6\n' | ../../Src/bin/cube.x
printf '2\n4 4 4\n5 5 5\n'  > no00.inp
echo 5 | ../../Src/bin/put_vct.x
```

432 atoms become 430.

### A C15 cluster from the library

The 11-interstitial C15 cluster in a 20 x 20 x 20 bcc cell.

```sh
printf '2\n2.8553\n20 20 20\n' | ../../Src/bin/cube.x
cp ../../Examples/C15_inputs/MagicNumbers/0011_c15a.inp no00.inp
echo 6 | ../../Src/bin/put_vct.x
```

16000 atoms become 16011, and the screen shows
`The no of defects SIAs ...: 11`.

### Moving a cluster

To shift the cluster of the previous recipe by (10, 10, 10), for a larger
cell, answer `2` (C15 list) and the three steps to `read_write_n00.x`:

```sh
printf '2\n10\n10\n10\n' | ../../Src/bin/read_write_n00.x
cp no00_new.inp no00.inp
```

### An icosahedral cluster in fcc

Two linked icosahedral units in a 10 x 10 x 10 fcc cell.

```sh
printf '1\n2.8553\n10 10 10\n' | ../../Src/bin/cube.x
printf '2\n3 3 3\n5 3 3\n'     > no00.inp
echo 7 | ../../Src/bin/put_vct.x
```

4000 atoms become 4013.

### A dumbbell jump

The two end states of one `<111>` dumbbell jumping next to a second one that
stays, in a 6 x 6 x 6 bcc cell.

```sh
printf '2\n2.8553\n6 6 6\n'          | ../../Src/bin/cube.x
printf '2 1 1\n4 4 4 1\n6 4 4 0\n'   > no00m.inp
echo 22 | ../../Src/bin/put_vct.x
```

The states are in `deb_loop.gin_0` and `fin_loop.gin_0`, and for LAMMPS in
`deb_loop.lmp` and `fin_loop.lmp`, 434 atoms each.

## The example directories

### Examples/001_SIA_110

The introductory example: one dumbbell in a bcc cell, driven by a script.

| File | Role |
|---|---|
| `insert.sh` | Writes `cube.inp`, `no00.inp` and `put.inp`, then runs `cube.x` and `put_vct.x` |
| `clean.sh` | Removes the files produced by a run |
| `README` | Two-line instructions |
| other files | The result converted to other formats (POSCAR, xyz, LAMMPS, CFG, gin) with external tools, kept as illustrations |

Run it with `sh insert.sh`. The script contains three small blocks to edit:
the cell (lattice type, lattice parameter, number of cells), the list of sites,
and the defect type. Keep the lattice parameter at 2.8553.

### Examples/004a_SIA_C15, 011a_SIA_C15, 017a_SIA_C15, 101_SIA_c15

Complete runs for C15 clusters of 4, 11, 17 and 101 interstitials. Each
directory holds the reference crystal (`formd.str`), the list of units
(`no00.inp`) and every output file.

| Directory | Reference cell | Atoms after insertion |
|---|---|---|
| `004a_SIA_C15` | 4 x 4 x 4 bcc, 128 atoms | 132 |
| `011a_SIA_C15` | 20 x 20 x 20 bcc, 16000 atoms | 16011 |
| `017a_SIA_C15` | 20 x 20 x 20 bcc, 16000 atoms | 16017 |
| `101_SIA_c15` | 20 x 20 x 20 bcc, 16000 atoms | 16101 |

To reproduce one, run `echo 6 | ../../Src/bin/put_vct.x` in its directory. The
current sources regenerate the stored `formd_new.str` exactly, which makes
these directories a convenient check of a new build. The run adds the LAMMPS
files `formd.lmp` and `formd_new.lmp`, which the stored results predate.

### The C15 library

Ready-made unit lists for C15 clusters, to be copied to `no00.inp` and used
with type `6`.

| Location | Files |
|---|---|
| `Examples/C15_inputs/` | `0004_c15a.inp` and `0101_c15a.inp`, the clusters of 4 and 101 interstitials |
| `Examples/C15_inputs/MagicNumbers/` | Clusters of 11, 17, 23, 29, 35, 37 and 47 interstitials |
| `Examples/C15_inputs/Draft/` | About forty clusters of 30 to 90 interstitials, several variants per size |

The file names read `NNNN_c15x.inp`: NNNN is the number of interstitials and
the letter distinguishes different shapes of the same size.

The centres in these files lie between 5 and 23, which places the clusters
near the middle of a 20 x 20 x 20 bcc cell. For a larger cell, shift them with
`read_write_n00.x`.

### Examples/fcc_big_loops

A loop in fcc. The directory holds a 12 x 12 x 12 fcc reference crystal
(`formd.str`, 6912 atoms, built with a0 = 2.8553 Å), a list of 19 sites in a
(111) plane (`no00.inp`), and a `README` with the three steps: `cube.x`,
`pre_put_vct.x`, then `put_vct.x` after copying `no00_new.inp` to `no00.inp`.

Running `put_vct.x` there with type `1`, `2` or `210` gives a 19-interstitial
loop, 6931 atoms.
