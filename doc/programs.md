# The programs

All programs are interactive: they print a question and read the answer from
standard input. They work in the current directory and use fixed file names.
The usual way to run them is to pipe the answers in, as the examples do:

```sh
printf '2\n2.8553\n8 8 8\n' | cube.x
echo 6 | put_vct.x
```

Or keep the answers in a small file (`cube.inp`, `put.inp`) and redirect it:
`cube.x < cube.inp`.

Because the file names are fixed, run each case in its own directory. A second
run overwrites the results of the first.

## cube.x: the reference crystal

| Question | Answer |
|---|---|
| `FCC=1/BCC=2` | `1` for fcc, `2` for bcc |
| `Put the size of the cube Cu 3.615/ Fe 2.8553:` | The lattice parameter in Å. Use **2.8553** ([why](concepts.md#the-compiled-lattice-parameter)) |
| `n1,n2,n3?` | Number of cubic unit cells along x, y and z |

It writes `formd.str` (the crystal, read by the other programs), together with
`unit_cell`, `unit_vect` and `fort.11`, which hold the same crystal in other
layouts.

The cell holds 2 x n1 x n2 x n3 atoms for bcc and 4 x n1 x n2 x n3 for fcc.

## put_vct.x: inserting the defects

| Question | Answer |
|---|---|
| `Choice the sia (1 - 110)(2 - 111) ...` | The code of the defect type |

The codes are listed in the [defect catalogue](defects.md).

| Reads | When |
|---|---|
| `formd.str` | always |
| `no00.inp` | all types except `22` and `221` |
| `no00m.inp` | types `22` and `221` |
| `gao_unit.in` | type `4` |

It writes the defective crystal and a set of visualisation files, described in
[files.md](files.md). The file to carry on with is `formd_new.lmp` for LAMMPS,
`formd_new.str` (Cartesian) or `formd_new_ndm_scl.str` (reduced coordinates).
The reference crystal is written to `formd.lmp` as well, so that both
structures are at hand in the same format.

On the screen it reports the number of atoms of the reference crystal, the
number of entries read, the largest site indices found (`n1MAX,n2MAX,n3MAX`),
the LAMMPS files written and, for the cluster types, the net number of
interstitials created:

```
 The no of defects SIAs ...................:          11
LAMMPS data file written ..........: formd_new.lmp (16011 atoms)
```

## pre_put_vct.x: the sites of a dislocation loop

This program cuts a planar patch of sites out of the reference crystal and
writes it as a list for `put_vct.x`. It reads `formd.str`.

| Question | Answer |
|---|---|
| `1/plane, 2/bi-plane` | `1` for a loop one atomic plane thick, `2` for two adjacent planes |
| `The habit plane ... 1 - 100/ 2 - 111 / 3 - 110` | The plane the loop lies in |
| `The shape of the loop ... (1/rectangle, 2/circle, 3/hexa)` | The outline |
| size | See below |
| `Choice the sia ... (1/<110>/2/<111>/3<100>` | Dumbbell type for the program's own output |

The size questions depend on the outline. All lengths are in Å.

| Outline | Questions | Sites kept |
|---|---|---|
| rectangle | two sides | within half a side of the centre along each of two in-plane directions |
| circle | radius | closer to the centre than the radius |
| hexa | one size (asked as "radius") | within half the size of the centre along each edge normal |

Details of the outlines:

- The rectangle is aligned with [100] and [010] in the (100) plane, and with
  [1 -1 0] and [1 1 -2] in the (111) plane. In the (110) plane the program asks
  for the side along z first and the other side second, and both edges come out
  √2 times longer than the numbers entered.
- "hexa" gives a hexagon in the (111) and (110) planes and an octagon in the
  (100) plane.

The plane passes through the atom nearest to the middle of the box, and the
outline is centred on that atom. The program prints how many sites it found:

```
 In this plane the number of atoms....:         430
 The size of the loop (no of atoms)...:          55
```

The second number is the number of interstitials of the loop (here a circle of
radius 8 Å in the (111) plane of a 12 x 12 x 12 fcc cell). Adjust the size
until it is the one you want.

The result is `no00_new.inp`. To build the loop, copy it to `no00.inp` and run
`put_vct.x` with the dumbbell type of your choice:

```sh
cp no00_new.inp no00.inp
echo 2 | put_vct.x
```

Things to know:

- `pre_put_vct.x` also builds the loop itself and writes `formd_new.str`,
  `formd_new_ndm_scl.str`, `geom_new.data` and `sia.xyz`, with its own dumbbell
  length and a cell expanded to keep the volume per atom constant. `put_vct.x`
  overwrites these files, and the two-step route above is the one the examples
  use.
- The limits it prints before the size questions ("must be smaller than") are
  indications only. For the hexagon in the (111) plane the printed limit is 0
  and should be ignored.
- The bi-plane option works for fcc only; see [troubleshooting.md](troubleshooting.md#known-limitations).
- `fort.12` receives the sites of the cut plane and of the loop, for plotting.

## read_write_n00.x: moving a list of sites

Site lists are tied to a place in the box. This program shifts a list by a
constant integer vector, for instance to move a cluster designed in a small
cell to the middle of a large one. It reads `no00.inp` and writes
`no00_new.inp`.

It first prints the largest, smallest and middle value of each index, then
asks:

| Question | Answer |
|---|---|
| `Type of no00 file ? (1 - traditional loops, 2 - C15 clusters, 3 - migration with care)` | `1` for a plain list of sites, `2` for a C15 list (the unit types that follow the sites are copied). `3` is not implemented |
| `The step on X?` | Shift of n1 |
| `The step on Y?` | Shift of n2 |
| `The step on Z?` | Shift of n3 |

The shift must be a lattice vector: three even or three odd numbers in bcc, an
even sum in fcc.

## cube_shell.x: a reference crystal in two blocks

A variant of `cube.x` for work that embeds a small cell in a large one. It
writes the atoms of an inner block first and the rest of the crystal after
them.

| Question | Answer |
|---|---|
| `FCC=1/BCC=2/TT=3` | `1` fcc, `2` bcc, `3` the 24-atom cubic cell of the C15 Laves structure |
| `Put the size of the cube ...` | Lattice parameter in Å |
| `the outter shell n1,n2,n3?` | Number of unit cells of the whole crystal |
| `the inner shell n1, n2,n3?` | Number of unit cells of the inner block |

The inner block is the corner of the crystal that starts at the origin.

| File | Content |
|---|---|
| `formd.str` | The whole crystal, inner block first |
| `formd.str_inner` | The inner block alone, with its own cell |
| `formd.gin` | The whole crystal in reduced coordinates |
