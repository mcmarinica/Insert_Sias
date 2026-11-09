# File formats

All files are plain text and have fixed names. Lengths are in Å.

## Input files

### formd.str: the reference crystal

Written by `cube.x`, read by `put_vct.x` and `pre_put_vct.x`.

```
     128
  11.4212000000000   0.0000000000000   0.0000000000000
   0.0000000000000  11.4212000000000   0.0000000000000
   0.0000000000000   0.0000000000000  11.4212000000000
Fe    0.0000000000000   0.0000000000000   0.0000000000000
Fe    1.4276500000000   1.4276500000000   1.4276500000000
...
```

| Lines | Content | Fortran format |
|---|---|---|
| 1 | number of atoms | `i8` |
| 2-4 | the three cell vectors, one per line | `3f18.13` |
| then | one line per atom: label and Cartesian coordinates | `a3,3f18.13` |

The file is read with these **fixed formats**: the label occupies the first
three columns and each number a field of 18 columns. A file written by hand or
by another tool must respect the columns exactly.

### no00.inp: the defect sites

Read in free format. The first line holds the number of entries; the sites
follow, one per line, as three integers in units of a0/2. Some defect types
append further lines.

| Defect type | Layout |
|---|---|
| `1`, `2`, `21`, `210`, `3`, `31`, `4`, `5`, `7` | `N`, then N lines `n1 n2 n3` |
| `6`, `67` | `N`, then N lines `n1 n2 n3`, then N lines with the unit type (`1`, `2` or `4`) |
| `77` | `Nico Nextra`, then Nico lines `n1 n2 n3`, then Nextra lines `n1 n2 n3 t` |

Details and examples are in the [defect catalogue](defects.md).

### no00m.inp: dumbbell jumps

Used by types `22` and `221` in place of `no00.inp`: a first line
`N Nmig Nfix`, then N lines `n1 n2 n3 m` with m = 1 for a dumbbell that jumps
and 0 for one that stays.

### gao_unit.in: template for type 4

Six lines `x y z`, the positions of the six inserted atoms relative to the
listed site, in units of a0.

## Output of put_vct.x

### The defective crystal

These four files hold the same structure. The inserted atoms come first, the
remaining atoms of the crystal after them.

| File | Layout |
|---|---|
| `formd_new.str` | As `formd.str`: number of atoms, cell, `Fe x y z` in Cartesian coordinates |
| `formd_new_ndm_scl.str` | The `.gin` layout below, in reduced coordinates |
| `formd_new.lmp` | LAMMPS data file, described below |
| `geom_new.data` | Number of atoms, a line `1.0`, cell, then `x y z` in Cartesian coordinates |

The `.gin` layout, also used by the `.gin` files of the examples:

```
1   1   1
  11.4212000000000   0.0000000000000   0.0000000000000
   0.0000000000000  11.4212000000000   0.0000000000000
   0.0000000000000   0.0000000000000  11.4212000000000
     129
      0.4365404291391   0.4365404291391   0.5000000000000   1
      0.5634595708609   0.5634595708609   0.5000000000000   1
...
```

| Lines | Content |
|---|---|
| 1 | `1 1 1` |
| 2-4 | the three cell vectors |
| 5 | number of atoms |
| then | reduced coordinates and a species index, always `1` |

### LAMMPS data files

```
# Insert_Sias: crystal with defects (formd_new.str), lengths in Angstrom

129 atoms
1 atom types

      0.0000000000000000     11.4212000000000007   xlo xhi
      0.0000000000000000     11.4212000000000007   ylo yhi
      0.0000000000000000     11.4212000000000007   zlo zhi

Atoms # atomic

         1     1      4.9858155492837888      4.9858155492837888      5.7106000000000003
         2     1      6.4353844507162119      6.4353844507162119      5.7106000000000003
...
```

| File | Content | Written for |
|---|---|---|
| `formd.lmp` | the reference crystal, as `formd.str` | all types |
| `formd_new.lmp` | the crystal with the defects, as `formd_new.str` | all types |
| `deb_loop.lmp`, `fin_loop.lmp` | initial and final states of a jump, as the `.gin_0` files | types `22`, `221` |

The files use atom style `atomic` and one atom type; each atom line holds the
atom number, the type and the Cartesian coordinates in Å. The atoms are in the
order of the file named in the second column, with full double precision. The
box starts at the origin. A line `xy xz yz` is added only for a tilted cell.

How to read and rescale them in LAMMPS is shown in
[concepts.md](concepts.md#using-the-result-in-lammps); for other formats see
[concepts.md](concepts.md#converting-to-other-formats).

### The reference crystal again

| File | Layout | Written for |
|---|---|---|
| `formd.lmp` | LAMMPS data file | all types |
| `grid.gin` | `.gin` | dumbbells, vacancies, jumps |
| `formd_new_ndm_org.str` | `.gin` | types `4`, `6`, `67`, `7`, `77` |
| `geom_org.data` | as `geom_new.data` | all types |

For the cluster types (`6`, `67`, `7`, `77`), `geom_org.data` lists the
reference crystal *without* the removed atoms, while its first line still gives
the full number. Use `formd.lmp` or `formd_new_ndm_org.str` when you need the
perfect crystal.

### Visualisation files

These are `xyz` files (number of atoms, a blank line, then `label x y z`). The
labels are markers, not chemistry:

| Label | Meaning |
|---|---|
| `Fe` | an atom of the defect (in the vacancy files, a vacant site) |
| `Cu` | an atom of the surrounding crystal |
| `He`, `H` | a lattice site emptied by the defect |

For the cluster types, `He` marks the removed first neighbours of a C15 unit
or the removed neighbours of an icosahedral unit, and `H` the removed second
neighbours of a C15 unit.

The letters of the file names say what each file contains: `u` only the defect,
`v` with the emptied sites, `b` with the bulk, `c` cut to a box around the
defect, `j` without the trailing cell lines.

| File | Content | Written for |
|---|---|---|
| `sia.xyz` | The whole defective crystal: defect atoms as `Fe`, the rest as `Cu`. The cell vectors are appended after the atoms | all types |
| `jsia.xyz` | The same without the appended cell vectors | all types |
| `usia.xyz` | The defect atoms only, cell vectors appended | all types |
| `jusia.xyz` | The defect atoms only | all types |
| `juvsia.xyz` | Dumbbell atoms and, as `He`, the sites they replace, each followed by its integer coordinates | dumbbells, jumps |
| `uvsia.xyz` | Inserted atoms and emptied sites | cluster types |
| `uvbsia.xyz` | Inserted atoms, emptied sites and the complete reference crystal | cluster types |
| `uvbcsia.xyz` | As `uvbsia.xyz`, restricted to a box around the defect | cluster types |
| `uvbcsia_sim.xyz` | The defect and the crystal atoms around it, emptied sites left out: the real local structure | dumbbells, jumps, cluster types |
| `uvbcsia_cube.xyz` | As `uvbcsia_sim.xyz`, in a cubic box | dumbbells, jumps |

"Cluster types" are `6`, `67`, `7` and `77`. For these, `uvbcsia_sim.xyz` has a
faulty header: the first line overstates the number of atoms and the blank
second line is missing. Correct the header before loading the file in a viewer.

For vacancies (type `5`), `usia.xyz` and `jusia.xyz` list the vacant sites, and
`sia.xyz` and `jsia.xyz` list all the sites of the reference crystal with the
vacant ones labelled `Fe`. The crystal with the atoms actually removed is in
`formd_new.str`, `formd_new_ndm_scl.str`, `formd_new.lmp` and `geom_new.data`.

### Other files

| File | Content |
|---|---|
| `unit_cell`, `unit_vect` | The defective crystal and its cell vectors in an older layout (atom index, coordinates, three flags `1 1 1`) |
| `deb_loop.gin_0`, `fin_loop.gin_0` | Initial and final states of a jump, `.gin` layout (types `22`, `221`); the same states are in `deb_loop.lmp` and `fin_loop.lmp` |
| `k_dfct.gin`, `k_bulk.gin`, `k_dfct_vide.gin` | The rotated cut of type `67` |
| `fort.77` | List of the inserted atoms (cluster types) |
| `fort.567` | For each entry, the number of other entries at the first to fifth neighbour distance of bcc. Meaningful for types `2` and `3` only |
| `fort.42` | A header only; not used (cluster types) |

## Output of the other programs

| Program | File | Content |
|---|---|---|
| `cube.x` | `formd.str` | the reference crystal |
| | `unit_cell`, `unit_vect` | the same crystal in the older layout |
| | `fort.11` | one line per atom: `1 x y z 0 0 0` |
| `pre_put_vct.x` | `no00_new.inp` | the sites of the loop, in the plain-list layout of `no00.inp` |
| | `fort.12` | the sites of the cut plane, then of the loop |
| | `formd_new.str`, `formd_new_ndm_scl.str`, `geom_new.data`, `sia.xyz`, `unit_cell`, `unit_vect` | its own build of the loop ([programs.md](programs.md#pre_put_vctx-the-sites-of-a-dislocation-loop)) |
| `read_write_n00.x` | `no00_new.inp` | the shifted list |
| `cube_shell.x` | `formd.str`, `formd.str_inner`, `formd.gin`, `fort.11` | see [programs.md](programs.md#cube_shellx-a-reference-crystal-in-two-blocks) |

## Cleaning up

A run leaves more than twenty files behind. The `clean.sh` scripts of the
example directories remove them; copy one into your working directory.
