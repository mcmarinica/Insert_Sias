# Troubleshooting

## Check every result

The programs do few consistency checks, so look at each result before using
it.

1. **Count the atoms.** The first line of `formd_new.str` must be the number
   you expect: the reference count plus the net number of interstitials, or
   minus the number of vacancies.
2. **Look at it.** Open `uvbcsia_sim.xyz` (the defect with its surroundings) or
   `usia.xyz` (the defect alone) in a viewer.
3. **Check the shortest distance.** With the reader of
   [concepts.md](concepts.md#converting-to-other-formats):

   ```python
   from ase.neighborlist import neighbor_list
   atoms = read_gin("formd_new_ndm_scl.str", "Fe")
   print(neighbor_list("d", atoms, 3.0).min())
   ```

   At a0 = 2.8553 Å a sound configuration gives between about 1.3 Å (the two
   atoms of a dumbbell) and the first-neighbour distance (2.47 Å in bcc,
   2.02 Å in fcc). A value near zero means that atoms overlap.

## Error messages

| Message | Cause | Remedy |
|---|---|---|
| `there is a big problem in input` / `Check check_site`, followed by two numbers (sites found, sites asked) | A listed triplet is not a lattice site of the reference crystal | Check the parity rule of the lattice, that the site is inside the box, and that the crystal was built with a0 = 2.8553 |
| `Check check_site 3*` | Type `4`: a neighbour site of the di-interstitial is missing | Move the site away from the faces of the box |
| `Big problem in read_config_sia` / `The cell is not big enough. Increase the size!` | An index is larger than the box | Use a larger cell or shift the list |
| `There are atoms which exceed the box .... change the input` | A C15 or icosahedral unit reaches below index 0 | Move the cluster towards the middle of the box |
| `Check <check_site> 10*` | Some sites that a C15 unit should empty do not exist, usually near an upper face. **The program does not stop** and the result is wrong | Move the cluster towards the middle of the box |
| `There is no implementation for this octa type` | Type `6`: a unit type other than 1, 2 or 4 | Check the lines after the list of centres |
| `all the line should have the same parity` | Types `7`, `77`: a centre mixes even and odd indices | Use three odd indices |
| `as100 too small compared to a0. Increase as100` | Type `31` with the default dumbbell length | Raise `as100` above 1.90 Å in the source and recompile |
| `as111 too small ...`, `Reduce as111 ...`, `No solution for this value of as111` | Types `21`, `221`: dumbbell length outside the allowed window | See the window in [defects.md](defects.md#dumbbells-with-a-relaxed-row-types-21-and-31) |
| `nfix + nmig should be equal to nsia` and other messages naming `no00m.inp` | Types `22`, `221`: inconsistent header or flags | Check the [layout of `no00m.inp`](defects.md#dumbbell-jumps-types-22-and-221) |
| `This software it isi not designed for tsia` | Unknown defect code | Use a code of the [catalogue](defects.md#overview) |
| `end-of-file during read` (Intel) or `End of file` (gfortran) | An input is shorter than announced: fewer sites than the count on line 1, missing unit types for type `6`, missing answers on standard input, or a missing input file | Complete the input. If `formd.str` or `no00.inp` did not exist, the program has left an empty file of that name; delete it |
| `n1MAX-1 and n2MAX-1 should be equal and multiple of 4`, `At least the n3MAX should be increased` | Type `67`: unsuitable reference cell | See [type 67](defects.md#type-67) |
| `ERROR <lammps_end>: formd_new.lmp is inconsistent` | The number of atoms written differs from the number announced. The other output files have the same flaw: the structure itself is wrong | Look for overlapping atoms, typically an icosahedral unit on an even centre |
| `WARNING <lammps_begin>: ... is not written` | The cell of `formd.str` is one LAMMPS cannot hold without a rotation: its first vector is not along x or its second not in the xy plane. The other files are written | Use a reference crystal from `cube.x` |

One message can be ignored:
`WARNING ... the analyze_config SIA .... Not yet implemented for tsia` only
concerns the neighbour statistics of `fort.567`.

## Known limitations

These were observed with the current sources.

**The lattice parameter is not checked.** A reference crystal built with a
value other than 2.8553 Å is accepted, and the defect may be put at the wrong
site without any message.

**Repeated sites are not detected.** A site listed twice receives two defects
on top of each other.

**Clusters are not periodic.** Types `4`, `6`, `67`, `7` and `77` do not wrap
through the cell boundaries. Dumbbells and vacancies can sit anywhere, and the
row relaxation of types `21`, `31` and `221` does follow the periodic images.

**Dumbbell lengths are those of bcc.** In fcc the compiled 2.05 Å of the
`<110>` dumbbell is longer than the first-neighbour distance, and each of its
atoms ends up 0.99 Å from the next atom of the crystal along [110]. With the
fcc values this distance is 1.24 Å. Compile the fcc values for fcc loops
([install.md](install.md#compiled-parameters)).

**Type 31 does not run with the default parameters.** It needs a `<100>`
dumbbell length above 1.90 Å, and the compiled value is 1.70 Å.

**Type 7 needs odd centres.** A centre with three even indices passes the
parity test, but only the central atom is removed while thirteen are inserted,
and atoms overlap. The run ends with
`ERROR <lammps_end>: formd_new.lmp is inconsistent`.

**Some files of the cluster types have inconsistent headers.** In
`uvbcsia_sim.xyz` the number of atoms is overstated and the blank second line
is missing. In `geom_org.data` the first line gives the full reference count
although the removed atoms are not listed. For type `77`, the sites of the
extra dumbbells are not marked in `uvsia.xyz`, `uvbsia.xyz` and `uvbcsia.xyz`,
whose first line overstates the number of atoms. The structure files
(`formd_new.str`, `formd_new_ndm_scl.str`, `formd_new.lmp`) are not affected.

**Bi-plane loops are for fcc.** The shift between the two planes in
`pre_put_vct.x` is the fcc vector (1, 0, 1) a0/2. In bcc the program writes
`no00_new.inp` and then stops with `there is a big problem in input`; the file
is not usable. Set `deca` to (1, 1, 1) in the source for bcc.

**Some printed limits are wrong.** For a hexagon in the (111) plane,
`pre_put_vct.x` announces a maximum size of 0.

**`fort.567` is meaningful for types 2 and 3 only.** The neighbour distances it
uses are not set for the other types.

**File names are fixed.** Two runs in the same directory overwrite each other.
