# Concepts

## The three steps

1. **Reference crystal.** `cube.x` writes a perfect bcc or fcc supercell to
   `formd.str`. The cell is orthogonal, its edges are along [100], [010] and
   [001], and one atom sits at the origin.
2. **Defect centres.** `no00.inp` lists the lattice sites that will carry a
   defect.
3. **Insertion.** `put_vct.x` replaces the atom at each listed site by the
   defect you chose, and writes the new crystal.

"Inserting a dumbbell at a site" means: the atom of that site is removed and
two atoms are put in its place, symmetrically about the site. One entry in
`no00.inp` therefore adds one atom. The larger objects (C15 unit, icosahedral
unit, Gao di-interstitial) work the same way, with several atoms removed and
several inserted around the listed site.

A cluster or a loop is a list of sites. For a loop of 55 parallel `<111>`
dumbbells, `no00.inp` holds the 55 sites of a disc cut in a (111) plane, and
`put_vct.x` puts one dumbbell on each.

## Integer site coordinates

Sites are named by three integers, their Cartesian coordinates in units of
a0/2:

```
(n1, n2, n3) = (x, y, z) / (a0/2)
```

In a cell of N1 x N2 x N3 cubic unit cells, each index runs from 0 to 2N - 1.

| Lattice | Which triplets are lattice sites | Examples |
|---|---|---|
| bcc | all three even (cube corners) or all three odd (cube centres) | `4 4 4`, `3 3 3`, `5 3 1` |
| fcc | the sum n1 + n2 + n3 is even | `4 4 4`, `3 3 2`, `5 4 1` |

A triplet that is not a lattice site is rejected by most defect types (see
[troubleshooting.md](troubleshooting.md)). Two types use off-lattice triplets on
purpose: the centre of an icosahedral unit in fcc is an octahedral site (all
three indices odd), and type `77` can add single atoms at arbitrary positions.

To put a defect in the middle of a 10 x 10 x 10 cell, use a site near
`10 10 10`.

## The compiled lattice parameter

`put_vct.x` and `pre_put_vct.x` find the integer coordinates of the atoms by
dividing their positions by a0/2, with a0 = 2.8553 Å fixed in the source. The
reference crystal must therefore be built with exactly this value:

- for bcc and for fcc,
- whatever the element you have in mind.

`cube.x` accepts any lattice parameter, and `put_vct.x` does not check it. With
another value the sites are misidentified: the program either stops with
`there is a big problem in input`, or runs and puts the defect in the wrong
place.

The same holds for the atom names: every output file labels the atoms `Fe` (and
`Cu`, `He`, `H` as markers in the visualisation files). These labels carry no
chemical meaning.

To obtain a crystal of another element, build the defect at 2.8553 Å and
rescale the result. The file `formd_new_ndm_scl.str` gives the positions in
reduced coordinates, so rescaling amounts to multiplying the three cell vectors
by `a_target / 2.8553`. This also scales the dumbbell lengths, which is harmless
for a configuration that will be relaxed.

## Using the result in LAMMPS

`put_vct.x` writes two LAMMPS data files at every run:

| File | Content |
|---|---|
| `formd.lmp` | the reference crystal, the same atoms as `formd.str` |
| `formd_new.lmp` | the crystal with the defects, the same atoms in the same order as `formd_new.str` |

They use atom style `atomic`, a single atom type and lengths in Å. In a LAMMPS
input:

```
units       metal
atom_style  atomic
boundary    p p p
read_data   formd_new.lmp
mass        1 183.84
```

The files hold no mass and no element: set the mass and the pair style for
your metal. The geometry is the one of the reference lattice parameter,
2.8553 Å. To bring it to the lattice parameter of your potential, rescale the
box and the atoms together, for instance for 3.165 Å:

```
variable    s equal 3.165/2.8553
change_box  all x scale ${s} y scale ${s} z scale ${s} remap
```

`change_box` scales the box about its centre. The pair `formd.lmp` and
`formd_new.lmp`, rescaled alike, gives the perfect and the defective crystal
that a formation energy needs.

Some atoms of a defect near a face of the box may lie slightly outside it;
LAMMPS puts them back through the periodic boundaries when it reads the file.

## Converting to other formats

Besides the LAMMPS files, the package writes its own text formats and some
`xyz` files ([files.md](files.md)). The POSCAR and CFG files found in some
example directories were made with separate tools that are not part of this
repository.

A few lines of Python with [ASE](https://wiki.fysik.dtu.dk/ase/) read
`formd_new_ndm_scl.str`, rescale it and write any format ASE knows:

```python
import numpy as np
from ase import Atoms
from ase.io import write

A0_CODE = 2.8553          # lattice parameter compiled into put_vct.x

def read_gin(path, symbol):
    """Read a file in the .gin layout, such as formd_new_ndm_scl.str."""
    with open(path) as f:
        lines = [l for l in f if l.strip()]
    cell = np.array([l.split() for l in lines[1:4]], dtype=float)
    natoms = int(lines[4])
    frac = np.array([l.split()[:3] for l in lines[5:5 + natoms]], dtype=float)
    return Atoms(symbol * natoms, cell=cell, scaled_positions=frac, pbc=True)

atoms = read_gin("formd_new_ndm_scl.str", "W")
a_target = 3.165                                   # the lattice parameter you want
atoms.set_cell(atoms.cell * a_target / A0_CODE, scale_atoms=True)
write("POSCAR", atoms, format="vasp")              # or "defect.xyz", "defect.cif", ...
```

In every output file the inserted atoms come first and the remaining bulk atoms
follow, so the first lines of the converted file are the defect.

## What the programs do not do

- **No relaxation.** Positions are geometric guesses.
- **No volume change.** `put_vct.x` keeps the cell of the reference crystal.
- **No chemistry.** One species only.
- **One defect type per run.** All entries of `no00.inp` get the same type and
  the same orientation. The only built-in mixture is type `77`.
- **One orientation per type.** The dumbbells are along [110], [111] or [001].
  Other variants of the same family require a change in the source
  (`read_sia` in `Src/lib/point_defects.f90`).
