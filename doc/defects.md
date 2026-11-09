# Defect catalogue

`put_vct.x` asks for one number, the defect type. This page describes each
type: what is built, which lattice it is meant for, and how to write the list
of sites.

Throughout, sites and offsets are written in integer coordinates (units of
a0/2, see [concepts.md](concepts.md#integer-site-coordinates)), and lengths in
Å refer to the compiled values of [install.md](install.md#compiled-parameters).

## Overview

| Code | Defect | Lattice | Input | Removed / inserted per entry |
|---|---|---|---|---|
| `1` | `<110>` dumbbell | bcc, fcc | `no00.inp`, plain list | 1 / 2 |
| `2` | `<111>` dumbbell | bcc, fcc | `no00.inp`, plain list | 1 / 2 |
| `21` | `<111>` dumbbell, row neighbours spread | bcc | `no00.inp`, plain list | 1 / 2 |
| `210` | `<111>` dumbbell for Frank loops | fcc | `no00.inp`, plain list | 1 / 2 |
| `3` | `<100>` dumbbell | bcc, fcc | `no00.inp`, plain list | 1 / 2 |
| `31` | `<100>` dumbbell, row neighbours spread | bcc | `no00.inp`, plain list | 1 / 2 |
| `22` | `<111>` dumbbell jump | bcc | `no00m.inp` | 1 / 2 |
| `221` | `<111>` dumbbell jump, row neighbours spread | bcc | `no00m.inp` | 1 / 2 |
| `4` | Gao di-interstitial | bcc | `no00.inp`, plain list, and `gao_unit.in` | 4 / 6 |
| `5` | vacancy | bcc, fcc | `no00.inp`, plain list | 1 / 0 |
| `6` | C15 unit | bcc | `no00.inp`, sites then unit types | 10 / 12 |
| `67` | C15 cluster plus a rotated cut | bcc | as `6` | as `6` |
| `7` | icosahedral unit | fcc | `no00.inp`, plain list of centres | 6 / 13 |
| `77` | icosahedral units plus extras | fcc | `no00.inp`, two blocks | see below |

## The plain list

Most types read `no00.inp` in this form: the number of entries, then one site
per line.

```
3
4 4 4
5 5 5
6 6 6
```

The numbers are read in free format. Every site must be a lattice site of the
reference crystal, and no site may appear twice.

## Dumbbells: types 1, 2 and 3

The atom of each listed site is replaced by two atoms placed symmetrically
about the site.

| Code | Axis | Distance between the two atoms |
|---|---|---|
| `1` | [110] | 2.05 Å |
| `2` | [111] | 1.90 Å |
| `3` | [001] | 1.70 Å |

One entry gives a single interstitial. A list of neighbouring sites in one
plane gives a dislocation loop of parallel dumbbells;
[`pre_put_vct.x`](programs.md#pre_put_vctx-the-sites-of-a-dislocation-loop)
writes such lists.

The compiled distances suit bcc. For loops in fcc, switch to the fcc values
first ([install.md](install.md#compiled-parameters)).

## Dumbbells with a relaxed row: types 21 and 31

Type `21` builds the `<111>` dumbbell of type `2` and also moves the four
nearest atoms of the same [111] row, at the site ± (1,1,1) and ± (2,2,2), away
from the dumbbell. The extra atom is thus spread along the row, which is closer
to a crowdion than a bare dumbbell is. The row is followed through the
periodic boundaries.

The displacements are derived from the dumbbell length, a0 and the number p of
sites along the row in the cell (printed as `periodicity in <111> direction`).
The dumbbell length must exceed 1.65 Å and stay below 2.47 Å x p/(p+1); the
default of 1.90 Å qualifies.

Type `31` does the same for the `<100>` dumbbell, moving the neighbours at
n3 ± 2 and n3 ± 4 along [001]. It needs a dumbbell length above a0/1.5 =
1.90 Å. **With the default 1.70 Å it stops** with
`as100 too small compared to a0`; raise `as100` in the source to use it.

## Frank loops: type 210

As type `2`, except that one of the two atoms is shifted by an extra
a0/6 [2 -1 -1]. The dumbbell atoms are at

```
site - 0.95 Å [111]/√3
site + 0.95 Å [111]/√3 + a0/6 [2 -1 -1]
```

Applied to a disc of sites in a (111) plane of fcc, this adds a layer in a
faulted stacking position, the starting point for a Frank loop. Give it the
sites of a (111) loop from `pre_put_vct.x`.

## Vacancies: type 5

The atoms of the listed sites are removed. A list of several sites gives a
vacancy cluster or, with a planar list, a vacancy loop.

```
2
4 4 4
5 5 5
```

This is a di-vacancy on first neighbours in bcc.

## Gao di-interstitial: type 4

Four atoms are removed, the listed site and its neighbours at the site
+ (1,-1,-1), (-1,-1,1) and (1,1,1), and six atoms are inserted. The positions
of the six atoms are not built in: they are read from `gao_unit.in`, six lines
of `x y z` relative to the listed site, **in units of a0**.

A template:

```
   0.711532006526317  -0.711529095517008  -0.517605603651929
   0.187231232167829  -0.187231821140325  -0.372279359148133
   0.187518230675831   0.370903277326304   0.186737734733335
  -0.370903866298800  -0.187518819648328   0.186737734733335
  -0.518477099121354  -0.711126597609444   0.710912509746847
   0.711126008636947   0.518476510148858   0.710912509746847
```

Since the six positions come from a file, this type also serves for any other
compact object that replaces these four sites by six atoms.

## C15 clusters: types 6 and 67

### The unit

A C15 unit is centred on a lattice site c. Ten atoms around c are removed and
twelve are inserted, a net gain of two: one unit alone is a di-interstitial.

| | Type 1 | Type 2 |
|---|---|---|
| Removed first neighbours (a tetrahedron) | c + (-1,-1,-1), (1,1,-1), (1,-1,1), (-1,1,1) | c + (1,1,1), (-1,-1,1), (1,-1,-1), (-1,1,-1) |
| Removed second neighbours | c ± (2,0,0), (0,2,0), (0,0,2) | the same |
| Inserted | three atoms around each corner of the tetrahedron | the same, around the other tetrahedron |

The atom at c itself stays. The two types are the two possible orientations of
the tetrahedron; type 2 is type 1 rotated by 90° about z.

### Clusters

Larger clusters are made of linked units. Linked units alternate in type, and
the centre of the neighbour lies at a corner of the tetrahedron that is *not*
removed:

- a type 1 unit at c links to type 2 units at c + (1,1,1), (-1,-1,1),
  (1,-1,-1), (-1,1,-1);
- a type 2 unit at c links to type 1 units at c + (-1,-1,-1), (1,1,-1),
  (1,-1,1), (-1,1,1).

Linked units share removed sites and inserted atoms. The program removes the
duplicates and prints the result:

```
 The number of DIFFERENT VACANCIES left in the input ........:           14
 The number of DIFFERENT GAOS atoms left in input ............:           18
 The no of defects SIAs ...................:           4
```

One unit gives 2 interstitials and two linked units give 4. Beyond that the
count depends on the shape of the cluster; closed rings share more atoms.

### Input

After the list of centres come the unit types, one per line, in the same
order:

```
2
5 5 5
6 6 6
1
2
```

This is the four-interstitial cluster. The `no00.inp` of the example
directories `004a_SIA_C15`, `011a_SIA_C15`, `017a_SIA_C15` and `101_SIA_c15`
are such files for 4, 11, 17 and 101 interstitials
([examples.md](examples.md#examples004a_sia_c15-011a_sia_c15-017a_sia_c15-101_sia_c15)).
They use type 1 on odd centres and type 2 on even ones.

There is a third unit type, `4`: a partial unit that removes two sites and
inserts four atoms. The source describes it as the building block for clusters
with an odd number of interstitials. None of the supplied inputs use it.

### Restrictions

The cluster is not wrapped through the periodic boundaries. Every removed site
must exist in the box, so keep the centres at least three index units away
from all six faces; the middle of the box is the natural place. To move a
cluster, use [`read_write_n00.x`](programs.md#read_write_n00x-moving-a-list-of-sites)
with answer `2`.

### Type 67

Type `67` builds the same cluster as type `6`, then cuts a box out of the
result with its edges along [110], [-110] and [001], and writes it to
`k_dfct.gin`, with the matching perfect crystal in `k_bulk.gin` and the
cluster atoms alone in `k_dfct_vide.gin`.

This is a special-purpose option with compiled sizes. The reference crystal
must have the same odd number of cells along x and y (for example 49, 51 or
53) and be long enough along z, and the size of the cut is the constant
`pdist` in `Src/lib/kazuto_style.f90`. The program prints a suggested centre
for the cluster. Read the routine `write_out_kazuto` of that file before using
it.

## Icosahedral clusters in fcc: types 7 and 77

These types build the A15 clusters of A. M. Goryaeva *et al.*, Nat. Commun.
**14**, 3003 (2023): compact three-dimensional interstitial clusters with the
structure of the A15 Frank-Kasper phase.

### The unit

The centre c of a unit is an octahedral interstitial site of fcc: **three odd
indices**, for example `3 3 3`. Around it:

- the six nearest atoms, at c ± (1,0,0), (0,1,0), (0,0,1), are removed;
- each is replaced by a pair of atoms, giving twelve atoms at
  c + (±1, 0, ±½), (0, ±½, ±1), (±½, ±1, 0);
- one atom is added at c.

The twelve atoms form an icosahedral cage around the central one. Six atoms
out, thirteen in: one unit holds 7 interstitials.

### Clusters

Units whose centres are two index units apart along a cube axis share atoms.
The program merges them as for C15:

| Centres | Interstitials |
|---|---|
| `3 3 3` | 7 |
| `3 3 3`, `5 3 3` | 13 |
| `3 3 3`, `5 3 3`, `5 3 5` | 19 |
| `3 3 3`, `5 3 3`, `5 3 5`, `3 3 5` | 24 |

### Input for type 7

The plain list, holding the centres:

```
2
3 3 3
5 3 3
```

Restrictions:

- The three indices of a centre must have the same parity, otherwise the
  program stops. Use odd indices. Even indices are accepted but give
  overlapping atoms ([troubleshooting.md](troubleshooting.md#known-limitations)).
- As for C15, the cluster is not wrapped through the periodic boundaries; keep
  it away from the faces of the box.

### Type 77: icosahedral units with extras

Type `77` adds isolated `<100>` dumbbells and single atoms to the icosahedral
units, to decorate or extend a cluster. The header of `no00.inp` holds two
numbers, and two blocks follow:

```
1 3
5 5 5
8 5 5 3
5 5 2 2
7 5 5 4
```

- First line: the number of icosahedral units, then the number of extras.
- First block: the centres of the units, as for type `7`.
- Second block: one extra per line, `n1 n2 n3 t`.

| `t` | Effect | Requirement |
|---|---|---|
| `1`, `2`, `3` | The atom at the site is replaced by a `<100>` dumbbell along x, y, z (1.70 Å) | The site must be a lattice site |
| `4` | One atom is added at the position (n1, n2, n3) a0/2 | The position must be empty, for instance an octahedral site |

The example above gives one unit (+7), two dumbbells (+2) and one added atom
(+1): ten atoms more than the reference crystal. The number of extras may be 0.

## Dumbbell jumps: types 22 and 221

These types write the two end states of a `<111>` dumbbell that jumps to the
next site of its row, from (n1, n2, n3) to (n1+1, n2+1, n3+1). They read
`no00m.inp` instead of `no00.inp`:

```
2 1 1
4 4 4 1
6 4 4 0
```

- First line: the number of dumbbells, how many of them jump, how many stay.
  The last two must add up to the first.
- Then one line per dumbbell: the site and a flag, `1` if it jumps, `0` if it
  stays.

| Files | Content |
|---|---|
| `deb_loop.gin_0`, `deb_loop.lmp` | Initial state: every dumbbell on its listed site |
| `fin_loop.gin_0`, `fin_loop.lmp` | Final state: the flagged dumbbells one site further along [111] |

Each state is written twice, in reduced coordinates (`.gin_0`) and as a LAMMPS
data file (`.lmp`). All four files hold the same number of atoms in the same
order, so the two states can serve directly as the end points of a path
calculation. The other output files describe the initial state.

A jumping dumbbell may not have an index equal to 0 or to the largest index of
the box.

Type `221` adds the row relaxation of type `21` to both states.

## Names in the source

In the Fortran sources a defect type is never written as a number. Each code
is a named constant of the module `defect_types`
(`Src/lib/module_defect_types.f90`), whose value is the number typed at the
prompt.

| Code | Constant |
|---|---|
| `1` | `DEFECT_SIA_110` |
| `2` | `DEFECT_SIA_111` |
| `21` | `DEFECT_SIA_111_CARE` |
| `210` | `DEFECT_FRANK_LOOP` |
| `3` | `DEFECT_SIA_100` |
| `31` | `DEFECT_SIA_100_CARE` |
| `22` | `DEFECT_SIA_111_MIGRATION` |
| `221` | `DEFECT_SIA_111_MIGRATION_CARE` |
| `4` | `DEFECT_GAO` |
| `5` | `DEFECT_VACANCY` |
| `6` | `DEFECT_C15` |
| `67` | `DEFECT_C15_TWISTED` |
| `7` | `DEFECT_ICO` |
| `77` | `DEFECT_ICO_DUMBBELLS` |

The same module holds the functions that group the types the program treats
alike: `defect_is_dumbbell` (1, 2, 3, 21, 31, 210), `defect_is_migration`
(22, 221), `defect_is_c15` (6, 67), `defect_is_ico` (7, 77),
`defect_is_cluster` (6, 67, 7, 77) and `defect_is_known` (any code of the
table).

To add a defect type, give it a constant in this module, add it to
`defect_is_known` and to the group it belongs to, and handle it in
`Src/programs/put_vct.f90`.
