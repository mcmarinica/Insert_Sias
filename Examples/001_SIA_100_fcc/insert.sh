# This script inserts one <100> dumbbell in an fcc crystal.
# The three blocks below (cell, sites, defect type) are the ones to edit.

# Input of cube.x: fcc (option 1), lattice parameter a0 = 2.8553,
# 10 x 10 x 10 cells (4000 atoms).
cat > cube.inp <<EOF
1
2.8553
10 10 10
EOF

../../Src/bin/cube.x < cube.inp

# Sites of the dumbbells, in units of a0/2: here one dumbbell at (2,2,2).
# In fcc the three indices of a site have an even sum.
cat > no00.inp <<EOF
1
2 2 2
EOF

# Defect type: 3 is the <100> dumbbell.
cat > put.inp <<EOF
3
EOF

../../Src/bin/put_vct.x < put.inp

echo 'Result: formd_new.str (and formd_new.lmp for LAMMPS); the gin file is formd_new_ndm_scl.str'
