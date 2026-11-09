# This script inserts one <110> dumbbell in a bcc crystal.
# The three blocks below (cell, sites, defect type) are the ones to edit.

# Input of cube.x: bcc (option 2), lattice parameter a0 = 2.8553,
# 8 x 8 x 8 cells (1024 atoms).
cat > cube.inp <<EOF
2
2.8553
8 8 8
EOF

../../Src/bin/cube.x < cube.inp

# Sites of the dumbbells, in units of a0/2: here one dumbbell at (1,1,1).
cat > no00.inp <<EOF
1
1 1 1
EOF

# Defect type: 1 is the <110> dumbbell.
cat > put.inp <<EOF
1
EOF

../../Src/bin/put_vct.x < put.inp

echo 'Result: formd_new.str (and formd_new.lmp for LAMMPS); the gin file is formd_new_ndm_scl.str'
