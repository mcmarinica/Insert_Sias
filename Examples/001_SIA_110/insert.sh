#this script will insert an 110 SIAs at the position specifiend in the file no00.inp
#prepare the input for cube.x 
#which create a bcc crystal (option 2) 
#with the lattice parameter a0=2.8553 and 
#dimensions 4a0 x4a0 x4a0


cat > cube.inp <<EOF
2
2.8553
8 8 8
EOF

../../Src/bin/cube.x < cube.inp



#How many SIAs and in which lattice site.: here 1 SIA at (2,2,2) center inb a0/2 units
cat > no00.inp <<EOF
1
1 1 1
EOF



# the type of 110 SIA is option 1:
cat > put.inp <<EOF
4
EOF


# then the position of the SIAs  is read from  no00.inp by put_vct.x
../../Src/bin/put_vct.x < put.inp


echo 'Output is obvius ... the gin file is   formd_new_ndm_scl.str'
