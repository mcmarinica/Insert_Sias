#!/bin/sh
# Vacancies and vacancy clusters in a bcc crystal.
#
#   sh make_vacancies.sh                 all configurations, 6x6x6 cubic cells
#   sh make_vacancies.sh 8               all configurations, 8x8x8 cubic cells
#   sh make_vacancies.sh 8 02v1 03v4     only the configurations named
#
# A configuration NNvM is a list of vacant sites, INPS_VAC/no00.inp_NNvM.
# The results go to a directory named after the cell, for instance 06x06x06:
#
#   bulk.lmp  bulk.gin     the perfect crystal
#   NNvM.lmp  NNvM.gin     the crystal with the vacancies of configuration NNvM
#
# .lmp is a LAMMPS data file, .gin holds reduced coordinates.

A0=2.8553        # lattice parameter; it must be the one compiled in put_vct.x
NMIN=4           # smallest cell that holds every site of INPS_VAC

here=$(pwd)
bin=$(cd ../../Src/bin && pwd)
for exe in cube.x put_vct.x; do
    [ -x "$bin/$exe" ] || { echo "$bin/$exe is missing: build the programs first"; exit 1; }
done

# the cell
n=${1:-6}
[ $# -gt 0 ] && shift
case $n in
    ''|*[!0-9]*) echo "the first argument is the number of cubic cells per side, not '$n'"; exit 1 ;;
esac
[ "$n" -ge $NMIN ] || { echo "the cell must have at least $NMIN cubic cells per side"; exit 1; }
nbulk=$((2 * n * n * n))

# the configurations: those named on the command line, or all
if [ $# -gt 0 ]; then
    labels="$*"
else
    labels=$(ls INPS_VAC | sed 's/^no00\.inp_//')
fi
for label in $labels; do
    [ -f "INPS_VAC/no00.inp_$label" ] || { echo "no configuration named $label in INPS_VAC"; exit 1; }
done

out=$(printf '%02dx%02dx%02d' "$n" "$n" "$n")
work=$out/work
mkdir -p "$work" && cd "$work" || exit 1

printf '2\n%s\n%s %s %s\n' $A0 "$n" "$n" "$n" | "$bin/cube.x" > cube.log

for label in $labels; do
    input=$here/INPS_VAC/no00.inp_$label
    nvac=$(head -1 "$input" | tr -d ' ')

    cp "$input" no00.inp
    rm -f formd_new.lmp
    echo 5 | "$bin/put_vct.x" > put.log        # 5 = vacancies

    # put_vct.x must have removed exactly the listed atoms
    natoms=$(sed -n 3p formd_new.lmp 2> /dev/null | cut -d' ' -f1)
    if [ "$natoms" != "$((nbulk - nvac))" ]; then
        echo "$label failed, see $work/put.log"
        exit 1
    fi

    cp formd_new.lmp         "../$label.lmp"
    cp formd_new_ndm_scl.str "../$label.gin"
    echo "$label: $natoms atoms ($nvac removed)  ->  $out/$label.lmp  $out/$label.gin"
done

cp formd.lmp ../bulk.lmp
cp grid.gin  ../bulk.gin
echo "perfect crystal: $nbulk atoms  ->  $out/bulk.lmp  $out/bulk.gin"

cd "$here" && rm -r "$work"
