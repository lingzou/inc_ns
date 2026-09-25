#!/bin/bash

EXODIFF=/Users/lingzou/projects/SAM/moose/framework/contrib/exodiff/exodiff

# exodiff's default floor is 0.0, so it compares even physically-near-zero fields (a
# pinned pressure anchor row, vel_y in a straight channel) by pure relative difference,
# where roundoff-level absolute noise reads as a triple-digit percentage. Item 15: a
# numerically-neutral change to updateMassVeclocities() (item 12) flipped poiseuille's
# vel_y ~1e-7 by 9.5e-15 and that alone tripped this script. -Floor 1e-6 treats any pair
# of values that are both <= 1e-6 in magnitude as equal outright, which comfortably
# covers that class of noise (fields there sit at ~1e-7 to 1e-9) while leaving real
# differences undisturbed -- checked against the KSPCG->KSPGMRES switch (up to ~1e-2
# relative on vel_x), which is still flagged as different up to -Floor 1e-1.
FLOOR=1e-6

PROBLEMS=(cylinder cavity backward_facing_step cavity_buoyancy poiseuille)

for problem in "${PROBLEMS[@]}"
do
  log=$(mktemp)
  ./build/inc_exec -i input/${problem}.i > "$log" 2>&1

  temp=$($EXODIFF -Floor $FLOOR input/${problem}_out.e expected/${problem}_out.e)
  error=$?
  if [ $error -eq 0 ]
  then
    printf "%-10s %b\n" "$problem" "\e[0;32m[OK]\e[m"
  else
    echo "--- simulation output ($log) ---"
    cat "$log"
    echo "$temp"   # optional: show exodiff's own detail on failure
    printf "%-10s %b\n" "$problem" "\e[1;31m[DIFF ERROR]\e[m"
  fi
  rm -f "$log"
done
