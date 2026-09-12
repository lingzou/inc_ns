#!/bin/bash

EXODIFF=/Users/lingzou/projects/SAM/moose/framework/contrib/exodiff/exodiff

PROBLEMS=(cylinder cavity backward_facing_step cavity_buoyancy poiseuille)

for problem in "${PROBLEMS[@]}"
do
  log=$(mktemp)
  ./build/inc_exec -i input/${problem}.i > "$log" 2>&1

  temp=$($EXODIFF output/${problem}_out.e expected/${problem}_out.e)
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
