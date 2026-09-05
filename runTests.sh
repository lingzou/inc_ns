#!/bin/bash

EXODIFF=/Users/lingzou/projects/SAM/moose/framework/contrib/exodiff/exodiff

PROBLEMS=(cylinder cavity)

for problem in "${PROBLEMS[@]}"
do
  ./build/inc_exec -i input/${problem}.i

  temp=$($EXODIFF output/${problem}_out.e expected/${problem}_out.e)
  error=$?
  if [ $error -eq 0 ]
  then
    printf "%-10s %b\n" "$problem" "\e[0;32m[OK]\e[m"
  else
    echo "$temp"   # optional: show exodiff's own detail on failure
    printf "%-10s %b\n" "$problem" "\e[1;31m[DIFF ERROR]\e[m"
  fi
done
