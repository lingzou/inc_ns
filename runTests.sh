#!/bin/bash

# run all cases
./build/inc_exec -i cylinder.i

EXODIFF=/Users/lingzou/projects/SAM/moose/framework/contrib/exodiff/exodiff

temp=$($EXODIFF output/output.e expected/output_19_cylinder.e)
error=$?
if [ $error -eq 0 ]
then
  printf "%24b................" "\e[0;32m[OK]\e[m" "\n"
else
  echo "$temp"   # optional: show exodiff's own detail on failure
  printf "%24b........" "\e[1;31m[DIFF ERROR]\e[m" "\n"
fi