#!/bin/bash

# run all cases
"./inc_exec"

# do diff for each case
#temp=$(diff output/FluentTwoDMesh.vtu expected/FluentTwoDMesh.vtu)
#temp=$(diff output/U_STAR.vtu expected/U_STAR.vtu)
temp=$(diff output/output_19.vtu expected/output_19_cylinder.vtu)
error=$?
if [ $error -eq 0 ]
then
  printf "%24b................" "\e[0;32m[OK]\e[m" "\n"
else
  printf "%24b........" "\e[1;31m[DIFF ERROR]\e[m" "\n"
fi
