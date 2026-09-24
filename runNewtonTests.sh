#!/bin/bash
#
# Validation for the monolithic Newton-Krylov path ("-newton true"). Kept separate from
# runTests.sh, which regression-tests the default segregated solver against expected/*.e
# and must stay untouched.
#
# This does NOT exodiff the two modes against each other. They are different discrete
# problems: the segregated scheme builds the face mass flux from the predictor velocity
# u*, the monolithic residual from the velocity that actually solves momentum, and the
# two agree only to O(dt). See docs/newton_krylov_formulation.tex, "...but they do not
# share a fixed point". What is checked instead is the RATE: halve dt, and the
# difference between the modes should halve too.
#
# Run from the repo root.

set -u

PROBLEM=cavity            # input/cavity.i IS the lid-driven cavity (lid = wall with -bc_5_u 1.0)
END_TIME=0.2              # held fixed across the refinement so the runs are comparable
DTS=(0.01 0.005 0.0025)

if [ ! -x ./build/inc_exec ]; then
  echo "ERROR: ./build/inc_exec not found. Build first (see compile.sh)." >&2
  exit 1
fi

args=()
tmp_inputs=()
for dt in "${DTS[@]}"; do
  n=$(python3 -c "print(round($END_TIME/$dt))")
  tag="_newtontest_${PROBLEM}_dt${dt}"
  in="input/${tag}.i"
  sed -e "s|^-dt .*|-dt ${dt}|" -e "s|^-n_timestep .*|-n_timestep ${n}|" \
      "input/${PROBLEM}.i" > "$in"
  tmp_inputs+=("$in")

  for mode in segregated newton; do
    log=$(mktemp)
    if [ "$mode" = newton ]; then
      ./build/inc_exec -i "$in" -newton true > "$log" 2>&1
    else
      ./build/inc_exec -i "$in" > "$log" 2>&1
    fi
    if [ $? -ne 0 ]; then
      echo "--- $mode run failed at dt=$dt ---"
      tail -20 "$log"
      rm -f "$log" "${tmp_inputs[@]}"
      exit 1
    fi
    rm -f "$log"
  done
  printf "%-12s dt=%-8s %b\n" "$PROBLEM" "$dt" "both modes ran \e[0;32m[ok]\e[m"

  args+=("$dt" "input/${tag}_out.e" "input/${tag}_out_newton.e")
done

echo
python3 input/check_newton_order.py "${args[@]}"
status=$?

rm -f "${tmp_inputs[@]}"
for dt in "${DTS[@]}"; do
  rm -f input/_newtontest_${PROBLEM}_dt${dt}_out*.e \
        input/_newtontest_${PROBLEM}_dt${dt}_out*_convergence.csv
done

exit $status
