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
# Two cases, and the second one is not redundant:
#
#   cavity      closed box, no pressure B.C. anywhere. Exercises the homogeneous-Neumann
#               boundary branch of computeFaceFluxes() and the pressure anchor.
#
#   poiseuille  has a pressure_outlet, which is the ONLY way to reach the Dirichlet
#               branch of that function. Added after to-do item 1, where the Newton path
#               was found still carrying the pre-e13e2fb boundary flux: cavity alone
#               cannot see that class of bug, and did not. Measured on the unfixed code
#               the ratios below decay to 1.67 / 1.63 / 1.43; with the fix they sit at
#               3.53 / 3.33 / 3.84. Hence MIN_RATIO 2.0 for this case rather than the
#               default 1.5 -- it separates the two states from both sides with room to
#               spare, where the default would have caught the break by only 5%.
#
# Do NOT refine poiseuille past dt = 0.00125. At 6.25e-4 the mode-to-mode difference
# turns around and grows (vel_x 4.9e-4 -> 1.2e-2), identically in the fixed and unfixed
# code, so it is some other dt floor and not what this test is about.
#
# Run from the repo root.

set -u

if [ ! -x ./build/inc_exec ]; then
  echo "ERROR: ./build/inc_exec not found. Build first (see compile.sh)." >&2
  exit 1
fi

# run_case <problem> <end_time> <min_ratio> <dt> [<dt> ...]
# end_time is held fixed across the refinement so the runs are comparable.
run_case()
{
  problem=$1
  end_time=$2
  min_ratio=$3
  shift 3
  dts=("$@")

  args=()
  tmp_inputs=()
  status=0

  for dt in "${dts[@]}"; do
    n=$(python3 -c "print(round($end_time/$dt))")
    tag="_newtontest_${problem}_dt${dt}"
    in="input/${tag}.i"
    # -output_interval is forced to n so the last Exodus record is the final step at
    # every dt, which is the record check_newton_order.py reads. cavity.i does not set
    # the key at all, so the expression is a no-op there and it keeps its default of 1.
    sed -e "s|^-dt .*|-dt ${dt}|" -e "s|^-n_timestep .*|-n_timestep ${n}|" \
        -e "s|^-output_interval .*|-output_interval ${n}|" \
        "input/${problem}.i" > "$in"
    tmp_inputs+=("$in")

    for mode in segregated newton; do
      log=$(mktemp)
      if [ "$mode" = newton ]; then
        ./build/inc_exec -i "$in" -newton true > "$log" 2>&1
      else
        ./build/inc_exec -i "$in" > "$log" 2>&1
      fi
      if [ $? -ne 0 ]; then
        echo "--- $problem: $mode run failed at dt=$dt ---"
        tail -20 "$log"
        status=1
      fi
      rm -f "$log"
    done
    if [ $status -ne 0 ]; then
      break
    fi
    printf "%-12s dt=%-8s %b\n" "$problem" "$dt" "both modes ran \e[0;32m[ok]\e[m"

    args+=("$dt" "input/${tag}_out.e" "input/${tag}_out_newton.e")
  done

  if [ $status -eq 0 ]; then
    echo
    python3 input/check_newton_order.py --min-ratio "$min_ratio" "${args[@]}"
    status=$?
  fi

  rm -f "${tmp_inputs[@]}"
  for dt in "${dts[@]}"; do
    rm -f input/_newtontest_${problem}_dt${dt}_out*.e \
          input/_newtontest_${problem}_dt${dt}_out*_convergence.csv
  done

  return $status
}

overall=0

# input/cavity.i IS the lid-driven cavity (lid = wall with -bc_5_u 1.0).
run_case cavity     0.2 1.5 0.01  0.005  0.0025  || overall=1
echo
run_case poiseuille 0.1 2.0 0.005 0.0025 0.00125 || overall=1

exit $overall
