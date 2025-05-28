#!/bin/sh -e

set -x

QRY="${DATADIR}/small.bam"
N_READS=6160

$TRANSIGNER -o $RESULTS -p $THREADS --use-psw $QRY
"${SCRIPTS}/eval.py" -y "${DATADIR}/small_ground_truth.csv" \
    -y-hat "${RESULTS}" -n $N_READS > "${RESULTS}/eval.out"

set +x

awk '
BEGIN {
  # ANSI color codes
  GREEN = "\033[0;32m"
  RED = "\033[0;31m"
  NC = "\033[0m"

  # Expected thresholds
  expected["pearson'\''s"] = 0.999
  expected["spearman'\''s"] = 1.00
  expected["rmse"] = 9.250
  expected["precision"] = 0.986
  expected["recall"] = 0.985

  all_pass = 1
}
{
  metric = $1
  value = $2 + 0
  threshold = expected[metric]

  if (metric == "rmse") {
    if (value <= threshold) {
      status = GREEN "PASS" NC
    } else {
      status = RED "FAIL" NC
      all_pass = 0
    }
    comp = "<="
  } else {
    if (value >= threshold) {
      status = GREEN "PASS" NC
    } else {
      status = RED "FAIL" NC
      all_pass = 0
    }
    comp = ">="
  }

  printf "%-12s actual: %7.3f  expected: %7.3f (%s)  %s\n", metric, value, threshold, comp, status
}
END {
  if (all_pass)
    print "GOOD"
  else
    print "BAD"
}
' "${RESULTS}/eval.out" > "${RESULTS}.report"