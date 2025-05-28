#!/bin/sh

SCRATCH=$1
TRANSIGNER=$2
THREADS=$3

mkdir -p $SCRATCH

export SCRIPTS="./test_scripts"
export DATADIR="./test_data"
export THREADS
export TRANSIGNER

fn="${DATADIR}/big.bam"
url="ftp://ftp.ccb.jhu.edu/pub/data/NA12878_simulated_ONT/dRNA_simulated_1.bam"
if [ ! -e "$fn" ]; then
    echo "$fn not found, downloading..."
    cd $DATADIR || { echo "failed to cd into $DATADIR"; exit 1; }
    wget $url -O "big.bam" || { echo "download failed"; exit 1; }
fi

TESTS=""
run_test() {
    if [ "$#" -lt 2 ]; then
        echo "A test requires >= 2 params! exiting..."
        exit 1
    fi
    NAME="$1"
    FILE="$2"
    shift
    shift
    TESTS="${TESTS} ${NAME}"
    export RESULTS="${SCRATCH}/${NAME}"
    mkdir -p $RESULTS
    START="$(date +%s)"
    "${SCRIPTS}/${FILE}" "$@"
    STATUS="$?"
    END="$(date +%s)"
    if [ "${STATUS}" = "0" ]; then
        if [ -f "${RESULTS}.report" ] && [ "$(echo $(tail -n 1 "${RESULTS}.report"))" = "GOOD" ]; then
            rm -rf "${RESULTS}"
        fi
    fi
    eval "${NAME}_TIME"="$((END-START))"
}

set +e

run_test SMALL_DEFAULT "run_small_default.sh"
run_test SMALL_PSW "run_small_psw.sh"

# TODO: add BIG test cases

set -e

printf "\n"
ERR=0
for i in ${TESTS}; do
    VAL="${i}_TIME"
    eval TIME="\$$VAL"
    printf "\033[1m$i (Time: %ss)\033[0m\n" "${TIME}"
    STATUS="$(tail -n 1 "${SCRATCH}/${i}.report")"
    if [ "$STATUS" != "GOOD" ]; then
        printf "\033[31mTEST FAILED\033[0m\n"
        ERR=$((ERR+1))
    else
        printf "\033[32mTEST SUCCESS\033[0m\n"
    fi
    head -n -1 "${SCRATCH}/${i}.report"
    printf "\n"
done