#!/bin/sh -e

set -x

QRY="${DATADIR}/small_reads.fq"
TGT="${DATADIR}/small_txes.fa"

$TRANSIGNER -o $RESULTS -t $TGT -p $THREADS $QRY