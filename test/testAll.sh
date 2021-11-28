#!/bin/bash
set -e -o pipefail

MATRIX_DIR=~/Documenti/data/prj/
TEST_BIN=./test_SpGEMV_OMP.o
export GRID_ROWS=8 GRID_COLS=2; for m in $( find $MATRIX_DIR  -name "*mtx.gz" | shuf);do echo "#$m" ; ( $TEST_BIN $m RNDVECT | tee "log/testAll"$GRID_ROWS"X"$GRID_COLS"_"$OMP_SCHEDULE ) || exit 22 ;done; echo echo "done $GRID_ROWS X $GRID_COLS"
