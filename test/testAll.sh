#!/bin/bash
#WRAP TEST SCRIPT. EXPORT  OMP or CUDA for different set of tests
#export DEBUG for debug calls, and the target matrix to use in $1

mkdir -p log

MATRIX_DIR=~/data/prj/main/
matrixes=( $(find $MATRIX_DIR  -name "*mtx" | shuf) )
matrixesList=( $(cat "$1" | shuf ) )	#matrixes given in input file list
trgtMatrixes=${matrixes[@]}

set -e -o pipefail
if [ $CUDA ];then
	TEST_BIN=./test_SpMV_CUDA_Stats.o
	if [ $DEBUG ];then TEST_BIN=./test_SpMV_CUDA_DBG.o;	trgtMatrixes=${matrixesList[@]};fi
	for m in ${matrixes[@]};do
	    echo "#$m" ; ( $TEST_BIN $m RNDVECT | tee -a "log/testAllCUDA" ) || echo $m | tee -a log/testAllCUDA_errors
	done; echo echo "done CUDA";
elif [ $OMP ];then
	TEST_BIN=./test_SpMV_OMP_Stats.o;
	if [ $DEBUG ];then TEST_BIN=./test_SpMV_OMP_DBG.o;	trgtMatrixes=${matrixesList[@]};fi
	export GRID_ROWS=8 GRID_COLS=2
	for m in ${matrixes[@]};do
		echo "#$m" ; ( $TEST_BIN $m RNDVECT | tee "log/testAll"$GRID_ROWS"X"$GRID_COLS"_"$OMP_SCHEDULE ) || echo $m | tee -a log/testAllOMP_errors  
	done; echo echo "done $GRID_ROWS X $GRID_COLS"
else echo "EXPORT EITHER CUDA or OMP"
fi
