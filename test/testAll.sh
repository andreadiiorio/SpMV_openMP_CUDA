#!/bin/bash
set -e -o pipefail

MATRIX_DIR=~/data/prj/main/
MATRIX_LIST=$1	#log/testAllCUDA_errors

mkdir -p log
#CUDA
#TEST_BIN=./test_SpMV_CUDA_Stats.o
TEST_BIN=./test_SpMV_CUDA_DBG.o

if [ $CUDA ];then
	#for m in $( find $MATRIX_DIR  -name "*mtx" | shuf);do 
	for m in $( cat $MATRIX_LIST | shuf );do
	    echo "#$m" ; ( $TEST_BIN $m RNDVECT | tee -a "log/testAllCUDA" ) || echo $m | tee -a log/testAllCUDA_errors
	done; echo echo "done CUDA";
elif [ $OMP ];then
	#OMP
	TEST_BIN=./test_SpMV_OMP_DBG.o;
	#TEST_BIN=./test_SpMV_OMP_Stats;
	export GRID_ROWS=8 GRID_COLS=2
	#for m in $( find $MATRIX_DIR  -name "*mtx.gz" | shuf);do
	for m in $( cat $MATRIX_LIST | shuf );do
		echo "#$m" ; ( $TEST_BIN $m RNDVECT | tee "log/testAll"$GRID_ROWS"X"$GRID_COLS"_"$OMP_SCHEDULE ) || echo $m | tee -a log/testAllOMP_errors  
	done; echo echo "done $GRID_ROWS X $GRID_COLS"
else echo "EXPORT EITHER CUDA or OMP"
fi
