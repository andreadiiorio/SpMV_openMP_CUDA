#!/bin/bash
#WRAP TEST SCRIPT. EXPORT  OMP or CUDA for different set of tests
#export DEBUG for debug calls, and the target matrix to use in $1

mkdir -p log

MATRIX_DIR=~/data/prj/main/
matrixes=( $(find $MATRIX_DIR  -name "*mtx" ) )
matrixesList=( $(cat "$1" ) )	#matrixes given in input file list
trgtMatrixes=${matrixes[@]}

set -e -o pipefail
if [ $CUDA ];then
	TEST_BIN=./test_SpMV_CUDA_Stats.o
	if [ $DEBUG ];then TEST_BIN=./test_SpMV_CUDA_DBG.o;	trgtMatrixes=${matrixesList[@]};fi
	for m in ${matrixes[@]};do
	    echo "#$m" ; ( $TEST_BIN $m RNDVECT | tee -a "log/testAllCUDA.log" ) || echo $m | tee -a log/testAllCUDA_errors
	done; echo echo "done CUDA";
fi
if [ $OMP ];then
	TEST_BIN=./test_SpMV_OMP_Stats_DECREASING_T_NUM.o
	if [ $DEBUG ];then TEST_BIN=./test_SpMV_OMP_DBG.o;	trgtMatrixes=${matrixesList[@]};fi

	export GRID_ROWS=8 GRID_COLS=5;ompConf="$GRID_ROWS"X"$GRID_COLS"_"$OMP_SCHEDULE"
	for m in ${matrixes[@]};do		echo "#$m" ; ( $TEST_BIN $m RNDVECT | tee -a "log/testAll_$ompConf.log" ) || echo $m | tee -a log/testAllOMP_"$ompConf"_errors; done; echo echo "done $GRID_ROWS X $GRID_COLS"
	if [ $DEBUG ];then exit $?;fi
	export GRID_ROWS=5 GRID_COLS=8;ompConf="$GRID_ROWS"X"$GRID_COLS"_"$OMP_SCHEDULE"
	for m in ${matrixes[@]};do		echo "#$m" ; ( $TEST_BIN $m RNDVECT | tee -a "log/testAll_$ompConf.log" ) || echo $m | tee -a log/testAllOMP_"$ompConf"_errors; done; echo echo "done $GRID_ROWS X $GRID_COLS"
	export GRID_ROWS=10 GRID_COLS=4;ompConf="$GRID_ROWS"X"$GRID_COLS"_"$OMP_SCHEDULE"
	for m in ${matrixes[@]};do		echo "#$m" ; ( $TEST_BIN $m RNDVECT | tee -a "log/testAll_$ompConf.log" ) || echo $m | tee -a log/testAllOMP_"$ompConf"_errors; done; echo echo "done $GRID_ROWS X $GRID_COLS"
	export GRID_ROWS=4 GRID_COLS=10;ompConf="$GRID_ROWS"X"$GRID_COLS"_"$OMP_SCHEDULE"
	for m in ${matrixes[@]};do		echo "#$m" ; ( $TEST_BIN $m RNDVECT | tee -a "log/testAll_$ompConf.log" ) || echo $m | tee -a log/testAllOMP_"$ompConf"_errors; done; echo echo "done $GRID_ROWS X $GRID_COLS"
	export GRID_ROWS=14 GRID_COLS=3;ompConf="$GRID_ROWS"X"$GRID_COLS"_"$OMP_SCHEDULE"
	for m in ${matrixes[@]};do		echo "#$m" ; ( $TEST_BIN $m RNDVECT | tee -a "log/testAll_$ompConf.log" ) || echo $m | tee -a log/testAllOMP_"$ompConf"_errors; done; echo echo "done $GRID_ROWS X $GRID_COLS"
	export GRID_ROWS=13 GRID_COLS=3;ompConf="$GRID_ROWS"X"$GRID_COLS"_"$OMP_SCHEDULE"
	for m in ${matrixes[@]};do		echo "#$m" ; ( $TEST_BIN $m RNDVECT | tee -a "log/testAll_$ompConf.log" ) || echo $m | tee -a log/testAllOMP_"$ompConf"_errors; done; echo echo "done $GRID_ROWS X $GRID_COLS"
#else echo "EXPORT EITHER CUDA or OMP"
fi
