"""
parse output log over a large set of matrixes groups into a CSV on stdout
source matrixes lines starts with ##, computing and configuration lines start with @
expected prefixed lines in this order + template
##matrix
@sparse matrix: ROWSXCOLS-nnzNNZ grid: RXC
@compute ... sparse matrix: ROWSXCOLS-nnzNNZ with func:X timeAvg:XXX timeVar:XXX timeInternalAvg:XXX timeInternalVar:XXX 
TEMPLATE:
SpGEMV_OMP_test.c       AVG_TIMES_ITERATION:5, sparse matrix: 144649x144649-2148786NNZ  conf grid: 8x8 
ompGetAllICV:
omp sched gather:       kind: 2 - explicitly_NONmonotonic       omp chunkSize: 1
@computing SpGEMV with sparse matrix: 144649x144649-2148786NNZ  with func:0 at:0x403600      timeAvg:4.227846e-03 timeVar:3.075045e-08       timeInternalAvg:0.000000e+00 timeInternalVar:0.000000e+00 
...
usage <logFile>
"""
from collections import namedtuple
from re import finditer
from sys import argv,stderr

FIELDS = "source,funcN,timeAvg,timeVar,internalTimeAvg,internalTimeVar,size,NNZ,gridSize,sampleSize"
Execution = namedtuple("Execution",FIELDS)

getReGroups=lambda pattern,string:\
    finditer(pattern,string).__next__().groups()
getReMatch=lambda pattern,string:\
    finditer(pattern,string).__next__().groups()[0]
def getReMatchFull(pattern,string):
    #return first occurence of @pattern in @string or ""
    try:
        a=finditer(pattern,string).__next__()
        x,y=a.span()
        out=a.string[x:y]
    except: out=""
    return out

GRID_PATTERN="[0-9]+x[0-9]+"
SIZE_PATTERN=GRID_PATTERN+"-[0-9]+NNZ"
FP_PATTERN="[-+]?\d+\.?\d+e[-+]\d+"
parseSizes=lambda s:  s #[int(x) for x in s.split("x")] #TODO not good CSV PARSED

def parseConfigSize(l):
    size  = parseSizes(getReMatch("sparse matrix:\s*("+SIZE_PATTERN+")",l))
    matSize,nnz = size.split("-")
    gridSize = parseSizes(getReMatch("grid:\s*("+GRID_PATTERN+")",l))
    sampleSize = getReMatch("AVG_TIMES_ITERATION:\s*(\d)",l)
    return matSize,nnz,gridSize,sampleSize

def parseComputeTimes(l):
    funcN   = int(getReMatch("func:\s*(\d)",l))
    timeAvg = float(getReMatch("timeAvg:\s*("+FP_PATTERN+")",l))
    timeVar = float(getReMatch("timeVar:\s*("+FP_PATTERN+")",l))
    timeInternalAvg = float(getReMatch("timeInternalAvg:\s*("+FP_PATTERN+")",l))
    timeInternalVar = float(getReMatch("timeInternalVar:\s*("+FP_PATTERN+")",l))
    return funcN,timeAvg,timeVar,timeInternalAvg,timeInternalVar

if __name__ == "__main__":
    if "-h" in argv[1] or len(argv)<2:  print(__doc__);exit(1)
    
    executionTimes = list() #Execution tups
    with open(argv[1]) as f:    log=f.read()
    linesGroup = [ g.split("\n") for g in log.split("##")]
    
    for i,g in enumerate(linesGroup):
        if len(g) < 4:  print("not complete group",i,g,file=stderr);continue
        #splitting log parts of computation
        header,configSiz,runtimeSchedule= g[0],g[1],g[2]
        computes = list(filter(lambda l:"@" in l,g[2:]))
        #parsing
        src = header.replace(" ","_")
        matSize,nnz,gridSize,sampleSize = parseConfigSize(configSiz)
        #preparingTime = float(getReMatch("preparing time:\s*("+FP_PATTERN+")",configSiz))
        for l in computes:
            funcN,timeAvg,timeVar,timeInternalAvg,timeInternalVar = parseComputeTimes(l)
            executionTimes.append(Execution(src,funcN,timeAvg,timeVar,timeInternalAvg,timeInternalVar,matSize,nnz,gridSize,sampleSize))
    
    print(FIELDS)
    for e in executionTimes:    
        for f in e: print(f,end=", ")
        print("")
    print("\n")
