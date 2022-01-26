"""
parse output log over a large set of matrixes groups into a CSV on stdout
main log files structure: 
#matrixName
config
omp runtime config
@compute func ID ....
STATS... 
STATS... 
STATS... 
@compute func ID .... 
where STATS may be more then one for evaluating different amount of parallelism used (OMP set_thread_num )
see template files in this folder

usage <logFile>
"""
from collections import namedtuple
from re import finditer
from sys import argv,stderr
import string

_FIELDS     = "source,funcID,timeAvg,timeVar,internalTimeAvg,internalTimeVar,size,NNZ,maxRowNNZ,OMPgridSize,sampleSize"
_FIELDS_OPT = "threadNum"
FIELDS      = _FIELDS + " " + _FIELDS_OPT
Execution = namedtuple("Execution",FIELDS,defaults=["0"])

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

GRID_PATTERN="\d+x\d+"
SIZE_PATTERN=GRID_PATTERN+"-\d+NNZ-\d+=MAX_ROW_NZ"
FP_PATTERN="[-+]?\d+\.?\d+e[-+]\d+"
parseSizes=lambda s:  s #[int(x) for x in s.split("x")] #TODO not good CSV PARSED

def parseConfigSize(l):
    size  = parseSizes(getReMatch("sparse matrix:\s*("+SIZE_PATTERN+")",l))
    matSize,nnz,_maxRowNNZ = size.split("-")
    maxRowNNZ = _maxRowNNZ.split("=")[0]
    gridSize = parseSizes(getReMatch("grid:\s*("+GRID_PATTERN+")",l))
    sampleSize = getReMatch("AVG_TIMES_ITERATION:\s*(\d+)",l)
    return matSize,nnz,maxRowNNZ,gridSize,sampleSize

def parseComputeFuncID(l):
    funcID   = getReMatch("func:\s*(\d.*) at",l)
    return funcID

def parseComputeTimes(l):
    try:        threadNum = getReMatch("threadNum:\s*(\d+)",l)
    except:     threadNum = 0
    timeAvg = float(getReMatch("timeAvg:\s*("+FP_PATTERN+")",l))
    timeVar = float(getReMatch("timeVar:\s*("+FP_PATTERN+")",l))
    timeInternalAvg = float(getReMatch("timeInternalAvg:\s*("+FP_PATTERN+")",l))
    timeInternalVar = float(getReMatch("timeInternalVar:\s*("+FP_PATTERN+")",l))
    return timeAvg,timeVar,timeInternalAvg,timeInternalVar,threadNum

isClean = lambda l: len(l.strip())>1 and any(c in string.printable for c in l)
if __name__ == "__main__":
    if "-h" in argv[1] or len(argv)<2:  print(__doc__);exit(1)
    
    executionTimes = list() #Execution tups
    with open(argv[1]) as f:    log=f.read()
    matrixGroup = log.split("#")
    linesGroup = [ g.split("@") for g in matrixGroup ]
    
    for i,mGroup in enumerate(matrixGroup):
        if len(mGroup) < 3:  print("not complete mGroup",i,mGroup,file=stderr);continue
        mg = mGroup.split("\n")  
        header,configSiz,ompRuntimeSchedule = mg[0],mg[1],mg[2]
        src = header.replace(" ","_")

        matSize,nnz,maxRowNNZ,gridSize,sampleSize = parseConfigSize(configSiz)

        for compLinesFuncGroup in linesGroup[i][1:]: 
            #@computing SpMV with func:ID FUNC_ID_STR at:..... 
            compLinesFuncG = compLinesFuncGroup.split("\n")
            computesFuncID,computesTimes = compLinesFuncG[0],compLinesFuncG[1:]
            funcID = parseComputeFuncID(computesFuncID)
            for l in filter(lambda l:isClean,computesTimes):
                #-[threadNum=..],timeAvg,....
                try: timeAvg,timeVar,timeInternalAvg,timeInternalVar,threadNum = parseComputeTimes(l)
                except: continue;   #print("not clean line",l,file=stderr)  #TODO DEBUG
                executionTimes.append(Execution(src,funcID,timeAvg,timeVar,timeInternalAvg,\
                  timeInternalVar,matSize,nnz,maxRowNNZ,gridSize,sampleSize,threadNum))
    
    print(FIELDS)
    for e in executionTimes:    
        for f in e: print(f,end=", ")
        print("")
    print("\n")
