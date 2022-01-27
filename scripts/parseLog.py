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

_FIELDS     = "source,funcID,timeAvg,timeVar,internalTimeAvg,internalTimeVar,matRows,matCols,NNZ,maxRowNNZ,sampleSize"
_FIELDS_OPT_OMP  = "ompSchedKind,ompChunkSize,ompMonotonic,threadNum,ompGridSize_x,ompGridSize_y"
_FIELDS_OPT_CUDA = "blockSize_x,blockSize_y,blockSize_z,gridSize_x,gridSize_y,gridSize_z"
_FIELDS_OPT      = _FIELDS_OPT_OMP + "," + _FIELDS_OPT_CUDA
FIELDS      = _FIELDS + "," + _FIELDS_OPT
Execution = namedtuple("Execution",FIELDS,defaults=[None]*len(_FIELDS_OPT.split(",")) )

hasFields = lambda l: len(l.strip().split()) > 2
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
parseGridPattern = lambda s: s.strip().split("x")
FP_PATTERN="[-+]?\d+\.?\d+e[-+]\d+"
parseSizes=lambda s:  s #[int(x) for x in s.split("x")] #TODO not good CSV PARSED

def parseConfigSize(l):
    size  = parseSizes(getReMatch("sparse matrix:\s*("+SIZE_PATTERN+")",l))
    matSizes,nnz,_maxRowNNZ = size.split("-")
    matSizes  = parseGridPattern(matSizes)
    maxRowNNZ = int(_maxRowNNZ.split("=")[0])
    nnz = int( nnz.replace("NNZ","") )
    sampleSize = getReMatch("AVG_TIMES_ITERATION:\s*(\d+)",l)
    return matSizes,nnz,maxRowNNZ,sampleSize

def parseComputeFuncID(l):
    funcID   = getReMatch("func:\s*(.*) at",l)
    return funcID

def parseCompute(l):
    try:        
            threadNum = getReMatch("threadNum:\s*(\d+)",l)
            ompGridSize = parseSizes(getReMatch("ompGridSize:\s*("+GRID_PATTERN+")",l))
            ompGridSize = parseGridPattern(ompGridSize)
    except: threadNum,ompGridSize=0,[None,None]     #expected 2D  ompGridSize
    try:        
            cudaBlkSize   = getReMatch("cudaBlockSize:\s*(\d+\s+\d+\s+\d+)",l).split()
            cudaGridSize  = getReMatch("cudaGridSize:\s*(\d+\s+\d+\s+\d+)",l).split()
    except: cudaBlkSize,cudaGridSize = [],[]
    timeAvg = float(getReMatch("timeAvg:\s*("+FP_PATTERN+")",l))
    timeVar = float(getReMatch("timeVar:\s*("+FP_PATTERN+")",l))
    timeInternalAvg = float(getReMatch("timeInternalAvg:\s*("+FP_PATTERN+")",l))
    timeInternalVar = float(getReMatch("timeInternalVar:\s*("+FP_PATTERN+")",l))
    return timeAvg,timeVar,timeInternalAvg,timeInternalVar,\
      threadNum,ompGridSize,cudaBlkSize,cudaGridSize

def parseOmpRuntimeSchedule(l):
    schedKind = getReMatch("kind:\s*(OMP_.+)\s+omp chunk",l)
    chunkSize = getReMatch("omp chunkSize:\s+(\d+)",l)
    monotonic = getReMatch("monotonic:\s(.)",l)
    return schedKind,chunkSize,monotonic
if __name__ == "__main__":
    if "-h" in argv[1] or len(argv)<2:  print(__doc__);exit(1)
    
    executionTimes = list() #Execution tups
    with open(argv[1]) as f:    log=f.read()
    matrixGroup = log.split("#")
    linesGroup = [ g.split("@") for g in matrixGroup ]
    
    for i,mGroup in enumerate(matrixGroup):
        if len(mGroup) < 3:  print("not complete mGroup",i,mGroup,file=stderr);continue
        mg = mGroup.split("\n")  
        header,configSiz,ompSched= mg[0],mg[1],mg[2]
        src = header.replace(" ","_")
        ompSched,ompChunkSize,ompMonotonic = parseOmpRuntimeSchedule(ompSched)
        _ompSchedCfg = [ ompSched, ompChunkSize, ompMonotonic ]
        matSizes,nnz,maxRowNNZ,sampleSize = parseConfigSize(configSiz)

        for compLinesFuncGroup in linesGroup[i][1:]: 
            #@computing SpMV with func:ID FUNC_ID_STR at:..... 
            compLinesFuncG = compLinesFuncGroup.split("\n")
            computesFuncID,computesTimes = compLinesFuncG[0],compLinesFuncG[1:]
            funcID = parseComputeFuncID(computesFuncID)
            for l in filter(hasFields,computesTimes):
                #-[threadNum=..],timeAvg,....
                try: tAvg,tVar,tIntAvg,tIntVar,threadN,ompGridSize,cudaBlkSize,cudaGridSize=parseCompute(l)
                except Exception as e: 
                    print("not clean line",l,e,file=stderr);continue

                #obfuscate default prints for easier reading tables
                isCudaEntry      = None in ompGridSize
                if isCudaEntry:
                    _ompSchedCfg = [None] * len(_ompSchedCfg)
                    threadN      = None
                #else:  #OMP ENTRY, ... 

                #insert parsed compute entries infos
                executionTimes.append(Execution(src,funcID,tAvg,tVar,tIntAvg,\
                  tIntVar,*matSizes,nnz,maxRowNNZ,sampleSize,\
                  *_ompSchedCfg,threadN,*ompGridSize,*cudaBlkSize,*cudaGridSize))
    
    print(FIELDS)
    for e in executionTimes:    
        for f in e: print(f,end=", ")
        print("")
    print("\n")
