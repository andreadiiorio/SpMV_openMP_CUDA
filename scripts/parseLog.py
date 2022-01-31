#Copyright Andrea Di Iorio 2022
#This file is part of SpMV_OMP_CUDA
#SpMV_OMP_CUDA is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#SpMV_OMP_CUDA is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with SpMV_OMP_CUDA.  If not, see <http://www.gnu.org/licenses/>.

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
==============================================================================================================
export: GROUP_IMPLEMENTATIONS=[false] -> group several compute entries of the same funcID and ComputeConf
                                         adding (in order) implementations compute times as new csvFields like
                                         timeAvg_funcID0_Conf0, timeAvg_funcID0_Conf1, ... timeAvg_funcID1_Conf0, ...
                                        
                                         PARSED LOG CAN HAVE LESS NUM OF COMPUTE LINES, 
                                         BUT THE SAME INITIAL funcID,ComputeConf lines (smaller groups are padded)
        GROUP_IMPLEMENTATIONS_KEEP_CONST_CCONF: assuming ComputeConf being costant across the imput
                                                ComputeConf fields will not be replaced by GROUPD Macro
        FLOAT_PRECISION_PY=[e.g. .17e -> precision of double to output in the csv]
        
usage <logFile>
"""
from collections import namedtuple
from re          import finditer
from sys         import argv,stderr
from os          import environ as env

GROUP_IMPLEMENTATIONS= "T" in env.get("GROUP_IMPLEMENTATIONS","F").upper()
GROUP_IMPLEMENTATIONS_KEEP_CONST_CCONF = "T" in env.get("GROUP_IMPLEMENTATIONS_KEEP_CONST_CCONF","F").upper()
FLOAT_PRECISION_PY=env.get("FLOAT_PRECISION_PY",".17e")

_FIELDS_MAIN     = "source,funcID,timeAvg,timeVar,internalTimeAvg,internalTimeVar,matRows,matCols,NNZ,maxRowNNZ,sampleSize"
MAIN_FIELDS = _FIELDS_MAIN.split(",")
#compute config for iterative run as optional fields (requires new python)
_FIELDS_OPT_OMP  = "ompSchedKind,ompChunkSize,ompMonotonic,threadNum,ompGrid"
_FIELDS_OPT_CUDA = "blockSize_x,blockSize_y,blockSize_z,gridSize_x,gridSize_y,gridSize_z"
_FIELDS_OPT      = _FIELDS_OPT_OMP + "," + _FIELDS_OPT_CUDA
FIELDS           = _FIELDS_MAIN + "," + _FIELDS_OPT
Execution   = namedtuple("Execution",FIELDS) #,defaults=[None]*len(_FIELDS_OPT.split(",")) ) #require new python
ComputeConf = namedtuple("ComputeConf",_FIELDS_OPT)

GROUPD = "GRPD_ENTRY"   #entry that has been groupped (e.g. because of GROUP_IMPLEMENTATIONS=T)
PADD   = "PADD" #None
GROUP_IMPLEMENTATIONS_TRGT_FIELDS = ["timeAvg"] #,"timeVar"] #fields to "multiplex"
#aux ComputeConf fields selection for GROUP_IMPLEMENTATIONS 
_none           = lambda l: []
_identity       = lambda l: l
_ompGrid        = lambda l: [l[4]]
selectFieldsToAdd = _none

filterCompLines = lambda l: "OMP CSR 0" not in l

hasFields = lambda l,fNum=2: len(l.strip().split()) > fNum or len(l.strip().split("/")) > fNum
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


def parseOmpRuntimeSchedule(l):
    schedKind = getReMatch("kind:\s*(OMP_.+)\s+omp chunk",l)
    chunkSize = getReMatch("omp chunkSize:\s+(\d+)",l)
    monotonic = getReMatch("monotonic:\s(.)",l)
    return schedKind,chunkSize,monotonic

###main parse from here
def parseComputeL(l):
    try:        
            threadNum = getReMatch("threadNum:\s*(\d+)",l)
            ompGridSize = parseSizes(getReMatch("ompGridSize:\s*("+GRID_PATTERN+")",l))
            ompGridSize = parseGridPattern(ompGridSize)
    except: threadNum,ompGridSize=0,[None,None]     #expected 2D  ompGridSize
    try:        
            cudaBlkSize   = getReMatch("cudaBlockSize:\s*(\d+\s+\d+\s+\d+)",l).split()
            cudaGridSize  = getReMatch("cudaGridSize:\s*(\d+\s+\d+\s+\d+)",l).split()
    except: cudaBlkSize,cudaGridSize = [None]*3,[None]*3
    timeAvg = float(getReMatch("timeAvg:\s*("+FP_PATTERN+")",l))
    timeVar = float(getReMatch("timeVar:\s*("+FP_PATTERN+")",l))
    timeInternalAvg = float(getReMatch("timeInternalAvg:\s*("+FP_PATTERN+")",l))
    timeInternalVar = float(getReMatch("timeInternalVar:\s*("+FP_PATTERN+")",l))
    return timeAvg,timeVar,timeInternalAvg,timeInternalVar,\
      threadNum,ompGridSize,cudaBlkSize,cudaGridSize

def parseComputeLines(lGroup,configSiz,ompSched):
    """
    parse set of computation lines groups in @lGroup
    composed of a series of (at least 2) lines, an headerLine and compute times
    @computing SpMV with func:ID FUNC_ID_STR at:.....   (header)
    [threadNum=..,cudaBlkSize,...],timeAvg,....         (compute lines)
    [threadNum=..,cudaBlkSize,...],timeAvg,....
    
    @Returns:   list of Execution namedtuple of parsed lines @lGroup
    """
    out = list()
    matSizes,nnz,maxRowNNZ,sampleSize = parseConfigSize(configSiz)
    ompSched,ompChunkSize,ompMonotonic = parseOmpRuntimeSchedule(ompSched)
    ompSchedCfg = [ ompSched, ompChunkSize, ompMonotonic ]

    for compLinesFuncGroup in filter(filterCompLines,lGroup):
        #@computing SpMV with func:ID FUNC_ID_STR at:..... 
        compLinesFuncG = compLinesFuncGroup.split("\n")
        computesFuncID,computesTimes = compLinesFuncG[0],compLinesFuncG[1:]
        funcID = parseComputeFuncID(computesFuncID)
        for l in filter(hasFields,computesTimes):
            #-[threadNum=..],timeAvg,....
            try: tAvg,tVar,tIntAvg,tIntVar,threadN,ompGridSize,cudaBlkSize,cudaGridSize=parseComputeL(l)
            except Exception as e: 
                print("not clean line",l,e,file=stderr);continue

            #obfuscate default prints for easier reading tables
            isCudaEntry      = None in ompGridSize
            if isCudaEntry:
                ompSchedCfg    = [None] * len(ompSchedCfg)
                threadN         = None
                tIntAvg,tIntVar = None,None #TODO mesured kernel time only
            #else:  #OMP ENTRY, ... 

            #insert parsed compute entries infos
            out.append(Execution(src,funcID,tAvg,tVar,tIntAvg,\
              tIntVar,*matSizes,nnz,maxRowNNZ,sampleSize,\
              *ompSchedCfg,threadN,"x".join(ompGridSize),*cudaBlkSize,*cudaGridSize))
    return out
def groupImplementations(executionTimes,trgtFields=GROUP_IMPLEMENTATIONS_TRGT_FIELDS):
    """
       given the Execution named tuples in @executionTimes
       group them in a single entry, keeping the info of the funcID and computeConf
       @trgtFields: list of times fields to "multiplex" for each groupped field
       (group by every field not in trgtFields, in particular (optional) computeConf fields)

       @grupdCConfFieldsSelect: select what fields preserve from computeConf of each
         entry in @executionTimes 

       :Returns an Execution namedtuple, where the given @trgtFields will be
        lists of ( funcID,ComputeConf,time )
    """
    trgtFieldsGroups = { trgtF:list() for trgtF in trgtFields }
    mainFixdFields = executionTimes[0][:len(MAIN_FIELDS)]
    optCConfFieldsFixFirst = executionTimes[0][len(MAIN_FIELDS):]
    for e in executionTimes:
        groupdFields = ( e.funcID, ComputeConf(*e[len(MAIN_FIELDS):]) )
        #gather trgtFields of each compute, along with its context info
        for trgtF in trgtFields:    
            trgtFieldsGroups[trgtF].append( [*groupdFields, getattr(e,trgtF)] )
    #get the output namedtuple as a "blank" entry, setting the main fields of the first entry
    #andy a the constant GROUPD for every other fields, that has been groupped in the out entry
    optCConfFields = [GROUPD]*len(_FIELDS_OPT.split(","))
    if GROUP_IMPLEMENTATIONS_KEEP_CONST_CCONF: optCConfFields = optCConfFieldsFixFirst
    out = Execution(*mainFixdFields,*optCConfFields)
    out = out._replace(funcID=GROUPD)
    #set the target,groupped fields in the output entry
    for trgtF in trgtFields:    out = out._replace( **{trgtF:trgtFieldsGroups[trgtF]} )
    
    return out

if __name__ == "__main__":
    if len(argv) < 2 or "-h" in argv[1]:  print(__doc__);exit(1)
    
    executionTimes = list() #Execution tups
    with open(argv[1]) as f:    log=f.read()
    matrixGroup = log.split("#")
    linesGroup = [ g.split("@") for g in matrixGroup ]
    
    for i,mGroup in enumerate(matrixGroup):
        if len(mGroup) < 3:  print("not complete mGroup",i,mGroup,file=stderr);continue
        mg = mGroup.split("\n")  
        header,configSiz,ompSched= mg[0],mg[1],mg[2]
        src = header.replace(" ","_").split("/")[-1]
        computeEntries = parseComputeLines(linesGroup[i][1:],configSiz,ompSched)
        if GROUP_IMPLEMENTATIONS: #merge all computeEntries in a single one
            computeEntries = [ groupImplementations(computeEntries) ]
        executionTimes += computeEntries
   
    #audit data as CSV
    if GROUP_IMPLEMENTATIONS:
        ##RBcomputeEntries have for each trgtFields lists of ( funcID,ComputeConf,time )

        #get trgtFields (one of them) context informations for CSV header
        _trgtF_0 = GROUP_IMPLEMENTATIONS_TRGT_FIELDS[0]
        _largestGrppdEntryLen,largestGrppdEntryIdx \
            = max((len(getattr(e,_trgtF_0)),i) for i,e in enumerate(executionTimes))
        largestGrppdEntry = executionTimes[largestGrppdEntryIdx]
        csvMultiplexedFiledsSufx = list()
        for funcID,computeConf,_t in getattr(largestGrppdEntry,_trgtF_0):
            cconfCSVfields = selectFieldsToAdd(computeConf)
            csvMultiplexedFiledsSufx.append("_".join([str(x) for x in (funcID,*cconfCSVfields)]))
        multiplexedCSVHeader = str()
        for f in MAIN_FIELDS:
            if f in GROUP_IMPLEMENTATIONS_TRGT_FIELDS:
                f = ", ".join([f+"_"+suffx for suffx in csvMultiplexedFiledsSufx ])
            multiplexedCSVHeader    += f+", "
        print(multiplexedCSVHeader[:-2],_FIELDS_OPT,sep=",   ")    #remove last ","
        #get a dummy entry to pad Execution entries with target fields with less (multiplexed) values then the max
        #padEntryV = ((funcID,cconf,None) for funcID,cconf,_t in getattr(largestGrppdEntry,_trgtF_0))
        padEntryV = getattr(largestGrppdEntry,_trgtF_0)
        for i in range(len(padEntryV)): padEntryV[i][-1] = PADD
        #dump csv rows, padding entries with less multiplexed entries
        for e in executionTimes:    
            for f,x in e._asdict().items():
                if f in GROUP_IMPLEMENTATIONS_TRGT_FIELDS:
                    toPadN = _largestGrppdEntryLen - len(x)
                    if toPadN > 0:  x += padEntryV[len(x):]
                    #here select only the times of the groupped(context-ed) values
                    x = ",".join([str(xx[-1]) for xx in x])
                print(x,end=", ")
                #TODO TODO
            print("")

    else:       #not GROUP_IMPLEMENTATIONS
        print(FIELDS)
        for e in executionTimes:       #dump csv rows
            for f in e: print(f,end=", ")
            print("")
            #if type(f) == float:    print(format(f,FLOAT_PRECISION_PY),end=", ")
            #else:                   print(f,end=", ")
        print("\n")
