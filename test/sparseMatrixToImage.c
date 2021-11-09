/*
 * dev Andrea Di Iorio
 * draw a [sparse] matrix as an image with a black square dot for each nonzero elem
 * scaling the size of the matrix, "oring" the nz elem in a square of the original matrix into 
 * the smaller destination pixel grid.
 * Highlight nz elem into square of pixel of custom size, that will incorporate near elements
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <unistd.h>
#include <errno.h>
#include <sys/stat.h>

#include "macros.h"
#include "sparseMatrix.h"
#include "sparseMatrixToImage.h"
#include "utils.h"
#include "parser.h"

//inline export 
int IS_NNZ(spmat* smat,uint i,uint j);

///aux
double* parseDenseMatrixPatternFromFile(char* fpath,uint* M,uint* N){
    ERRPRINT("TODO IMPLEMENT");
    return NULL;
}
//set @i,@j in the scaled pixel grid into raw PPM -> RGB data as a black dot
static inline void _setBlackDotRGB(uint w,uint h, uchar rawData[][w],long i,long j){
    if  ( i<0 || i> (int)h )        return;
    for ( int jj=j; jj<j+PPM_CHAN_NUM; jj++ ){
        //skip outOfBorder enlarged dots
        if ( jj<0 || jj > (int) w)    continue; 
        //set RGB black -> no intensity in every channel
        rawData[i][jj] = NNZ_PIXEL_COLOR;
    }
}

/*
 * set a black dot pixel in MATRIX coordinate @i,@j , that will be enlarged
 * to be more visible to a square of @unifyNearW^2 near pixels 
 * (discarding previous 0 -> white dots, 
 *  NB never setted a white dot, memset at the start)
 *  TODO TODO FIX unifyNearW!!
 */ 
static void setBlackNZPixel(ppmData* data,long i,long j,ushort unifyNearW){
    //set the first dot
    _setBlackDotRGB(data->width*PPM_CHAN_NUM,data->height, 
      (uchar (*)[data->width*PPM_CHAN_NUM]) data->data,i,PPM_CHAN_NUM*(j));
    //set the enlarged dots
    for (short ww,w=0;w<unifyNearW;w++){
        for (short zz,z=0;z<unifyNearW;z++){
            //make the highlight unify square centered in (@i,@j)
            ww = INT_DIV_CEIL(w,2); 
            zz = INT_DIV_CEIL(z,2); 
            if (!(w % 2))   ww *= -1;
            if (!(z % 2))   zz *= -1;

            _setBlackDotRGB(data->width*PPM_CHAN_NUM,data->height, 
              (uchar (*)[data->width*PPM_CHAN_NUM]) data->data, 
              i+ww,PPM_CHAN_NUM*(j+zz));
        }
    }
}

/////TOPLEVEL IMAGE CONVERT FUNCTIONS
void denseMatrixToPPM(ppmData* data,
  uint M, uint N, double mat[][N],ushort step, ushort unifyNearW){
    char nz;    //flag to mark a founded nz
    for (uint i=0;i<M;i+=step){
        for (uint j=0;j<N;j+=step){
            //i,j point to the first: top,left element in the search square
            nz = 0;
            for (uint w=0; w<step && !nz; w++){
                for (uint z=0; z<step && !nz; z++){
                    if (mat[i+w][j+z]){
                        nz = (!0);
                        setBlackNZPixel(data,i/step,j/step,unifyNearW);
                        break;
                    }
                }
            }
        }
    }
}
void sparseMatrixToPPM(ppmData* data,spmat* sparseMat,
    ushort step, ushort unifyNearW){
    
    char nz;    //flag to mark a founded nz
    #pragma omp parallel for schedule(static) private(nz)
    for (uint i=0;      i<sparseMat->M; i+=step){
        for (uint j=0;  j<sparseMat->N; j+=step){
            //i,j point to the first: top,left element in the search square
            nz = 0;
            for (uint w=0; w<step && !nz; w++){
                for (uint z=0; z<step && !nz; z++){
                    if ((nz = IS_NNZ(sparseMat,i+w,j+z))){
                        setBlackNZPixel(data,i/step,j/step,unifyNearW);
                        break;
                    }
                }
            }
        }
    }
}

//check if MXN dense matrix @mat has a black pixel corresponding for each nonzero element
//TODO add support for step, unifyNearW, NOW CONSIDERED AS dflt... ->1a1 mat
static int checkDenseMatrixToPPM(uint M, uint N,
  unsigned char rawData[][3*N],double mat[][N],ushort step, ushort unifyNearW){
    for (uint i=0; i<N; i++){ 
        for (uint j=0; j<M; j++){
            if (mat[i][j] && (rawData[i][3*j+0] != NNZ_PIXEL_COLOR || 
              rawData[i][3*j+1] != NNZ_PIXEL_COLOR || rawData[i][3*j+2] != NNZ_PIXEL_COLOR)){
                fprintf(stderr,"not matching NNZ at (%u,%u)\n",i,j);
                return EXIT_FAILURE;
            }
            else if (!mat[i][j] && (rawData[i][3*j+0] != Z_PIXEL_COLOR || 
              rawData[i][3*j+1] != Z_PIXEL_COLOR || rawData[i][3*j+2] != Z_PIXEL_COLOR)){
                fprintf(stderr,"not matching NNZ at (%u,%u)\n",i,j);
                return EXIT_FAILURE;
            }
        }
    }
    return EXIT_SUCCESS;
}
#ifdef MAIN_SPMAT_IMG
/*
 * build this file as a stand alone executable, that get the input matrix from a serialized file
 * along with the parameters from argv defining MAT2PPM_STANDALONE
 * otherwise just export the next function
 */

#include <limits.h>
//inline export 
void freeSpmatInternal(spmat*);
void freeSpmat(spmat*);

int main(int argc,char** argv){
    void* map = NULL;
    int out=EXIT_FAILURE,outFd=0,mode=S_IRWXU;
    if (argc < 3 )    {
        ERRPRINT("usage: inputMatrixFile, "PATTERN_DENSE" || "MM_COO","
         " [collapseSquareW=1, unifyNearW=1, outPPMFPath="DFLT_OUTPATH"]\n");
        return EXIT_FAILURE;
    }
    double* denseMat=NULL;
    spmat*  sparseMat=NULL;
    ppmData* data=NULL;
    uint M,N;
    if (!strncmp(argv[2],PATTERN_DENSE,strlen(PATTERN_DENSE))){
        if(!(denseMat = parseDenseMatrixPatternFromFile(argv[1],&M,&N))){
            ERRPRINT("dense matrix pattern parse failed\n");
            return EXIT_FAILURE;
        }
    } else if (!strncmp(argv[2],MM_COO,strlen(MM_COO))){
        if (!(sparseMat = MMtoCSR(argv[1]))){
            ERRPRINT("sparse matrix MM coordinate parsing failed\n");
            return EXIT_FAILURE;
        }
        M = sparseMat -> M;
        N = sparseMat -> N;
        /*if (!(denseMat = CSRToDense(sparseMat))){ 
            ERRPRINT("CSRToDense FAILED\n");
            goto _free;
        }*/
    }else {ERRPRINT("INVALID IN MATRIX FORMAT!\n");return EXIT_FAILURE;}
    ///options
    uint collapseSquareW = 1, unifyNearW = 0;
    char *ptr;
    const char* outFname=DFLT_OUTPATH;
    if (argc > 3){
        collapseSquareW=strtoul(argv[3],&ptr,10);
        if (ptr==argv[3] ||  collapseSquareW== UINT_MAX){
            perror("strtol errd");
            return EXIT_FAILURE;
        }
    }
    if (argc > 4){
        unifyNearW=strtoul(argv[4],&ptr,10);
        if (ptr==argv[4] ||  unifyNearW== UINT_MAX){
            perror("strtol errd");
            return EXIT_FAILURE;
        }
    }
    if (argc > 5)   outFname = argv[5];

    data=malloc(sizeof(*data));
    if (!data){
        ERRPRINT("ppmData alloc for dense matrix failed\n");
        return EXIT_FAILURE;
    }
    //uint pad=2*pixelsPaddingPerElem;data- width=ceil(N/collapseSquareW)*(1+pad)
    //set out image size considering both scaling and padding
    data -> width  = ceil(N / collapseSquareW);
    data -> height = ceil(M / collapseSquareW);
    int headerLen = snprintf(data->header,PPM_HEADER_MAX_LEN,
        "P6\n%lu %lu\n255\n",data->width,data->height); 
    if (headerLen < 0){
        ERRPRINT("snprintf error");
        goto _free; 
    }
    ulong pixelDataLen,dataLen;
    if(__builtin_umull_overflow(data->width,data->height*PPM_CHAN_NUM,&pixelDataLen)){
        ERRPRINT("pixelDataLen overflow\n");
        goto _free;
    }
    if(__builtin_uaddl_overflow(pixelDataLen,headerLen,&dataLen)){
        ERRPRINT("dataLen overflow\n");
        goto _free;
    }
    printf("building ppm image at:%s with header:\n%s\n=>pixelDataLen: %luMB,"
      "collapseSquare:%u,unify:%u\n",outFname,data->header,dataLen>>20,collapseSquareW,unifyNearW);
    
    ///out file mmap for easy write
    outFd=open(outFname, O_RDWR | O_CREAT | O_EXCL | O_TRUNC, mode);
    if (errno==EEXIST)     outFd=open(outFname, O_RDWR | O_TRUNC, mode);
    if (outFd<0){
        perror("open outFd failed");
        goto _free;
    }
    if (ftruncate(outFd,dataLen)<0){
        perror("ftruncate err");
        goto _free;
    }
    map = mmap(NULL, dataLen, PROT_WRITE, MAP_SHARED, outFd, 0);    
    if (map == MAP_FAILED){
        perror("mmap failed... ");
        goto _free;
    }

    memcpy(map,data->header,headerLen);        //write header
    data->data=map+headerLen; //directly write converted matrix to outfile via mmap
    //memset(data->data,Z_PIXEL_COLOR,pixelDataLen); //ZERO=>IMPLICIT BY HOLE CREATED IN TRUNC
    printf("Alloc and prepare over, starting conversion\n");
    if (!strncmp(argv[2],MM_COO,strlen(MM_COO))){
    //TODO CONVERT OVER A DENSE [CONVERTED] MATRIX
    //denseMatrixToPPM(data,M,N,(double (*)[N]) denseMat,collapseSquareW,unifyNearW);
    sparseMatrixToPPM(data,sparseMat,collapseSquareW,unifyNearW);
    } //TODO else PATTERN_DENSE
    
#ifdef TEST
    if (collapseSquareW!=1 || unifyNearW !=0)
        {ERRPRINT("TODO MAKE TEST CASE FOR THIS CONFIG");goto _free;}
    if (!(denseMat = CSRToDense(sparseMat))){ 
        ERRPRINT("CSRToDense FAILED\n");
        goto _free;
    }
    if (checkDenseMatrixToPPM(M,N, (uchar (*)[N]) data->data,
        (double (*)[N]) denseMat,collapseSquareW,unifyNearW))   goto _free;
    //TODO add support for step, unifyNearW, NOW CONSIDERED AS dflt... ->1a1 mat
    printf("DENSE MATCHING TEST PASSEED\n");
#endif
    out = EXIT_SUCCESS;
    
    _free:
    if(outFd)       close(outFd);
    if(data)        free(data);
    if(denseMat)    free(denseMat);
    if(map){
        if(munmap(map,dataLen) == -1){
            perror("Error un-mmapping the file");
        }
    //if(ftruncate(outFd,actualLen)<0){perror("ftruncate err ");goto _free;}//remove excess from mmapped
    }
    if (sparseMat)  freeSpmat(sparseMat);
    return out;
}
#endif
