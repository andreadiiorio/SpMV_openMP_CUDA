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

#include "sparseMatrix"


#define PPM_HEADER_MAX_LEN 100
#define PPM_CHAN_NUM       3
#define RGB_BLACK          0xFFFFFF

#define DFLT_OUTPATH       "/tmp/mat.ppm"
typedef struct{
    char* data;     //width*height*3 triples of RGB pixel
    uint  width;

    uint  height;
    char[PPM_HEADER_MAX_LEN] header;
}   ppmData;

//set @i,@j in the scaled pixel grid into raw PPM -> RGB data as a black dot
static void _setBlackDotRGB(uint w,uint h, char rawData[][w],int i,int j){
    if (i<0 || j<0 || i >= w/3 || j > h/3)  return; //skip outOfBorder enlarged dots
    //set RGB black -> no intensity in every channel
    rawData[i][j]   = 0;
    rawData[i][j+1] = 0;
    rawData[i][j+1] = 0;
}

/*
 * set a black dot in pixel coordinate (top,left=0) @i,@j , that will be enlarged
 * to be more visible to a square of @unifyNearW^2 near pixels 
 * (discarding previous 0 -> white dots, 
 *  NB never setted a white dot, memset at the start)
 */ 
static void setBlackNZPixel(ppmData* data,int i,int j,ushort unifyNearW){
    for (short ww,w=0;w<unifyNearW;w++){
        for (short zz,z=0;z<unifyNearW;z++){
            //make the highlight unify square centered in (@i,@j)
            ww = INTDIVCEIL(w,2); 
            zz = INTDIVCEIL(z,2); 
            if (!(w % 2))   ww *= -1;
            if (!(z % 2))   zz *= -1;
            _setBlackDotRGB(data->with*CHAN_NUM, data->height, data->data, i+ww,j+zz);
        }
    }
}

/*
 * convert dense matrix in @data into ppm RGB triple pixels with a black dot per NZ elem in data->data
 * assign each consecutive @elemSquareWPixel^2 matrix elements 
 * to a dot in the PPM image, if at least 1 nz is in the square
 * the dot will also include unifyNearW^2 image pixel to let the dot more visible
 */
static void denseMatrixToPPM(ppmData* data, uint M, uint N, double mat[][N],
  ushort step, ushort unifyNearW){
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
/*
 * map each nz eleme in @sparseMat into a pixel dots square
 * similarly as before in _denseMatrixToPPM
 */
static void sparseMatrixToPPM(ppmData* data,spmat* sparseMat,
    ushort step, ushort unifyNearW){
    //TODO 
}
/*
 * TODO MACRO SELECTION
 * build this file as a stand alone executable, that get the input matrix from a serialized file
 * along with the parameters from argv defining MAT2PPM_STANDALONE
 * otherwise just export the next function
 */
int main(int argc,char** argv){
    void* map = NULL;
    int outFd=0,mode=S_IRWXU;
    //int CHAN_NUM = 4 ;                 //default .pam with RBG_ALPHA
    if (argc < 3 )    {
        ERRPRINT("usage: inputMatrixFile, "PATTERN_DENSE" || "MM_COO","
         " [elemSquareWPixel=1, unifyNearW=1, outPPMFPath="DFLT_OUTPATH"]\n");
        return EXIT_FAILURE;
    }
    double* denseMat;
    spmat*  sparseMat;
    if (strncmp(argv[2],DENSE,strlen(PATTERN_DENSE))){
        if(!(denseMat = parseDenseMatrixPatternFromFile(argv[1]))){
            ERRPRINT("dense matrix pattern parse failed\n");
            return EXIT_FAILURE;
        }
    } else if (strncmp(argv[3],MM_COO,strlen(MM_COO))){
        if (!(sparseMat = MMtoCSR(argv[1]))){
            ERRPRINT("sparse matrix MM coordinate parsing failed\n");
            return EXIT_FAILURE;
        }
        if (!(denseMat = CSRToDense(sparseMat))){ //TODO REPLACE WITH _sparseMatrixToPPM LATER
            ERRPRINT("CSRToDense FAILED\n");
            goto _free;
        }
    }
    } else {ERRPRINT("INVALID IN MATRIX FORMAT!\n");return EXIT_FAILURE;}
        
    const char* outFname=DFLT_OUTPATH;
    if (argc >= 4)    outFname = argv[3];

    ppmData* data=malloc(sizeof(*data));
    if (!data){
        ERRPRINT("ppmData alloc for dense matrix failed\n");
        return NULL;
    }
    //uint pad=2*pixelsPaddingPerElem;data- width=ceil(N/elemSquareWPixel)*(1+pad)
    //set out image size considering both scaling and padding
    data -> width  = ceil(N / elemSquareWPixel);
    data -> height = ceil(M / elemSquareWPixel);
    int headerLen = snprintf(data->header,PPM_HEADER_MAX_LEN,
        "P6\n%u %u\n255\n",data->width,data->header); 
    if (headerLen < 0){
        ERRPRINT("snprintf error");
        goto _free; 
    }
    uint dataLen = data->width * data->height * PPM_CHAN_NUM + headerLen;
    
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
    map = mmap(NULL, outFileMaxLen, PROT_WRITE, MAP_SHARED, outFd, 0);    
    if (map == MAP_FAILED){
        perror("mmap failed... ");
        goto _free;
    }


    memcpy(map,data->header,headerLen);        //write header
    memset(map+headerLen,255,pixelDataLen);    //preset the img as white
    data->data=map; //directly write converted matrix to outfile via mmap
    denseMatrixToPPM(M,N,denseMat,elemSquareWPixel);//TODO REPLACE WITH _sparseMatrixToPPM LATER


    if(munmap(map,outFileMaxLen) == -1){
        perror("Error un-mmapping the file");
        goto _free;
    }
    //if(ftruncate(outFd,actualLen)<0){perror("ftruncate err ");goto _free;}//remove excess from mmapped

    out = EXIT_SUCCESS;
    _free;
    if(outFd)   close(outFd);
    free(data);
    return out;
}
