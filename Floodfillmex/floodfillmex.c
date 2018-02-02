#include "mex.h"

    
void floodfillmex(double *map, double *y, double *z, double *path,mwSize m,mwSize n)
{
 int sx,sy,tx,ty,value=1,i,steps,v,j=0;
 int ss,tt;     // start & target 1-D array in C
 int *cc;    // pointer

 sx=y[1];       // start point
 sy=y[2];
 ss=(sy-1)*m+sx-1; // matlab address 2-D in C 1-D
 tx=z[1];     // target point
 ty=z[2];
 tt=(ty-1)*m+tx-1; // matlab address 2-D in C 1-D
/******************* C environment**********************/             
 /*************************map******************/ 
 /* matlab matrix[x][y]     C matrix [N]
  11 12 13 14 15 16      0 6  12 18 24 30 
  21 22 23 24 25 26      1 7  13 19 25 31
  31 32 33 34 35 36     2 8(S)14 20 26(T) 32
  41 42 43 44 45 46      3 9  15 21 27 33
  51 52 53 54 55 56      4 10 16 22 28 34 
  61 62 63 64 65 66      5 11 17 23 29 35
 */
 cc[0]=ss;
 //////////////// map initalize ///////////////////////////////////
 for (j;j>=0;j--){
     if (cc[j]-1>=0){
       if (map[cc[j]-1]==0){
           map[cc[j]-1]=value;    // check up value
       }
     }
     if ((cc[j]+1)<=(n*m-1)){
       if (map[cc[j]+1]==0){
           map[cc[j]+1]=value;   // check down value
        }      
     }
     if (cc[j]-m>=0){
       if (map[cc[j]-m]==0){
           map[cc[j]-m]=value;   // check left value
       }       
     }
     if (cc[j]+m<=(n*m-1)){ 
     if (map[cc[j]+m]==0){
           map[cc[j]+m]=value;   // check right value
       }       
     }
 }
     for (i=0;i++;i<=n*m-1){
         if (map[i]==value){
             cc[j]=i;
             j=j+1;
         }
     }
 value=value+1;
 /////////// map operation //////////////////////////////
 while (1){
 for (j;j>=0;j--){
     if (cc[j]-1>=0){
       if (map[cc[j]-1]==0){
           map[cc[j]-1]=value;    // check up value
       }
     }
     if ((cc[j]+1)<=(n*m-1)){
       if (map[cc[j]+1]==0){
           map[cc[j]+1]=value;   // check down value
        }      
     }
     if (cc[j]-m>=0){
       if (map[cc[j]-m]==0){
           map[cc[j]-m]=value;   // check left value
       }       
     }
     if (cc[j]+m<=(n*m-1)){ 
     if (map[cc[j]+m]==0){
           map[cc[j]+m]=value;   // check right value
       }       
     }
 }
 //////////////////get value position then fill map ////////////////////////
 j=0;
 for (i=0;i++;i<=n*m-1){
         if (map[i]==value){
             cc[j]=i;
             j=j+1;
         }
     }
  value=value+1;   
  if (value>n*m-1)// jump loop
      break;
 }
 /////////////////// path  //////////////////////////////////
 steps=map[tt];  // target nummber in map
 for (v=steps;v>0;v--) {
       if (tt-1>=0) {
        if (map[tt-1]==(steps-1)) {   // up
            path[steps-1]=1;            
        }
       } 
       if (tt+1<=n*m-1) {
        if (map[tt+1]==(steps-1)) {   // down
            path[steps-1]=2;
          }
       }
      if (tt-m>=0) {
          if (map[tt-m]==( steps-1)) {  // left
            path[steps-1]=3;
        }
      }
      if (ty-1<=n*m-1) {
          if (map[tt+m]==(steps-1)) {  // right
            path[steps-1]=4;
        }
      }        
 ///////find best number from end to start///////////////////////////// 
       // replace position then find again
       if (path[steps-1]==4) {
            tt=tt+m;
          }
          else if (path[steps-1]==3) {
            tt=tt-m;
          } 
          else if (path[steps-1]==2) {
            tt=tt+1;
          }
          else
            tt=tt-1;
     }
}

/*****************************mex interface function******************************************/  
/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
/* variable declarations here */
    double *inMatrix_map;               /* mxn input matrix */
    double *inMatrix_start;               /* 1x2 input matrix */
    double *inMatrix_target;               /* 1x2 input matrix */
    size_t mrows;                   /* size of matrix_map */
    size_t ncols;                   /* size of matrix_map */
    double *outMatrix;              /* output matrix about path */

/*****************************check I/O numbers 3/1 ******************************************/      
    /* check three inputs start_location, target_location and map, all array */ 
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
                          "Three inputs required.");
    }
    /* check one outputs path array */
    if(nlhs!=1) {
    mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs",
                      "One output required.");
    }
/*****************************check 3 input array type and 1 output type 2D/1D *****************/  
    if( !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0]) ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:not double",
                          "Input multiplier must be double.");
    }
    if( !mxIsDouble(prhs[1]) || 
         mxIsComplex(prhs[1]) ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:not double",
                          "Input multiplier must be double.");
    }
    if( !mxIsDouble(prhs[2]) || 
         mxIsComplex(prhs[2]) ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:not double",
                          "Input multiplier must be double.");
    }
   
    /* make sure the first input argument is 2D array */
        if((mxGetM(prhs[0])==1||
		   (mxGetN(prhs[0])==1)) ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector",
                          "Input must be a 2D vector.");
    }

    /* check that number of rows in second input argument is 1 */
    if(mxGetM(prhs[1])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector",
                          "Input must be a row vector.");
    }
    /* check that number of rows in third input argument is 1 */
    if(mxGetM(prhs[1])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector",
                          "Input must be a row vector.");
    }
 /***********************call the function floodfill from within mexfunction *****************/    
    
    /* get the value of the map array input  */
    inMatrix_map = mxGetPr(prhs[0]);

    /* create a pointer to the real data in the input matrix  */
    inMatrix_start = mxGetPr(prhs[1]);

    /* create second pointer to the real data in the second input matrix  */
    inMatrix_target = mxGetPr(prhs[2]);
    
    /* get dimensions of the input matrix */
    ncols = mxGetN(prhs[0]);
    mrows = mxGetM(prhs[0]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,mrows*ncols,mxREAL);

    /* get a pointer to the real data in the output matrix */
    outMatrix = mxGetPr(plhs[0]);

    /* call the computational routine */
    floodfillmex(inMatrix_map, inMatrix_start, inMatrix_target, outMatrix, (mwSize)mrows, (mwSize)ncols);       
/* code here */    
}        
