/***********************************************************************************************
 *
 * MATLAB:  mem2 = Detection2(A);
 * 
 * 
 *
 *  Stefan Schulte (stefan.schulte@Ugent.be)
 *  Last modified: 30/05/06
 *
 ************************************************************************************************/

#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
//   printf("\t%d %s %f %s %c\n", x, str, pi, "WOW", c);


double LARGE (double x, double p1, double p2) {
   double res = 0.0;
   if ((x > p1) && (x < p2))   res = (x-p1)/(p2-p1);
   else if (x > p2)            res = 1.0;
   else                        res = 0.0;
   return res;
}

double absol(double a) {
   double b;
   if(a<0)   b=-a;
   else   b=a;
   return b;
}

double minimum(double a, double b) {
   double x;
   if(a<=b)   x = a;
   else   x=b;
   return x;
}

double maximum(double a, double b) {
   double x;
   if(a<=b)   x = b;
   else   x=a;
   return x;
}

/**************************************************************************************
*  The main function for the calculation of the shrinkage method
*
***************************************************************************************/

void CALL_MEM(double **A1,double W,double W2,double *mem,int M, int N) { 
   int i,j,k,l,pos,Ymax,ind;   
   double tmp, tmp2,hlp,p1,p2,m1,m2,Xmax;
   int *H1, *H2, *H3;
   
   int rand1a = 0;
   int rand1b = 0;
   int rand2a = 0;
   int rand2b = 0;
   
   int rand3a = 0;
   int rand3b = 0;
   int rand4a = 0;
   int rand4b = 0;
   
   H1 = malloc(26*sizeof(int));
   H2 = malloc(26*sizeof(int));
   H3 = malloc((2*W2+1)*(2*W2+1)*sizeof(int));
   
   for (i=0; i<26;i++){
       H1[i] = 0;
       H2[i] = 0;
   }
   for (i=0; i<(2*W2+1)*(2*W2+1);i++){
       H3[i] = 0;
   }
   for (j=0; j<N;j++)
       for (i=0; i<N;i++){
          pos =  i+j*N;
          mem[pos] = 0;
       }
//   med = malloc(8*sizeof(double));
//   histo = malloc(256*sizeof(double));
//   sorted = malloc(256*sizeof(double));
//   maxi = malloc(10*sizeof(double));
      
   for(i=0; i<M; i++){
      for(j=0; j<N; j++){
         /* step 1. Determine the local window*/      
         if(i < W) {
            rand1a = i;
            rand1b = W;
         }
         else {
              if (i>M-W-1){
               rand1a = W;
               rand1b = M-i-1;
            }
            else{
               rand1a = W;
               rand1b = W;
            }
         }

         
         if(j < W) {
            rand2a = j;
            rand2b = W;
         }
         else {
            if (j > N-W-1){
               rand2a = W;
               rand2b = N-j-1;
            }
            else{
               rand2a = W;
               rand2b = W;
            }
         }
          
          
          
          
         for(k=0; k<26; k++){ 
             H1[k] = 0; H2[k] = 0;
         }
         for(k=-rand1a; k<=rand1b; k++){ 
            for(l=-rand2a; l<=rand2b; l++){ 
               tmp = A1[i+k][j+l];
               H1[(int)maximum((int)(ceil(tmp/10.0)-1),0)]++;
            }
         }
         Xmax = 0;
         for(k=0; k<26; k++) 
             if (H1[k]>Xmax){
                 Xmax = H1[k];
                 Ymax = k;
             }
         H2[Ymax] = 1;
         for(k=Ymax+1; k<26; k++){ 
             if (H1[k] != 0){
                 H2[k] = 1;
             }
             else break;
         }
         for(k=Ymax-1; k>= 0; k--){ 
             if (H1[k] != 0){
                 H2[k] = 1;
             }
             else break;
         }

         for(k=-rand1a; k<=rand1b; k++){ 
            for(l=-rand2a; l<=rand2b; l++){ 
               tmp = (int)ceil((A1[i+k][j+l])/10.0);
               mem[(int)((i+k)+(j+l)*M)] += (1.0 - H2[(int)maximum((int)tmp-1,0)]);
            }
         }
      }
   }
   
   
   // 2nd iteration!
//   p1 = (12.0/25.0)*(2.0*(double)W+1)*(2*(double)W+1);
//   p2 = (18.0/25.0)*(2.0*(double)W+1)*(2*(double)W+1);
   p1 = 0.0;
   p2 = (18.0/25.0)*(2.0*(double)W+1)*(2*(double)W+1);
         
   for(i=0; i<M; i++){
      for(j=0; j<N; j++){
          
         /* step 1. Determine the local window*/      
         if(i < W) {
            rand1a = i;
            rand1b = W;
         }
         else {
              if (i>M-W-1){
               rand1a = W;
               rand1b = M-i-1;
            }
            else{
               rand1a = W;
               rand1b = W;
            }
         }

         
         if(j < W) {
            rand2a = j;
            rand2b = W;
         }
         else {
            if (j > N-W-1){
               rand2a = W;
               rand2b = N-j-1;
            }
            else{
               rand2a = W;
               rand2b = W;
            }
         }
          
         /* step 1. Determine the local window*/      
         if(i < W2) {
            rand3a = i;
            rand3b = W2;
         }
         else {
              if (i>M-W2-1){
               rand3a = W2;
               rand3b = M-i-1;
            }
            else{
               rand3a = W2;
               rand3b = W2;
            }
         }

         
         if(j < W2) {
            rand4a = j;
            rand4b = W2;
         }
         else {
            if (j > N-W2-1){
               rand4a = W2;
               rand4b = N-j-1;
            }
            else{
               rand4a = W2;
               rand4b = W2;
            }
         }
          
        p1 = 0.0;
        p2 = (18.0/25.0)*(rand1a+rand1b+1)*(rand2a+rand2b+1);
          
         mem[i+j*M] = LARGE ((double)mem[i+j*M], p1, p2);
         
         ind = 0;
         for(k=-rand3a; k<=rand3b; k++){ 
            for(l=-rand4a; l<=rand4b; l++){ 
               tmp = A1[i+k][j+l];
               for (pos=0; pos<ind; pos++){ 
                   if  (tmp < H3[pos]){
                       tmp2 = H3[pos];
                       H3[pos] = tmp;
                       tmp = tmp2;
                   }
               }
               H3[ind] = tmp;
               ind++;
            }
         }
         
         pos = (int)((rand3a+rand3b+1)*(rand4a+rand4b+1)/2);
         m1 = maximum(H3[pos]-H3[pos-1],H3[pos+1]-H3[pos]);
         m2 = (double)absol(H3[pos]-A1[i][j]);
         hlp = LARGE (m2, 1.09*m1, 1.35*m1);
         if (mem[i+j*M]>0) mem[i+j*M] = minimum(mem[i+j*M],hlp);
      }
   }
          
   free(H1); 
   free(H2); 
   free(H3); 
}  /* End of callFuzzyShrink */


#define Im1      prhs[0]
#define WSIZE    prhs[1]
#define WSIZE2   prhs[2]

#define OUT plhs[0]

/**
*  The interaction with Matlab (mex):
*        nlhs = amount of output arguments (= 1)
*        nrhs = amount of input arguments (= 2)
*     *plhs[] = link to the output 
*     *prhs[] = link to the input 
*
**/
void mexFunction( int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[] ) {
    int row, col, i, M, N,W,W2;
    double **mem, **A1;
    
    if (nlhs!=1)
        mexErrMsgTxt("It requires one output arguments only [M1].");
    if (nrhs!=3)
       mexErrMsgTxt("this method requires three input argument [Im1,WSIZE]");

    /* Get input values */  
    M = mxGetM(Im1);
    N = mxGetN(Im1);
    W = mxGetScalar(WSIZE);
    W2 = mxGetScalar(WSIZE2);

    /**
    * Allocate memory for return matrices 
    **/
    OUT = mxCreateDoubleMatrix(M, N, mxREAL);  
    mem = mxGetPr(OUT);

    /**
    * Dynamic allocation of memory for the input array
    **/
    A1 = malloc(M*sizeof(int));
    for(i=0;i<M;i++)
      A1[i] = malloc(N*sizeof(double));

     /**
     * Convert ARRAY_IN and INPUT_MASK to 2x2 C arrays (MATLAB stores a two-dimensional matrix 
     * in memory as a one-dimensional array) 
     ***/
     for (col=0; col < N; col++)
         for (row=0; row < M; row++) {
             A1[row][col] = (mxGetPr(Im1))[row+col*M];
	      }
	      
    /* Call callFuzzyShrink function */ 
    CALL_MEM(A1,W,W2,mem,M,N);

    for(i=0;i<M;i++)  free(A1[i]);
    free(A1); 
}
/* end mexFcn*/