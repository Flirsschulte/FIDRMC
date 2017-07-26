/**************************************************************************
% Fuzzy Fixed-valued Impulse Noise Reduction Method for Colour Images
%
%  The paper of the FIDRMC is proposed in: 
%
%  Stefan Schulte, Valérie De Witte, , Mike Nachtegael
%  Dietrich Van Der Weken and  Etienne E. Kerre:
%  A Fuzzy Two-step Filter for Impulse Noise Reduction From Colour Images,
%  IEEE Transactions on Image Processing 15(11), 2006, 3567 - 3578
%  
% Stefan Schulte (stefan.schulte@Ugent.be):
% Last modified: 15/01/06
%
%**************************************************************************/

#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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

double NOISE(double x, double ** noise, int amount) {
   double ret, tmp;
   int i;
   ret = 0.0; tmp = 0.0;
   
   for (i = 0; i<amount; i++){
      if((x == noise[i][0]) || ((x < noise[i][0]) && (x >= noise[i][1])) || ((x > noise[i][0]) && (x<=noise[i][3])))
          tmp = 1.0;
      else if((x<noise[i][1]) && (x>noise[i][2]))  tmp = ((noise[i][1]-x)/(noise[i][1]-noise[i][2]));
      else if ((x>noise[i][3])&&(x<noise[i][4])) tmp = ((noise[i][4]-x)/(noise[i][4]-noise[i][3]));
      else tmp = 0.0;
      ret = maximum(ret,tmp);
   }
   return ret;
}

double Median(double * x, int len) {
   double tmp;
   int i,j;
   for (i = 0; i<len; i++){   
       for (j = i+1; j<len; j++){
           if (x[j] < x[i]){
               tmp = x[i];
               x[i] = x[j];
               x[j] = tmp;
           }
       }
   }
   return x[(int)(len/2)];
}

double EQUALI(double x, double cen, double a, double b) {
   double res, diff;
   diff = absol(x-cen);
   if (diff <= a) res = 1.0;
   else if ((diff > a) && (diff <= b)) res = (b-diff)/(b-a);
   else res = 0;
   return res;
}


/**************************************************************************************
*  The main progam of the FIDRMC noise reduction method
***************************************************************************************/
void callDenoise(double **A1,double **A2,double **A3,double **fixed1,double **fixed2,double **fixed3,double *filt1,double *filt2,double *filt3,int M, int N,int W) { 
   int i,j,k,l, amount1, amount2, amount3,it, loop;   
   double som, wei;
   int *over_i, *over_j, OVER, OVER2, ok,teller;
   double noise1, noise2, noise3,cen1, cen2, cen3,glob;
   double noise1a, noise2a, noise3a, nR, nG, nB;
   double *RG, *RB, *GB, *NRG, *NRB, *NGB, *cor;
   double coR, coG, coB, soR, soG, soB, resR, resG, resB, resR2, resG2, resB2;
   double centRG, centRB, centGB, cenRG, cenRB, cenGB, memRG, memRB, memGB, mem ;
   double corr1, corr2;
   
   int rand1a = 0;
   int rand1b = 0;
   int rand2a = 0;
   int rand2b = 0;
        
   over_i = malloc(M*N*sizeof(double));
   over_j = malloc(M*N*sizeof(double));
   RG = malloc((2*W+1)*(2*W+1)*sizeof(double));
   RB = malloc((2*W+1)*(2*W+1)*sizeof(double));
   GB = malloc((2*W+1)*(2*W+1)*sizeof(double));
   NRG = malloc((2*W+1)*(2*W+1)*sizeof(double));
   NRB = malloc((2*W+1)*(2*W+1)*sizeof(double));
   NGB = malloc((2*W+1)*(2*W+1)*sizeof(double));
   cor = malloc((2*W+1)*(2*W+1)*sizeof(double));

   for (i=0; i<(2*W+1)*(2*W+1); i++) {
       RG[i] = 0;  GB[i] = 0; RB[i] = 0; cor[i] = 0;
       NRG[i] = 0;  NGB[i] = 0; NRB[i] = 0; 
   }
     
   for (i=0; i < M; i++)
         for (j=0; j < N; j++) {
             filt1[i+j*M] = A1[i][j];
             filt2[i+j*M] = A2[i][j];
             filt3[i+j*M] = A3[i][j];
	      }

   
   amount1 = fixed1[10][0];
   amount2 = fixed2[10][0];
   amount3 = fixed3[10][0];
   it = 1;
   OVER = 0;
   OVER2 = 0;
   while((it == 1)||(OVER != 0)) {
      // First iteration  
      if(it==1){
          for(i=2; i<M-2; i++){
              for(j=2; j<N-2; j++){
                   ok = 0;
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
                 /* end step 1. */      
         
                 cen1 = A1[i][j];
                 cen2 = A2[i][j];
                 cen3 = A3[i][j];
                 noise1a = NOISE(cen1, fixed1, amount1);
                 noise2a = NOISE(cen2, fixed2, amount2);
                 noise3a = NOISE(cen3, fixed3, amount3);
                
                 
                 glob = minimum(minimum(noise1a,noise2a),noise3a);

                 /************************************************
                 *     THE RED COMPONENT
                 ***********************************************/
                 noise1 = NOISE(cen1, fixed1, amount1);
                 if(noise1 > 0){
                    teller =0;
                    for (k=i-rand1a; k<=i+rand1b; k++){
                       for (l=j-rand2a; l<=j+rand2b; l++){
                          nR = NOISE(A1[k][l], fixed1, amount1);
                          nG = NOISE(A2[k][l], fixed2, amount2);
                          nB = NOISE(A3[k][l], fixed3, amount3);
                          
                          RG[teller] = A1[k][l] - A2[k][l];
                          RB[teller] = A1[k][l] - A3[k][l];
                          
                          NRG[teller] = 1.0 - maximum(nR,nG);
                          NRB[teller] = 1.0 - maximum(nR,nB);
                          cor[teller] = A1[k][l];
                          teller++;
                       }
                    }
                    
                    // one of the three components is noise-free
                    if (glob == 0) {
                        if ((noise2a==0) && (noise3a==0)) {
                            coG = 0; coB = 0; soG = 0; soB = 0;
                            for (k=0; k<teller; k++){
                                coG +=  NRG[k]*RG[k];
                                coB +=  NRB[k]*RB[k];
                                soG +=  NRG[k];
                                soB +=  NRB[k];
                            }
                            if (soB == 0) {
                                if (soG == 0) resR2 = cen1; 
                                else resR2 = minimum (maximum(cen2+(coG/soG),0),255); 
                            }
                            else if (soG == 0) resR2 = minimum (maximum(cen3+(coB/soB),0),255); 
                            else  resR2 = minimum (maximum((((cen3+(coB/soB))+(cen2+(coG/soG)))/2.),0),255); 
                        }
                        else if (noise2a==0) {
                            coG = 0; soG = 0;
                            for (k=0; k<teller; k++){
                                coG +=  NRG[k]*RG[k];
                                soG +=  NRG[k];
                            }
                            if (soG == 0) resR2 = cen1; 
                            else  resR2 = minimum (maximum(cen2+(coG/soG),0),255); 
                        }
                        else {
                            coB = 0; soB = 0;
                            for (k=0; k<teller; k++){
                                coB +=  NRB[k]*RB[k];
                                soB +=  NRB[k];
                            }
                            if (soB == 0) resR2 = cen1; 
                            else  resR2 = minimum (maximum(cen3+(coB/soB),0),255); 
                        }
                    }
                    else resR2 = -300;
                    
                    cenRG = Median(RG,teller); cenRG = RG[2];
                    cenRB = Median(RB,teller); cenRB = RB[2];
                    
                    centRG = absol(cen1-cen2);
                    centRB = absol(cen1-cen3);
                    
                    memRG = EQUALI(centRG,cenRG,10.0,20.0);
                    memRB = EQUALI(centRB,cenRB,10.0,20.0);
                    
                    corr1 = EQUALI(absol(cen1-cor[2]),absol(cor[5]-cor[3]),10.0,20.0);
                    corr2 = EQUALI(absol(cen1-cor[6]),absol(cor[5]-cor[3]),10.0,20.0);
                    mem = minimum(minimum(memRG,memRB), maximum(maximum(corr1,corr2),1.0-glob));
                    
                    if (resR2 == -300){
                       som = 0.0; wei = 0.0; 
                       for (k=i-rand1a; k<=i+rand1b; k++){
                          for (l=j-rand2a; l<=j+rand2b; l++){
                             wei += 1.0 - NOISE(A1[k][l], fixed1, amount1); 
                             som +=  A1[k][l]*(1.0 - NOISE(A1[k][l], fixed1, amount1)); 
                          }
                        }
                       
                        if (wei == 0) resR = Median(cor,teller);
                        else resR = (som/wei);
                    }
                    else resR = resR2;
                    
                    noise1 = NOISE(resR, fixed1, amount1);
                    if ((noise1 > 0) && (ok==0)){
                       over_i[OVER2] = i;
                       over_j[OVER2] = j;
                       OVER2++;

                       ok = 1;
                    }
                 }
                 else resR = cen1;

                 /************************************************
                 *     THE GREEN COMPONENT
                 ***********************************************/
                 noise2 = NOISE(cen2, fixed2, amount2);
                 if(noise2 > 0){
                    teller =0;
                    for (k=i-rand1a; k<=i+rand1b; k++){
                       for (l=j-rand2a; l<=j+rand2b; l++){
                          nR = NOISE(A1[k][l], fixed1, amount1);
                          nG = NOISE(A2[k][l], fixed2, amount2);
                          nB = NOISE(A3[k][l], fixed3, amount3);
                          
                          RG[teller] = A2[k][l] - A1[k][l];
                          GB[teller] = A2[k][l] - A3[k][l];
                          
                          NRG[teller] = 1.0 - maximum(nR,nG);
                          NGB[teller] = 1.0 - maximum(nG,nB);
                          cor[teller] = A2[k][l];
                          teller++;
                       }
                    }
                    
                    // one of the three components is noise-free
                    if (glob == 0) {
                        if ((noise1a==0) && (noise3a==0)) {
                            coR = 0; soR = 0; coB = 0; soB = 0;
                            for (k=0; k<teller; k++){
                                coR +=  NRG[k]*RG[k];
                                coB +=  NGB[k]*GB[k];
                                soR +=  NRG[k];
                                soB +=  NGB[k];
                            }
                            if (soB == 0) {
                                if (soR == 0) resG2 = cen2; 
                                else resG2 = minimum (maximum(cen1+(coR/soR),0),255); 
                            }
                            else if (soR == 0) resG2 = minimum (maximum(cen3+(coB/soB),0),255); 
                            else  resG2 = minimum ( maximum(((cen3+(coB/soB))+(cen1+(coR/soR)))/2.0 ,0) ,255); 
                        }
                        else if (noise1a==0) {
                            coR = 0; soR = 0;
                            for (k=0; k<teller; k++){
                                coR +=  NRG[k]*RG[k];
                                soR +=  NRG[k];
                            }
                            if (soR == 0) resG2 = cen2; 
                            else  resG2 = minimum (maximum(cen1+(coR/soR),0),255); 
                        }
                        else {
                            coB = 0; soB = 0;
                            for (k=0; k<teller; k++){
                                coB +=  NGB[k]*GB[k];
                                soB +=  NGB[k];
                            }
                            if (soB == 0) resG2 = cen2; 
                            else  resG2 = minimum (maximum(cen3+(coB/soB),0),255); 
                        }
                    }
                    else resG2 = -300;
                    
                    cenRG = Median(RG,teller); cenRG = RG[2];
                    cenGB = Median(GB,teller); cenGB = GB[2];
                    
                    centRG = absol(cen2-cen1);
                    centGB = absol(cen2-cen3);
                    
                    memRG = EQUALI(centRG,cenRG,10.0,20.0);
                    memGB = EQUALI(centGB,cenGB,10.0,20.0);
                    
                    corr1 = EQUALI(absol(cen2-cor[2]),absol(cor[5]-cor[3]),10.0,20.0);
                    corr2 = EQUALI(absol(cen2-cor[6]),absol(cor[5]-cor[3]),10.0,20.0);
                    mem = minimum(minimum(memRG,memGB), maximum(maximum(corr1,corr2),1.0-glob));
                    
                    if (resG2 == -300){
                       som = 0.0; wei = 0.0; 
                       for (k=i-rand1a; k<=i+rand1b; k++){
                          for (l=j-rand2a; l<=j+rand2b; l++){
                             wei += 1.0 - NOISE(A2[k][l], fixed2, amount2); 
                             som +=  A2[k][l]*(1.0 - NOISE(A2[k][l], fixed2, amount2)); 
                          }
                        }
                        if (wei == 0) resG = Median(cor,teller);
                        else resG = (som/wei);
                    }
                    else resG = resG2;
                    
                    noise2 = NOISE(resG, fixed2, amount2);
                    if ((noise2 > 0) && (ok==0)){
                       over_i[OVER2] = i;
                       over_j[OVER2] = j;
                       OVER2++;
                       ok = 1;
                    }
                 }
                 else resG = cen2;

                 /************************************************
                 *     THE BLUE COMPONENT
                 ***********************************************/
                 noise3 = NOISE(cen3, fixed3, amount3);
                 if(noise3 > 0){
                    teller =0;
                    for (k=i-rand1a; k<=i+rand1b; k++){
                       for (l=j-rand2a; l<=j+rand2b; l++){
                          nR = NOISE(A1[k][l], fixed1, amount1);
                          nG = NOISE(A2[k][l], fixed2, amount2);
                          nB = NOISE(A3[k][l], fixed3, amount3);
                          
                          RB[teller] = A3[k][l] - A1[k][l];
                          GB[teller] = A3[k][l] - A2[k][l];
                          
                          NRB[teller] = 1.0 - maximum(nR,nB);
                          NGB[teller] = 1.0 - maximum(nG,nB);
                          cor[teller] = A3[k][l];
                          teller++;
                       }
                    }
                    
                    // one of the three components is noise-free
                    if (glob == 0) {
                        if ((noise1a==0) && (noise2a==0)) {
                            coR = 0; soR = 0; coG = 0; soG = 0; 
                            for (k=0; k<teller; k++){
                                coG +=  NGB[k]*GB[k];
                                coR +=  NRB[k]*RB[k];
                                soG +=  NGB[k];
                                soR +=  NRB[k];
                            }
                            if (soR == 0) {
                                if (soG == 0) resB2 = cen3; 
                                else resB2 = minimum (maximum(cen2+(coG/soG),0),255); 
                            }
                            else if (soG == 0) resB2 = minimum (maximum(cen1+(coR/soR),0),255); 
                            else  resB2 = minimum ( maximum(((cen1+(coR/soR))+(cen2+(coG/soG)))/2.0 ,0) ,255); 
                        }
                        else if (noise2a==0) {
                            coG = 0; soG = 0; 
                            for (k=0; k<teller; k++){
                                coG +=  NGB[k]*GB[k];
                                soG +=  NGB[k];
                            }
                            if (soG == 0) resB2 = cen3; 
                            else  resB2 = minimum (maximum(cen2+(coG/soG),0),255); 
                        }
                        else {
                            coR = 0; soR = 0; 
                            for (k=0; k<teller; k++){
                                coR +=  NRB[k]*RB[k];
                                soR +=  NRB[k];
                            }
                            if (soR == 0) resB2 = cen3; 
                            else  resB2 = minimum (maximum(cen1+(coR/soR),0),255); 
                        }
                    }
                    else resB2 = -300;
                    
                    cenGB = Median(GB,teller); cenGB = GB[2];
                    cenRB = Median(RB,teller); cenRB = RB[2];
                    
                    centGB = absol(cen3-cen2);
                    centRB = absol(cen3-cen1);
                    
                    memGB = EQUALI(centGB,cenGB,10.0,20.0);
                    memRB = EQUALI(centRB,cenRB,10.0,20.0);
                    
                    corr1 = EQUALI(absol(cen3-cor[2]),absol(cor[5]-cor[3]),10.0,20.0);
                    corr2 = EQUALI(absol(cen3-cor[6]),absol(cor[5]-cor[3]),10.0,20.0);
                    mem = minimum(minimum(memGB,memRB), maximum(maximum(corr1,corr2),1.0-glob));
                    
                    if (resB2 == -300){
                       som = 0.0; wei = 0.0; 
                       for (k=i-rand1a; k<=i+rand1b; k++){
                          for (l=j-rand2a; l<=j+rand2b; l++){
                             wei += 1.0 - NOISE(A3[k][l], fixed3, amount3); 
                             som +=  A3[k][l]*(1.0 - NOISE(A3[k][l], fixed3, amount3)); 
                          }
                        }
                        if (wei == 0) resB = Median(cor,teller);
                        else resB = (som/wei);
                    }
                    else resB = resB2;
                    
                    noise3 = NOISE(resB, fixed3, amount3);
                    if ((noise3 > 0) && (ok==0)){
                       over_i[OVER2] = i;
                       over_j[OVER2] = j;
                       OVER2++;
                       ok = 1;
                    }
                 }
                 else resB = cen3;
                 filt1[i+j*M] = resR;
                 filt2[i+j*M] = resG;
                 filt3[i+j*M] = resB;
              }
          }
      }

       // Secound iteration  
      else {
         OVER2 = 0;
         for (loop = 0; loop<OVER;loop++){
             ok = 0;
                          
            i = over_i[loop];
            j = over_j[loop];

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
            /* end step 1. */      
            
            cen1 = A1[i][j];
            cen2 = A2[i][j];
            cen3 = A3[i][j];
            noise1a = NOISE(cen1, fixed1, amount1);
            noise2a = NOISE(cen2, fixed2, amount2);
            noise3a = NOISE(cen3, fixed3, amount3);
                 
            glob = minimum(minimum(noise1a,noise2a),noise3a);
                 

            /************************************************
            *     THE RED COMPONENT
            ***********************************************/
            noise1 = NOISE(cen1, fixed1, amount1);
            if(noise1 > 0){
               teller =0;
                for (k=i-rand1a; k<=i+rand1b; k+=W){
                    for (l=j-rand2a; l<=j+rand2b; l+=W){
                            nR = NOISE(filt1[k+l*M], fixed1, amount1);
                            nG = NOISE(filt2[k+l*M], fixed2, amount2);
                            nB = NOISE(filt3[k+l*M], fixed3, amount3);
                          
                            RG[teller] = filt1[k+l*M] - filt2[k+l*M];
                            RB[teller] = filt1[k+l*M] - filt3[k+l*M];
                          
                            NRG[teller] = 1.0 - maximum(nR,nG);
                            NRB[teller] = 1.0 - maximum(nR,nB);
                            cor[teller] = filt1[k+l*M];
                            teller++;
                    }
                }

                // one of the three components is noise-free
                if (glob == 0) {
                    if ((noise2a==0) && (noise3a==0)) {
                        coG = 0; soG = 0;  coB = 0; soB = 0;
                        for (k=0; k<teller; k++){
                            coG +=  NRG[k]*RG[k];
                            coB +=  NRB[k]*RB[k];
                            soG +=  NRG[k];
                            soB +=  NRB[k];
                        }
                        if (soB == 0) {
                            if (soG == 0) resR2 = cen1; 
                            else resR2 = minimum (maximum(cen2+(coG/soG),0),255); 
                        }
                        else if (soG == 0) resR2 = minimum (maximum(cen3+(coB/soB),0),255); 
                        else  resR2 = minimum ( maximum(((cen3+(coB/soB))+(cen2+(coG/soG)))/2.0 ,0) ,255); 
                    }
                    else if (noise2a==0) {
                        coG = 0; soG = 0;
                        for (k=0; k<teller; k++){
                           coG +=  NRG[k]*RG[k];
                           soG +=  NRG[k];
                        }
                        if (soG == 0) resR2 = cen1; 
                        else  resR2 = minimum (maximum(cen2+(coG/soG),0),255); 
                    }
                    else {
                        coB = 0; soB = 0;
                        for (k=0; k<teller; k++){
                            coB +=  NRB[k]*RB[k];
                            soB +=  NRB[k];
                        }
                        if (soB == 0) resR2 = cen1; 
                        else  resR2 = minimum (maximum(cen3+(coB/soB),0),255); 
                    }
                }
                else resR2 = -300;
            
                centRG = absol(cen1-cen2);
                centRB = absol(cen1-cen3);
                  
                cenRG = Median(RG,teller); cenRG = RG[2];
                cenRB = Median(RB,teller); cenRB = RB[2];
                memRG = EQUALI(centRG,cenRG,10.0,20.0);
                memRB = EQUALI(centRB,cenRB,10.0,20.0);
                mem = minimum(memRG,memRB);
                if (resR2 == -300){
                    som = 0.0; wei = 0.0; 
                    for (k=i-rand1a; k<=i+rand1b; k+=W){
                      for (l=j-rand2a; l<=j+rand2b; l+=W){
                             wei += 1.0 - NOISE(filt1[k+l*M], fixed1, amount1); 
                             som +=  filt1[k+l*M]*(1.0 - NOISE(filt1[k+l*M], fixed1, amount1)); 
                      }
                   }
                   if (wei == 0) resR = Median(cor,teller);
                   else resR = mem*cen1 + (1.0 - mem)*(som/wei);
                }
                else resR = resR2;
                noise1 = NOISE(resR, fixed1, amount1);
                if ((noise1 > 0) &&  (ok==0)){
                   over_i[OVER2] = i;
                   over_j[OVER2] = j;
                   OVER2++;
                   ok = 1;
                }
            }
            else resR = cen1;

                 
            /************************************************
            *     THE GREEN COMPONENT
            ***********************************************/
            noise2 = NOISE(cen2, fixed2, amount2);
            if(noise2 > 0){
               teller =0;
                for (k=i-rand1a; k<=i+rand1b; k+=W){
                    for (l=j-rand2a; l<=j+rand2b; l+=W){
                         nR = NOISE(filt1[k+l*M], fixed1, amount1);
                         nG = NOISE(filt2[k+l*M], fixed2, amount2);
                         nB = NOISE(filt3[k+l*M], fixed3, amount3);
                          
                         RG[teller] = filt2[k+l*M] - filt1[k+l*M];
                         GB[teller] = filt2[k+l*M] - filt3[k+l*M];
                          
                         NRG[teller] = 1.0 - maximum(nR,nG);
                         NGB[teller] = 1.0 - maximum(nG,nB);
                         cor[teller] = filt2[k+l*M];
                         teller++;
                    }
                }

                // one of the three components is noise-free
                if (glob == 0) {
                    if ((noise1a==0) && (noise3a==0)) {
                        coR = 0; soR = 0;  coB = 0; soB = 0;
                        for (k=0; k<teller; k++){
                            coR +=  NRG[k]*RG[k];
                            coB +=  NGB[k]*GB[k];
                            soR +=  NRG[k];
                            soB +=  NGB[k];
                        }
                        if (soB == 0) {
                            if (soR == 0) resG2 = cen2; 
                            else resG2 = minimum (maximum(cen1+(coR/soR),0),255); 
                        }
                        else if (soR == 0) resG2 = minimum (maximum(cen3+(coB/soB),0),255); 
                        else  resG2 = minimum ( maximum(((cen3+(coB/soB))+(cen1+(coR/soR)))/2.0 ,0) ,255); 
                    }
                    else if (noise1a==0) {
                        coR = 0; soR = 0;
                        for (k=0; k<teller; k++){
                           coR +=  NRG[k]*RG[k];
                           soR +=  NRG[k];
                        }
                        if (soR == 0) resG2 = cen2; 
                        else  resG2 = minimum (maximum(cen1+(coR/soR),0),255); 
                    }
                    else {
                        coB = 0; soB = 0;
                        for (k=0; k<teller; k++){
                            coB +=  NGB[k]*GB[k];
                            soB +=  NGB[k];
                        }
                        if (soB == 0) resG2 = cen2; 
                        else  resG2 = minimum (maximum(cen3+(coB/soB),0),255); 
                    }
                }
                else resG2 = -300;
            
                centRG = absol(cen2-cen1);
                centGB = absol(cen2-cen3);
                  
                cenRG = Median(RG,teller); cenRG = RG[2];
                cenGB = Median(GB,teller); cenGB = GB[2];
                memRG = EQUALI(centRG,cenRG,10.0,20.0);
                memGB = EQUALI(centGB,cenGB,10.0,20.0);
                mem = minimum(memRG,memGB);
                if (resG2 == -300){
                    som = 0.0; wei = 0.0; 
                    for (k=i-rand1a; k<=i+rand1b; k+=W){
                      for (l=j-rand2a; l<=j+rand2b; l+=W){
                          wei += 1.0 - NOISE(filt2[k+l*M], fixed2, amount2); 
                          som +=  filt2[k+l*M]*(1.0 - NOISE(filt2[k+l*M], fixed2, amount2)); 
                      }
                   }
                   if (wei == 0) resG = Median(cor,teller);
                   else resG = mem*cen2 + (1.0 - mem)*(som/wei);
                }
                else resG = resG2;
                noise2 = NOISE(resG, fixed2, amount2);
                if ((noise2 > 0) && (ok==0)){
                   over_i[OVER2] = i;
                   over_j[OVER2] = j;
                   OVER2++;
                   ok = 1;
                }
            }
            else resG = cen2;
                 
            /************************************************
            *     THE BLUE COMPONENT
            ***********************************************/
           noise3 = NOISE(cen3, fixed3, amount3);
            if(noise3 > 0){
               teller =0;
                for (k=i-rand1a; k<=i+rand1b; k+=W){
                    for (l=j-rand2a; l<=j+rand2b; l+=W){
                            nR = NOISE(filt1[k+l*M], fixed1, amount1);
                            nG = NOISE(filt2[k+l*M], fixed2, amount2);
                            nB = NOISE(filt3[k+l*M], fixed3, amount3);
                          
                            RB[teller] = filt3[k+l*M] - filt1[k+l*M];
                            GB[teller] = filt3[k+l*M] - filt2[k+l*M];
                          
                            NRB[teller] = 1.0 - maximum(nR,nB);
                            NGB[teller] = 1.0 - maximum(nG,nB);
                            cor[teller] = filt3[k+l*M];
                            teller++;
                    }
                }

               // one of the three components is noise-free
                if (glob == 0) {
                    if ((noise2a==0) && (noise1a==0)) {
                        coG = 0; soG = 0;  coR = 0; soR = 0;
                        for (k=0; k<teller; k++){
                            coG +=  NGB[k]*GB[k];
                            coR +=  NRB[k]*RB[k];
                            soG +=  NGB[k];
                            soR +=  NRB[k];
                        }
                        if (soR == 0) {
                            if (soG == 0) resB2 = cen3; 
                            else resB2 = minimum (maximum(cen2+(coG/soG),0),255); 
                        }
                        else if (soG == 0) resB2 = minimum (maximum(cen1+(coR/soR),0),255); 
                        else  resB2 = minimum ( maximum(((cen1+(coR/soR))+(cen2+(coG/soG)))/2.0 ,0) ,255); 
                    }
                    else if (noise2a==0) {
                        coG = 0; soG = 0;
                        for (k=0; k<teller; k++){
                           coG +=  NGB[k]*GB[k];
                           soG +=  NGB[k];
                        }
                        if (soG == 0) resB2 = cen3; 
                        else  resB2 = minimum (maximum(cen2+(coG/soG),0),255); 
                    }
                    else {
                        coR = 0; soR = 0;
                        for (k=0; k<teller; k++){
                            coR +=  NRB[k]*RB[k];
                            soR +=  NRB[k];
                        }
                        if (soR == 0) resB2 = cen3; 
                        else  resB2 = minimum (maximum(cen1+(coR/soR),0),255); 
                    }
                }
                else resB2 = -300;
            
                centGB = absol(cen3-cen2);
                centRB = absol(cen3-cen1);
                  
                cenGB = Median(GB,teller); cenGB = GB[2];
                cenRB = Median(RB,teller); cenRB = RB[2];
                memGB = EQUALI(centGB,cenGB,10.0,20.0);
                memRB = EQUALI(centRB,cenRB,10.0,20.0);
                mem = minimum(memGB,memRB);
                    
                if (resB2 == -300){
                    som = 0.0; wei = 0.0; 
                    for (k=i-rand1a; k<=i+rand1b; k+=W){
                      for (l=j-rand2a; l<=j+rand2b; l+=W){
                          if ((k>0) && (l>0) && (k<M) && (l<N)){
                             wei += 1.0 - NOISE(filt3[k+l*M], fixed3, amount3); 
                             som +=  filt3[k+l*M]*(1.0 - NOISE(filt3[k+l*M], fixed3, amount3)); 
                          }
                      }
                   }
                   if (wei == 0) resB = Median(cor,teller);
                   else resB = mem*cen3 + (1.0 - mem)*(som/wei);
                }
                else resB = resB2;
                noise3 = NOISE(resB, fixed3, amount3);
                if ((noise3 > 0) && (ok==0)){
                   over_i[OVER2] = i;
                   over_j[OVER2] = j;
                   OVER2++;
                   ok = 1;
                }
            }
            else resB = cen3;

            filt1[i+j*M] = resR;
            filt2[i+j*M] = resG;
            filt3[i+j*M] = resB;
           
         }// end for
      }//end else
      
      it++;
      W++;
      if ((OVER == OVER2)||(it==15)) break;
      OVER = OVER2;
   }
   free(RG); free(RB); free(GB); 
   free(NRG); free(NRB); free(NGB); 
   free(cor); 
   free(over_i); 
   free(over_j); 
}  


#define Im1      prhs[0]
#define Im2      prhs[1]
#define Im3      prhs[2]
#define FIXED1   prhs[3]
#define FIXED2   prhs[4]
#define FIXED3   prhs[5]
#define WSIZE    prhs[6]

#define OUT1 plhs[0]
#define OUT2 plhs[1]
#define OUT3 plhs[2]

/**
*  The interaction with Matlab (mex):
*        nlhs = amount of output arguments (= 3)
*        nrhs = amount of input arguments (= 7)
*     *plhs[] = link to the output 
*     *prhs[] = link to the input 
*
**/
void mexFunction( int nlhs, mxArray  *plhs[], int nrhs, const mxArray  *prhs[] ) {
    int row, col, i, M, N, M2, N2,W;
    double **fixed1,**fixed2,**fixed3, **A1, **A2, **A3, **filt1,**filt2,**filt3;
    
    if (nlhs!=3)
        mexErrMsgTxt("It requires three output arguments only [OUT1 OUT2 OUT3].");
    if (nrhs!=7)
       mexErrMsgTxt("this method requires 7 input argument [R,G,B, FIXED1,FIXED2,FIXED3, WSIE]");

    /* Get input values */  
    M = mxGetM(Im1);
    N = mxGetN(Im1);

    M2 = mxGetM(FIXED1);
    N2 = mxGetN(FIXED1);
    
    W = mxGetScalar(WSIZE);

    /**
    * Allocate memory for return matrices 
    **/
    OUT1 = mxCreateDoubleMatrix(M, N, mxREAL);  
    filt1 = mxGetPr(OUT1);

    OUT2 = mxCreateDoubleMatrix(M, N, mxREAL);  
    filt2 = mxGetPr(OUT2);
    
    OUT3 = mxCreateDoubleMatrix(M, N, mxREAL);  
    filt3 = mxGetPr(OUT3);

    /**
    * Dynamic allocation of memory for the input array
    **/
    A1 = malloc(M*sizeof(int));
    for(i=0;i<M;i++)
      A1[i] = malloc(N*sizeof(double));

    A2 = malloc(M*sizeof(int));
    for(i=0;i<M;i++)
      A2[i] = malloc(N*sizeof(double));

    A3 = malloc(M*sizeof(int));
    for(i=0;i<M;i++)
      A3[i] = malloc(N*sizeof(double));

    fixed1 = malloc(M2*sizeof(int));
    for(i=0;i<M2;i++)
      fixed1[i] = malloc(N2*sizeof(double));

    fixed2 = malloc(M2*sizeof(int));
    for(i=0;i<M2;i++)
      fixed2[i] = malloc(N2*sizeof(double));

    fixed3 = malloc(M2*sizeof(int));
    for(i=0;i<M2;i++)
      fixed3[i] = malloc(N2*sizeof(double));

     /**
     * Convert ARRAY_IN and INPUT_MASK to 2x2 C arrays (MATLAB stores a two-dimensional matrix 
     * in memory as a one-dimensional array) 
     ***/
     for (col=0; col < N; col++)
         for (row=0; row < M; row++) {
             A1[row][col] = (mxGetPr(Im1))[row+col*M];
	      }
	      
     for (col=0; col < N; col++)
         for (row=0; row < M; row++) {
             A2[row][col] = (mxGetPr(Im2))[row+col*M];
	      }
	      
     for (col=0; col < N; col++)
         for (row=0; row < M; row++) {
             A3[row][col] = (mxGetPr(Im3))[row+col*M];
	      }
	      
     for (col=0; col < N2; col++)
         for (row=0; row < M2; row++) {
             fixed1[row][col] = (mxGetPr(FIXED1))[row+col*M2];
	      }
	      
     for (col=0; col < N2; col++)
         for (row=0; row < M2; row++) {
             fixed2[row][col] = (mxGetPr(FIXED2))[row+col*M2];
	      }
	      
     for (col=0; col < N2; col++)
         for (row=0; row < M2; row++) {
             fixed3[row][col] = (mxGetPr(FIXED3))[row+col*M2];
	      }
	     	      
	      
    /* Call callFuzzyShrink function */ 
    callDenoise(A1,A2,A3,fixed1,fixed2,fixed3,filt1,filt2,filt3,M,N,W);

    for(i=0;i<M;i++)  free(A1[i]);
    free(A1); 
    for(i=0;i<M;i++)  free(A2[i]);
    free(A2); 
    for(i=0;i<M;i++)  free(A3[i]);
    free(A3); 
    for(i=0;i<M2;i++)  free(fixed1[i]);
    free(fixed1); 
    for(i=0;i<M2;i++)  free(fixed2[i]);
    free(fixed2); 
    for(i=0;i<M2;i++)  free(fixed3[i]);
    free(fixed3); 
}
/* end mexFcn*/