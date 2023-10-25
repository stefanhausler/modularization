/*==========================================================
 * Implements the MATLAB function
 * 
 * function [rl] = minidx(logitDiff,pn1,pn2,nPostNodes);
 *
 *     rl         ... double (nPostNodesx1)
 *     nl         ... double (nPostNodesx1)
 *     logitDiff  ... double (NLx1)
 *     pn1        ... double (NLx1)
 *     pn2        ... double (NLx1)
 *     nPostNodes ... double 
 *     
 *========================================================*/
 
#include "mex.h"

#define min(a,b) (((a)<(b))?(a):(b))

/* The computational routine */
void minidx(double *logitDiff,double *pn1,double *pn2,double *rl,double *nl,int n, int m)
{
    int i,j,k;


    for (i=0; i<m; i++) {
       rl[i] = 1e20;
    }

    /* get minimum values for indices pn1 and pn2 */
    for (i=0; i<n; i++) {

       j = (int)pn1[i]-1;
       k = (int)pn2[i]-1;

       /* verify values */
       if ((j<0)|(k<0)|(j>=m)|(k>=m)) {
          mexErrMsgIdAndTxt("MyToolbox:minidx","Incorrect index.");
       }

       rl[j] = min(rl[j],logitDiff[i]);
       rl[k] = min(rl[k],logitDiff[i]);
       nl[j] = nl[j] + 1;
       nl[k] = nl[k] + 1;
    }

}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double *logitDiff;
    double *pn1;
    double *pn2;
    double *nPostNodes;

    double *rl;                  /* output vector */
    double *nl;                  /* output vector */

    int n,m;


/* function [rl] = minidx(logitDiff,pn1,pn2,nPostNodes);
 *
 *     rl         ... double (nPostNodesx1)
 *     nl         ... double (nPostNodesx1)
 *     logitDiff  ... double (NLx1)
 *     pn1        ... double (NLx1)
 *     pn2        ... double (NLx1)
 *     nPostNodes ... double 
 *     
 *========================================================*/


    /* check for proper number of arguments */
    if(nrhs!=4) {
        mexErrMsgIdAndTxt("MyToolbox:minidx:nrhs","Four inputs required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:minidx:nlhs","Two outputs required.");
    }
    
    /* make sure the first input argument is type double */
    if( !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:minidx:notDouble","Input matrix logitDiff must be type double.");
    }
    /* make sure the second input argument is type double */
    if( !mxIsDouble(prhs[1]) || 
         mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("MyToolbox:minidx:notDouble","Input matrix pn1 must be type double.");
    }
    /* make sure the third input argument is type double */
    if( !mxIsDouble(prhs[2]) || 
         mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("MyToolbox:minidx:notDouble","Input matrix pn2 must be type double.");
    }
    /* make sure the fourth input argument is type double */
    if( !mxIsDouble(prhs[3]) || 
         mxIsComplex(prhs[3])) {
        mexErrMsgIdAndTxt("MyToolbox:minidx:notDouble","Input nPostNodes must be type double.");
    }
    
    
   
    /* check that number of rows in first input argument is 1 */
    if(mxGetN(prhs[0])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:minidx:notColumnVector","Input logitDiff must be a column vector.");
    }
    if(mxGetN(prhs[1])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:minidx:notColumnVector","Input pn1 must be a column vector.");
    }
    if(mxGetN(prhs[2])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:minidx:notColumnVector","Input pn2 must be a column vector.");
    }
    if ((mxGetM(prhs[3])!=1)||(mxGetN(prhs[3])!=1)) {
        mexErrMsgIdAndTxt("MyToolbox:minidx:notColumnVector","Input nPostNodes must be a scalar.");
    }

    if(mxGetM(prhs[1])!=mxGetM(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:minidx:notColumnVector","Input pn1 must be of size logitDiff.");
    }
    if(mxGetM(prhs[2])!=mxGetM(prhs[0])) {
        mexErrMsgIdAndTxt("MyToolbox:minidx:notColumnVector","Input pn2 must be of size logitDiff.");
    }


    /* create a pointers to the real data in the input matrix  */
    logitDiff = mxGetPr(prhs[0]);
    pn1 = mxGetPr(prhs[1]);
    pn2 = mxGetPr(prhs[2]);
    nPostNodes = mxGetPr(prhs[3]);

    /* get size information */
    n = (int)mxGetM(prhs[0]);  
    m = (int)nPostNodes[0];

    /* create the output scalar */
    plhs[0] = mxCreateDoubleMatrix((int)nPostNodes[0],1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix((int)nPostNodes[0],1,mxREAL);

    /* get a pointer to the real data in the output matrix */
    rl = mxGetPr(plhs[0]);
    nl = mxGetPr(plhs[1]);

    /* call the computational routine */
    minidx(logitDiff,pn1,pn2,rl,nl,n,m);
}
