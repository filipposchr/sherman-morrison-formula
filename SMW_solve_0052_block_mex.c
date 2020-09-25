
#if !defined(_WIN32)
 #define dgesv dgesv_ 
#endif
#if !defined(_WIN32)
 #define dgemm dgemm_ 
#endif

#include "mex.h"
#include "lapack.h"
#include "blas.h"
#include <string.h>
void select_block(double *iarray,double *block,size_t m,size_t cps, size_t cpe,size_t rps, size_t rpe)
{
    size_t i;
    for (i=rps;i<rpe;i++)
    {
        memcpy(&block[(cpe-cps)*(i-rps)],&iarray[i*m+cps],sizeof(double)*(cpe-cps));
    }
}
void vector_update(double *Xblock,double *Yblock,size_t m, size_t p,size_t n) /* Xblock will contain results */
{
    mwSignedIndex dims[2];
    long int info;
    long int *iPivot;
    mxArray *mxPivot;
    dims[0] = m;
    dims[1] = p;
    mxPivot = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
    iPivot = (long int*)mxGetData(mxPivot);
    size_t i; 
    for (i=0;i<m;i++) Yblock[i + i * m] = Yblock[i + i * m] + 1; 
    dgesv(&m,&n,Yblock,&m,iPivot,Xblock,&p,&info);
}
void LapackSMWCompute(double *M, double *P, double *Bv,size_t m, size_t p,size_t blk, double *Output)
{
    mxArray *X,*b, *Y, *C;  /* in/out arguments to DGESV*/
    double *ptrX,*ptrb,*ptrY,*ptrC;
    long int *iPivot;
    mxArray *mxPivot;
    mwSignedIndex dims[2];
    long int info;
    
    X = mxCreateDoubleMatrix(m, p, mxREAL);
    b = mxCreateDoubleMatrix(p, 1, mxREAL);
    Y = mxCreateDoubleMatrix(m, p, mxREAL);
    C = mxCreateDoubleMatrix(m, p, mxREAL);
    
    ptrX=mxGetPr(X);
    ptrb=mxGetPr(b);
    ptrY=mxGetPr(Y);
    ptrC=mxGetPr(C);
    
    memcpy(ptrX, M, m*p*sizeof(double));
    memcpy(ptrb, Bv, p*sizeof(double));
    memcpy(ptrY, M, m*p*sizeof(double));
    memcpy(ptrC, P, m*p*sizeof(double));
    
    size_t n=1;
    dims[0] = m;
    dims[1] = p;
    mxPivot = mxCreateNumericArray(2, dims, mxINT32_CLASS, mxREAL);
    iPivot = (long int*)mxGetData(mxPivot);
    dgesv(&m,&n,ptrX,&m,iPivot,ptrb,&p,&info);
    dgesv(&m,&p,ptrY,&m,iPivot,ptrC,&p,&info);
    
    mxDestroyArray(X);
    mxDestroyArray(Y);
    
    /* Preparation for main loops */
    
    ptrY = ptrC;
    
    size_t nblk=m/blk;
    size_t k;
    size_t j;
    mxArray *Xblock,*Yblock,*Yswath,*Yswath2,*Xblock_vector;
    double *ptrXblock,*ptrYblock,*ptrXblock_vector,*ptrYswath,*ptrYswath2;
    
    /* Memory allocation */
    
    Xblock = mxCreateDoubleMatrix(blk,blk, mxREAL);
    Xblock_vector = mxCreateDoubleMatrix(blk, 1, mxREAL);
    Yblock = mxCreateDoubleMatrix(blk, blk, mxREAL);
    Yswath = mxCreateDoubleMatrix(m, blk, mxREAL);
    Yswath2 = mxCreateDoubleMatrix(m, blk, mxREAL);
    ptrXblock=mxGetPr(Xblock);
    ptrXblock_vector=mxGetPr(Xblock_vector);
    ptrYblock=mxGetPr(Yblock);
    ptrYswath=mxGetPr(Yswath);
    ptrYswath2=mxGetPr(Yswath2);
    
    /* main loop */
    
    for (k=0;k<nblk-1;k++)
    {
        size_t sI;
        size_t eI;
        sI = k * blk;
        eI = (k + 1) * blk;
        select_block(ptrb,ptrXblock_vector,1,sI,eI,0,1);
        select_block(ptrY,ptrYblock,m,sI,eI,sI,eI);
        select_block(ptrY,ptrYswath,m,0,m,sI,eI);
        vector_update(ptrXblock_vector,ptrYblock,blk,blk,1);
        const char chn='N';
        double alpha=-1;
        double beta=1;
        long int one=1;
        dgemm(&chn, &chn, &m, &one,&blk, &alpha, ptrYswath, &m, ptrXblock_vector, &blk, &beta, ptrb, &m);
        for (j = (k+1);j<nblk;j++)
        {
            size_t sJ;
            size_t eJ;
            sJ = j * blk;
            eJ = (j+1) * blk;
            select_block(ptrY,ptrXblock,m,sI,eI,sJ,eJ);
            select_block(ptrY,ptrYblock,m,sI,eI,sI,eI);
            select_block(ptrY,ptrYswath,m,0,m,sI,eI);
            select_block(ptrY,ptrYswath2,m,0,m,sJ,eJ);
            vector_update(ptrXblock,ptrYblock,blk,blk,blk);
            const char chn='N';
            double alpha=-1;
            double beta=1;
            long int one=1;
            dgemm(&chn, &chn, &m, &blk,&blk, &alpha, ptrYswath, &m, ptrXblock, &blk, &beta, ptrYswath2, &m);
            memcpy(&ptrY[m*sJ],ptrYswath2,sizeof(double)*blk*m);
        }
        
        
        
    }
    size_t sI;
    size_t eI;
    sI = (nblk-1) * blk;
    eI = (nblk) * blk;
    select_block(ptrb,ptrXblock_vector,1,sI,eI,0,1);
    select_block(ptrY,ptrYblock,m,sI,eI,sI,eI);
    select_block(ptrY,ptrYswath,m,0,m,sI,eI);
    vector_update(ptrXblock_vector,ptrYblock,blk,blk,1);
    const char chn='N';
    const char chn2='N';
    double alpha=-1;
    double beta=1;
    long int one=1;
    dgemm(&chn, &chn2, &m, &one,&blk, &alpha, ptrYswath, &m, ptrXblock_vector, &blk, &beta, ptrb, &m);
    memcpy(Output, ptrb, m*sizeof(double));
    
    /* free memory */
    mxDestroyArray(Xblock);
    mxDestroyArray(Xblock_vector);
    mxDestroyArray(Yblock);
    mxDestroyArray(Yswath);
    mxDestroyArray(Yswath2);
    mxDestroyArray(C);
    
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *A;    /* pointers to input matrice */
    double *b;    /* pointers to input vector */
    double *X;    /* pointers to output vector */
    
    double *M,*P,*Q;    /* pointers to additional matrices */
    size_t p,m,n;     /* matrix dimensions */
    double blk;              /* size of block */
    
    /* Check for proper number of arguments. */
    if ( nrhs != 6) {
        mexErrMsgIdAndTxt("MATLAB:matrixDivide:rhs",
                "This function requires 6 inputs.");
    }
    
    A = mxGetPr(prhs[0]); /* pointers input parameters */
    b = mxGetPr(prhs[1]);
    M = mxGetPr(prhs[2]);
    P = mxGetPr(prhs[3]);
    Q = mxGetPr(prhs[4]);
    blk = mxGetScalar(prhs[5]);
    
    /* dimensions of input matrices */
    m = mxGetM(prhs[0]);
    p = mxGetN(prhs[0]);
    n = mxGetM(prhs[1]);
    
    /* Validate input arguments */
    if (p != mxGetN(prhs[1])) {
        mexErrMsgIdAndTxt("MATLAB:matchdims",
                "Inner dimensions of matrices do not match.");
    }
    if (p != m) {
        mexErrMsgIdAndTxt("MATLAB:square",
                "Input matrix 1 must be square.");
    }
    if (n != 1) {
        mexErrMsgIdAndTxt("MATLAB:zerodivide",
                "Input matrix 2 mist be a column vector.");
    }
    if (p != mxGetN(prhs[2])) mexErrMsgIdAndTxt("MATLAB:matchdims","Incorrect dimensions of parameter 3.");
    if (m != mxGetM(prhs[2])) mexErrMsgIdAndTxt("MATLAB:matchdims","Incorrect dimensions of parameter 3.");
    if (p != mxGetN(prhs[3])) mexErrMsgIdAndTxt("MATLAB:matchdims","Incorrect dimensions of parameter 4.");
    if (m != mxGetM(prhs[3])) mexErrMsgIdAndTxt("MATLAB:matchdims","Incorrect dimensions of parameter 4.");
    if (m % (int)blk != 0 ) mexErrMsgIdAndTxt("MATLAB:","Incorrect parameter 6.");
        
    plhs[0] = mxCreateDoubleMatrix(m, 1, mxREAL);
    double *out = mxGetPr(plhs[0]);        /*output vector*/
    LapackSMWCompute(M, P, b,m, p,blk,out);
    
    
}
