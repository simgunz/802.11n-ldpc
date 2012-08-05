#include "mex.h"
#include "matrix.h"
#include "math.h"

void getL(int nmenok,int n,double *E, mxArray **mxA,double *r,double *L)
{
    int vn,i;
    double *A;
    double sum;
    int sizeA;
                                
    for(vn=0;vn<n;vn++) {
        A = mxGetPr(mxA[vn]);
        sizeA = mxGetN(mxA[vn]);                                            
        sum=0;
        for(i=0;i<sizeA;i++) {                
            sum += E[((int)A[i]-1) + vn*nmenok];                
        }
        L[vn] = sum  - r[vn];        
    }
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double *E,*r,*L;
    int n,nmenok,i;
    const mxArray* temp_cell;

    nmenok = mxGetM(prhs[0]);
    n = mxGetN(prhs[0]);

    E = mxGetPr(prhs[0]);
    r = mxGetPr(prhs[2]);

    mxArray *A[n];

    temp_cell = prhs[1];
    for(i=0;i<n;i++)
    {
        A[i]=mxGetCell(temp_cell,i);
    }

    plhs[0] = mxCreateDoubleMatrix(1,n, mxREAL);
    L = mxGetPr(plhs[0]);

    getL(nmenok,n,E,A,r,L);
}