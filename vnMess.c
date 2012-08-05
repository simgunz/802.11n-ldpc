#include "mex.h"
#include "matrix.h"
#include "math.h"

void vnmess(int nmenok,int n,double *E, mxArray **mxA,double *r,double *M)
{
    int vn,cn,i,j;
    double *A;
    double sum;
    int sizeA;

    
                            
    for(vn=0;vn<n;vn++) {
        A = mxGetPr(mxA[vn]);
        sizeA = mxGetN(mxA[vn]);                
        
        for(j=0;j<sizeA;j++) {
            cn = (int)A[j];
            sum=0;
            for(i=0;i<sizeA;i++) {
                if((int)A[i]!=cn){
                    sum += E[((int)A[i]-1) + vn*nmenok];
                }
            }
            M[(cn-1) + vn*nmenok] = sum  - r[vn];
        }
    }
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double *E,*r,*M;
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

    plhs[0] = mxCreateDoubleMatrix(nmenok,n, mxREAL);
    M = mxGetPr(plhs[0]);

    vnmess(nmenok,n,E,A,r,M);
}