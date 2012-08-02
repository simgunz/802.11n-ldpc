#include "mex.h"
#include "matrix.h"
#include "math.h"

#define M(a) ((a) == (0) ? (a-1) : (a))

double lntanh(double x)
{
    double result = -log(tanh(x/2));
    return isinf(result) ? 10000 : result;
}

void decoder(int k,int n,double *M, mxArray **mxA, mxArray **mxB,double *E)
{

    int vn,cn,i,j;
    double *B;
    int sizeB;
    for(cn=0;cn<n-k;cn++) {
        B = mxGetPr(mxB[cn]);
        sizeB = mxGetN(mxB[cn]);
        for(j=0;j<sizeB;j++) {
            vn = (int)B[j];
            double sum=0;
            int sig=1;
            for(i=0;i<sizeB;i++) {
                if(B[i]!=vn){
                    sum += lntanh(fabs(M[cn + ((int)B[i]-1)*(n-k)]));
                    sig *= (M[cn + ((int)B[i]-1)*(n-k)] > 0) ? 1 : -1;
                }
            }
            E[cn+(vn-1)*(n-k)] = sig*lntanh(sum);
        }
    }
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double *M,*E;

    int n,k,i;
    const mxArray* temp_cell;


    k = mxGetScalar(prhs[0]);
    n = mxGetScalar(prhs[1]);


    M = mxGetPr(prhs[2]);


    k = n-k;

    mxArray *A[n],*B[n-k];

    temp_cell=prhs[3];
    for(i=0;i<n;i++)
    {
        A[i]=mxGetCell(temp_cell,i);
    }

    temp_cell=prhs[4];
    for(i=0;i<n-k;i++)
    {
        B[i]=mxGetCell(temp_cell,i);
    }

    plhs[0] = mxCreateDoubleMatrix(n-k,n, mxREAL);
    E = mxGetPr(plhs[0]);

    decoder(k,n,M,A,B,E);

}