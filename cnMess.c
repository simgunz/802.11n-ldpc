#include "mex.h"
#include "matrix.h"
#include "math.h"

double lntanh(double x)
{
    double result = -log(tanh(x/2));
    return isinf(result) ? 10000 : result;
}

void decoder(int nmenok,int n,double *M, mxArray **mxB,double *E)
{
    int vn,cn,i,j;
    double *B;
    double sum;
    int sizeB,sig;
    
    for(cn=0;cn<nmenok;cn++) {
        B = mxGetPr(mxB[cn]);
        sizeB = mxGetN(mxB[cn]);
        
        for(j=0;j<sizeB;j++) {
            vn = (int)B[j];
            sum=0;
            sig=1;
            for(i=0;i<sizeB;i++) {
                if((int)B[i]!=vn){
                    sum += lntanh(fabs(M[cn + ((int)B[i]-1)*(nmenok)]));
                    sig *= (M[cn + ((int)B[i]-1)*(nmenok)] > 0) ? 1 : -1;
                }
            }
            E[cn+(vn-1)*(nmenok)] = sig*lntanh(sum);
        }
    }
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double *M,*E;
    int n,nmenok,i;
    const mxArray* temp_cell;

    nmenok = mxGetM(prhs[0]);    
    n = mxGetN(prhs[0]);    
    
    M = mxGetPr(prhs[0]);    

    mxArray *B[nmenok];

    temp_cell = prhs[1];
    for(i=0;i<nmenok;i++)
    {
        B[i]=mxGetCell(temp_cell,i);
    }

    plhs[0] = mxCreateDoubleMatrix(nmenok,n, mxREAL);
    E = mxGetPr(plhs[0]);

    decoder(nmenok,n,M,B,E);     
}