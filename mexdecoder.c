#include "mex.h"
#include "matrix.h"
#include "math.h"

#define Y(x) ((x) > (0) ? (0) : (1))

double lntanh(double x)
{
    double result = -log(tanh(x/2));
    return isinf(result) ? 10000 : result;
}

void cnmess(int nmenok,int n,double *M, mxArray **mxB,double *E)
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

int decoder(int nmenok,int n,int nCW,double sigmaw2, mxArray **mxA, mxArray **mxB, double *rr, double *H, double *u_out)
{
    int i,j,cw,vn,cn,ii,iteration;
    double L[n];
    double *M,*E;
    double *A,*B,*r;
    double sum;
    int sizeA,sizeB,binsum,checkOK,yCap[n];
    bool checkNOK;

    iteration = 30;
    checkOK = 0;

    M = malloc(sizeof(double)*nmenok*n);
    E = malloc(sizeof(double)*nmenok*n);
    r = malloc(sizeof(double)*n);

    for (cw=0;cw<nCW;cw++)
    {
        for (j=0;j<n;j++) {
            for (i=0;i<nmenok;i++) {
                M[i+j*nmenok] = -rr[j+cw*n];
            }
            r[j] = rr[j+cw*n];
        }

        for (ii=0;ii<iteration;ii++) {

            cnmess(nmenok,n,M,mxB,E);

            for(vn=0;vn<n;vn++) {
                A = mxGetPr(mxA[vn]);
                sizeA = mxGetN(mxA[vn]);
                sum=0;
                for(i=0;i<sizeA;i++) {
                    sum += E[((int)A[i]-1) + vn*nmenok];
                }
                yCap[vn] = Y(sum - r[vn]);
            }

            checkNOK=false;
            for(cn=0;cn<nmenok;cn++) {
                binsum=0;
                for(vn=0;vn<n;vn++) {
                    binsum += H[cn+vn*nmenok]*yCap[vn];
                }
                if(binsum & 1) {
                    checkNOK = true;
                    break;
                }
            }

            if(!checkNOK) {
                checkOK++;
                break;
            } else {
                vnmess(nmenok,n,E,mxA,r,M);
            }
        }

        for(i=0;i<(n-nmenok);i++) {
            u_out[cw*(n-nmenok)+i] = yCap[i];
        }

    }

    free(M);
    free(E);
    free(r);
    return checkOK;
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double *r,*u_out,*H,*checkOK;
    double sigmaw2;
    int n,nmenok,nCW,i;
    const mxArray* temp_cell;

    nmenok = (int)mxGetScalar(prhs[0]);
    n = (int)mxGetScalar(prhs[1]);
    nCW = (int)mxGetScalar(prhs[2]);
    sigmaw2 = (double)mxGetScalar(prhs[3]);

    mxArray *A[n];
    mxArray *B[nmenok];

    temp_cell = prhs[4];
    for(i=0;i<n;i++) {
        A[i]=mxGetCell(temp_cell,i);
    }
    temp_cell = prhs[5];
    for(i=0;i<nmenok;i++) {
        B[i]=mxGetCell(temp_cell,i);
    }

    r = mxGetPr(prhs[6]);
    H = mxGetPr(prhs[7]);

    plhs[0] = mxCreateDoubleMatrix(1,nCW*(n-nmenok), mxREAL);
    u_out = mxGetPr(plhs[0]);

    plhs[1] = mxCreateDoubleScalar(0);
    checkOK = mxGetPr(plhs[1]);

    *checkOK = decoder(nmenok,n,nCW,sigmaw2,A,B,r,H,u_out);
}