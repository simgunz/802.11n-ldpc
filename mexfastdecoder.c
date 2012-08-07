#include "mex.h"
#include "matrix.h"
#include "math.h"

#define Y(x) ((x) > (0) ? (0) : (1))

double lntanh(double x)
{
    double result = -log(tanh(x/2));
    return isinf(result) ? 10000 : result;
}

void cnmess(int nmenok,int n,double *M, double *B,double *E)
{
    int vn,cn,i,j;
    double sum;
    int sig;

    for(cn=0;cn<nmenok;cn++) {
        for(i=1;i<=B[cn];i++) {
            vn = (int)B[cn+i*nmenok];
            sum = 0;
            sig = 1;
            for(j=1;j<=B[cn];j++) {
                if((int)B[cn+j*nmenok]!=vn){
                    sum += lntanh(fabs(M[cn+((int)B[cn+j*nmenok]-1)*(nmenok)]));
                    sig *= (M[cn+((int)B[cn+j*nmenok]-1)*(nmenok)] > 0) ? 1 : -1;
                }
            }
            E[cn+(vn-1)*(nmenok)] = sig*lntanh(sum);
        }
    }
}

void vnmess(int nmenok,int n,double *E, double *A,double *r,double *M)
{
    int vn,cn,i,j;
    double sum;

    for(vn=0;vn<n;vn++) {
        for(i=1;i<=A[vn];i++) {
            cn = (int)A[vn+i*n];
            sum = 0;
            for(j=1;j<=A[vn];j++) {
                if((int)A[vn+j*n]!=cn){
                    sum += E[((int)A[vn+j*n]-1) + vn*nmenok];
                }
            }
            M[(cn-1)+vn*nmenok] = sum - r[vn];
        }
    }
}

int decoder(int k,int n,int nCW,double sigmaw2, double *A, double *B, double *H, double *rr, int iter, double *u_out)
{
    int i,j,cw,vn,cn,ii,nmenok;
    double L[n];
    double *M,*E;
    double *r;
    double sum;
    int sizeA,sizeB,binsum,checkOK,yCap[n];
    bool checkNOK;

    nmenok = n-k;
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

        for (ii=0;ii<iter;ii++) {

            cnmess(nmenok,n,M,B,E);

            for(vn=0;vn<n;vn++) {
                sum=0;
                for(i=1;i<=A[vn];i++) {
                    sum += E[((int)A[vn+i*n]-1) + vn*nmenok];
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
                vnmess(nmenok,n,E,A,r,M);
            }
        }

        for(i=0;i<k;i++) {
            u_out[cw*k+i] = yCap[i];
        }
    }

    free(M);
    free(E);
    free(r);
    return checkOK;
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    double *A,*B,*H,*r,*u_out,*checkOK;
    double sigmaw2;
    int k,n,nCW,i,iter;
    const mxArray* temp_cell;

    k = (int)mxGetScalar(prhs[0]);
    n = (int)mxGetScalar(prhs[1]);
    nCW = (int)mxGetScalar(prhs[2]);
    sigmaw2 = (double)mxGetScalar(prhs[3]);
    A = mxGetPr(prhs[4]);
    B = mxGetPr(prhs[5]);
    H = mxGetPr(prhs[6]);
    r = mxGetPr(prhs[7]);
    iter = (int)mxGetScalar(prhs[8]);

    plhs[0] = mxCreateDoubleMatrix(1,nCW*k, mxREAL);
    u_out = mxGetPr(plhs[0]);
    plhs[1] = mxCreateDoubleScalar(0);
    checkOK = mxGetPr(plhs[1]);

    *checkOK = decoder(k,n,nCW,sigmaw2,A,B,H,r,iter,u_out);
}