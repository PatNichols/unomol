
#include "SymmPack.hpp"

namespace unomol {
namespace SymmPack {

double TraceSymmPackProduct(
    double a[],double b[],int n) noexcept {
    double sum=0.0;
    int ij=0;
    for (int i=0; i<n; ++ij,++i) {
        for (int j=0; j<i; ++ij,++j) {
            sum+=a[ij]*b[ij];
        }
        sum+=0.5*a[ij]*b[ij];
    }
    return (sum+sum);
}

double SymmPackDiffNorm(
    double a[],double b[],int n) noexcept {
    double term=0.0;
    double sum=0.0;
    int ij=0;
    for (int i=0; i<n; ++i) {
        for (int j=0; j<i; ++ij,++j) {
            term=a[ij]-b[ij];
            sum+=term*term;
        }
        term=a[ij]-b[ij];
        sum+=0.5*term*term;
        ++ij;
    }
    sum=sum+sum;
    return (sqrt(sum)/n);
}



/////////////////////////////////////////////////////////////////////
//  This class finds the eigenvectors and eigenvalues of a
//   Symmetric packed matrix.
//   rsp Input:
//      Apack is SymmetricPacked matrix, it is destroyed on output
//      Evecs is the Eigenvector matrix.
//      Evals is the Eigenvalue vector.
//      n is the row/column dimension of Apack.
//   The only other method for public use is pythag
//      which finds the square root of (x*x+y*y).
/////////////////////////////////////////////////////////////////////
/**
 *  return the value sqrt(x*x+y*y) withput underflow
 *   or overflow error.
 */
inline double pythag(double x,double y) noexcept {
    double ax=fabs(x);
    double ay=fabs(y);
    if (ax>ay) {
        double t=ay/ax;
        return ax*sqrt(1.0+t*t);
    }
    if (ay==0.0) return 0.0;
    double t=ax/ay;
    return ay*sqrt(1.0+t*t);
}

void reduce(const int n,double* __restrict__ a,double* __restrict__ d,
            double* __restrict__ e) noexcept {
    int i,j,k;
    double scale,h,f,g,hh;
    int l,iz,jk;

    for (i = n; i--;) {
        l = i - 1;
        iz = i * (i + 1) / 2 - 1;
        h = scale = 0.0;
        for (k = 0; k < i ; k++) {
            d[k] = a[++iz];
            scale += fabs(d[k]);
        }
        if (scale == 0.0 || !i) {
            e[i] = 0.0;
            d[i] = a[++iz];
            a[iz] = 0.0;
            continue;
        }
        f = 1.0 / scale;
        for (k = 0; k <= l; k++) {
            d[k] *= f;
            h += d[k] * d[k];
        }
        f = d[l];
        g=-sqrt(h);
        if (f<0) g=-g;
        e[i] = scale * g;
        h -= f * g;
        d[l] = f - g;
        a[iz] = scale * d[l];
        if (i==1) {
            d[i] = a[++iz];
            a[iz] = scale * sqrt(h);
            continue;
        }
        jk=0;
        for (j=0; j<=l; j++) {
            f=d[j];
            g=0.0;
            for (k=0; k<j; k++) {
                g+=a[jk]*d[k];
                e[k]+=a[jk]*f;
                jk++;
            }
            e[j]=g+a[jk]*f;
            jk++;
        }
        /* form p */
        f=0.0;
        for (j=0; j<i; j++) {
            e[j]/=h;
            f+=e[j]*d[j];
        }
        hh=f/2.0/h;
        /* form q */
        for (j=0; j<i; j++) {
            e[j]-=hh*d[j];
        }
        jk=0;
        /* form reduced a */
        for (j=0; j<=l; j++) {
            f=d[j];
            g=e[j];
            for (k=0; k<=j; k++) {
                a[jk]-=(f*e[k]+g*d[k]);
                jk++;
            }
        }
        d[i]=a[++iz];
        a[iz]=scale*sqrt(h);
    }
    return;
}

int tql2(const int n,double* __restrict__ d,double* __restrict__ e,
         double* __restrict__ z,const int ldz) noexcept {
    int i,j,k,l,m,ii,l1, nm1;
    double b,c,c2,c3,dl1,ds1,el1,f,g,h,p,r,s,s2;
    double *zp0,*zp1;

    nm1=n-1;
    for (j=0; j<nm1; ++j) e[j]=e[j+1];
    e[nm1]=0.0;
    f=b=0.0;
    for (l=0; l<n; ++l) {
        j=0;
        h=fabs(d[l])+fabs(e[l]);
        if (b<h) b=h;
        m=l;
        while (m<(n-1) && (b+fabs(e[m]))!=b) ++m;
        if (m==l) {
            d[l]+=f;
            continue;
        }
        do {
            if (j==30) {
                fprintf(stderr,"Eigenvalue % d not found after 30 iterations\n",l+1);
                return (l+1);
            }
            ++j;
            l1=l+1;
            g=d[l];
            p=(d[l1]-g)/(2.0*e[l]);
            r=pythag(p,1.0);
            ds1=r;
            if (p<0) ds1=-r;
            ds1+=p;
            d[l]=e[l]/ds1;
            d[l1]=e[l]*ds1;
            dl1=d[l1];
            h=g-d[l];
            for (i=l1+1; i<n; ++i) d[i]-=h;
            f+=h;
            // ql transformation //
            p=d[m];
            c3=c2=c=1.0;
            el1=e[l1];
            s2=s=0.0;
            for (i=m-1; i>=l; i--) {
                c3=c2;
                c2=c;
                s2=s;
                g=c*e[i];
                h=c*p;
                r=pythag(p,e[i]);
                e[i+1]=s*r;
                s=e[i]/r;
                c=p/r;
                p=c*d[i]-s*g;
                d[i+1]=h+s*(c*g+s*d[i]);
                zp0=z+i*ldz;
                zp1=z+(i+1)*ldz;
                for (k=0; k<n; ++zp0,++zp1,++k) {
                    h=(*zp1);
                    (*zp1)=c*h+s*(*zp0);
                    (*zp0)=c*(*zp0)-s*h;
                }
            }
            p= (-s*s2*c3*el1*e[l])/dl1;
            e[l]=s*p;
            d[l]=c*p;
        } while ((b+fabs(e[l]))>b);
        d[l]+=f;
    }
    for (ii=1; ii<n; ++ii) {
        i=ii-1;
        k=i;
        p=d[i];
        for (j=ii; j<n; ++j) {
            if (d[j]<p) {
                k=j;
                p=d[j];
            }
        }
        if (k==i) continue;
        d[k]=d[i];
        d[i]=p;
        zp0=z+i*ldz;
        zp1=z+k*ldz;
        for (j=0; j<n; ++zp0,++zp1,++j) {
            p=(*zp0);
            (*zp0)=(*zp1);
            (*zp1)=p;
        }
    }
    return 0;
}

void
trback (const int n, double *__restrict__ a,
        double *__restrict__ z, const int ldz) noexcept {
    int ik, k, j, i;
    int l, iz;
    double s, h;
    double *zp;

    for (i = 1; i < n; i++) {
        l = i - 1;
        iz = (i * (l) / 2) + i;
        ik = iz + i;
        h = a[ik];
        if (h == 0.0)
            continue;
        for (j = 0; j < n; j++) {
            s = 0.0;
            ik = iz;
            zp = z + j*ldz;
            for (k = 0; k <= l; ++ik, k++) {
                s += a[ik] * (*zp);
                ++zp;
            }
            s = (s / h) / h;
            ik = iz;
            zp = z + j*ldz;
            for (k = 0; k <= l; ++ik, k++) {
                (*zp) -= a[ik] * s;
                ++zp;
            }
        }
    }
}


void rsp(const int n,
         double* __restrict__ a,
         double* __restrict__ z,
         double* __restrict__ evals,
         double* __restrict__ tmp) noexcept {
    int ierr;

    reduce(n,a,evals,tmp);
    copy_identity(n,n,z,n);
    ierr=tql2(n,evals,tmp,z,n);
    if (ierr) {
        fprintf(stderr,"Eigenvalue %d was not found \n",ierr+1);
        fatal_error("Singular matrix!\n");
    }
    trback(n,a,z,n);
    copy_trans(n,z);
}

/*
 assume Z*B*Z'=B(diagonal) this finds
  Z*A*Z'
*/
void sp_trans (int n, double *a, double *z, double *tmp) noexcept {
    int i, j, k;
    int ij, ik, ii, istr;
    double s, *zij, *zji, *zjk, *zik, *zi, *tjk;

    for (i = 0; i < n; i++) {
        zij = z + i * n;
        zji = z + i;
        for (j = 0; j < i; j++) {
            s = (*(zij));
            (*zij++) = (*zji);
            (*zji) = s;
            zji += n;
        }
    }
    for (i = 0; i < n; i++) {
        ii = i * (i + 1) / 2;
        for (j = 0; j < n; j++) {
            zjk = z + j * n;
            ik = ii;
            s = 0.0;
            for (k = 0; k <= i; ik++, k++)
                s += a[ik] * (*zjk++);
            --ik;
            for (k = (i + 1); k < n; istr++, k++) {
                ik=ik+k;
                s += a[ik] * (*zjk++);
            }
            (*(tmp + j * n + i)) = s;
        }
    }
    for (i = 0; i < n; i++) {
        zi = z + i * n;
        ij = i * (i + 1) / 2;
        for (j = 0; j <= i; ij++, j++) {
            zik = zi;
            tjk = tmp + j * n;
            s = 0.0;
            for (k = 0; k < n; k++)
                s += (*zik++) * (*tjk++);
            a[ij] = s;
        }
    }
    for (i = 0; i < n; i++) {
        zij = z + i * n;
        zji = z + i;
        for (j = 0; j < i; j++) {
            s = (*(zij));
            (*zij++) = (*zji);
            (*zji) = s;
            zji += n;
        }
    }
    return;
}

}
}
