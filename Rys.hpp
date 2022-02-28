#ifndef UNOMOL_RYS_HPP
#define UNOMOL_RYS_HPP
#include <cstdlib>
#include <iostream>
#include <cmath>
#include "Util.hpp"
#include <cstring>

namespace unomol {

#define MAXROOTS 12

class Rys {
    double **allocate_d2(int n1,int n2) {
        double ** p = new double*[n1];
        for (size_t i=0; i<n1; ++i) {
            p[i] = new double[n2];
            memset(p[i],0x0,sizeof(double)*n2);
        }
        return p;
    }
    double ***allocate_d3(int n1,int n2,int n3) {
        double *** p = new double**[n1];
        for (size_t k=0; k<n1; ++k) p[k] = allocate_d2(n2,n3);
        return p;
    }

    void delete_d2(double **p,int n1) {
        for (int i=n1; i;) {
            --i;
            delete [] p[i];
        }
        delete [] p;
    }
    void delete_d3(double ***p,int n1,int n2) {
        for (int i=n1; i;) {
            --i;
            delete_d2(p[i],n2);
        }
        delete [] p;
    }
    double *** Gx;
    double *** Gy;
    double *** Gz;
    double ** binomial;
    double *roots;
    double *weights;
    double B00,B1,B1p,C,Cp;
    int maxr;
  public:
    Rys(int maxlv = 2) {
        int maxl = (maxlv > 2) ? maxlv : 2;
        maxr = 2*maxl + 1;
        std::cerr << " Rys maxroots = " << maxr << "\n";
        Gx = allocate_d3(maxr,maxr,maxr);
        Gy = allocate_d3(maxr,maxr,maxr);
        Gz = allocate_d3(maxr,maxr,maxr);
        binomial = allocate_d2(maxr+1,maxr+1);
        roots = new double[maxr];
        weights = new double[maxr];
        double * fact = new double[maxr];
        fact[0] = 1.;
        for (int i=1; i<maxr; ++i) {
            fact[i] = double(i) * fact[i-1];
        }
        for (int i=0; i<=maxr; ++i) {
            for (int j=0; j<=i; ++j) {
                binomial[i][j] = fact[i]/fact[i-j]/fact[j];
            }
        }
        delete [] fact;
    }

    ~Rys() {
        delete [] weights;
        delete [] roots;
        delete_d2(binomial,maxr+1);
        delete_d3(Gz,maxr,maxr);
        delete_d3(Gy,maxr,maxr);
        delete_d3(Gx,maxr,maxr);
    }

    constexpr double Shift(
        const double * ab,const double *cd,
        const int * lv1,const int *lv2,
        const int * lv3,const int *lv4,
        const int nroots
    ) const noexcept {
        int l2 = lv2[0];
        int m2 = lv2[1];
        int n2 = lv2[2];
        int l12 = lv1[0] + l2;
        int m12 = lv1[1] + m2;
        int n12 = lv1[2] + n2;
        int l4 = lv4[0];
        int m4 = lv4[1];
        int n4 = lv4[2];
        int l34 = lv3[0] + l4;
        int m34 = lv3[1] + m4;
        int n34 = lv3[2] + n4;
        double sum{0.};
        for (int ir=0; ir<nroots; ++ir) {
            double xs = ShiftKernel(Gx[ir],ab[0],cd[0],l12,l2,l34,l4);
            double ys = ShiftKernel(Gy[ir],ab[1],cd[1],m12,m2,m34,m4);
            double zs = ShiftKernel(Gz[ir],ab[2],cd[2],n12,n2,n34,n4);
            sum += xs * ys * zs * weights[ir];
        }
        return sum;
    }

    void Recur(
        const double *p,const double *q,const double *a,const double *c,
        const double& pxp,const double& qxp,const double& txp,
        int lvt12,int lvt34,int nroots
    ) noexcept {
        double pa[3],qc[3],pq[3];
        pa[0] = p[0]-a[0];
        pa[1] = p[1]-a[1];
        pa[2] = p[2]-a[2];

        qc[0] = q[0]-c[0];
        qc[1] = q[1]-c[1];
        qc[2] = q[2]-c[2];

        pq[0] = p[0]-q[0];
        pq[1] = p[1]-q[1];
        pq[2] = p[2]-q[2];
        double pq2 = pq[0]*pq[0] + pq[1]*pq[1] + pq[2]*pq[2];
        double w = pxp * qxp / txp;
        double t = w * pq2;
        calculate_roots(t,nroots);
        for (int ir=0; ir<nroots; ++ir) {
            double r = roots[ir];
            double dr = r/(1.+r);
            double fff = dr/txp;
            B00 = 0.5 * fff;
            B1 = (0.5 - B00 * qxp) / pxp;
            B1p = (0.5 - B00 * pxp) / qxp;
            C = pa[0] - qxp * pq[0] * fff;
            Cp = qc[0] + pxp * pq[0] * fff;
            RecurKernel (Gx[ir], lvt12, lvt34);
            C = pa[1] - qxp * pq[1] * fff;
            Cp = qc[1] + pxp * pq[1] * fff;
            RecurKernel (Gy[ir], lvt12, lvt34);
            C = pa[2] - qxp * pq[2] * fff;
            Cp = qc[2] + pxp * pq[2] * fff;
            RecurKernel (Gz[ir], lvt12, lvt34);
        }
    }

    void calculate_roots(double t,int nroots) noexcept {
        switch (nroots) {
        case 0:
            std::cerr << "Rys cannot have zero roots\n";
            exit(-1);
        case 1:
            return root1(t);
        case 2:
            return root2(t);
        case 3:
            return root3(t);
        case 4:
            return root4(t);
        case 5:
            return root5(t);
        default:
            break;
        }
        rootN(nroots,t);
    }

    void root1(double t) noexcept;
    void root2(double t) noexcept;
    void root3(double t) noexcept;
    void root4(double t) noexcept;
    void root5(double t) noexcept;
    void rootN(int nroots,double t) noexcept;

    constexpr double ShiftKernel(
        const double * const * G,
        const double& abx,const double& cdx,
        const int l12,const int l2,const int l34,const int l4
    ) const noexcept {
        const double *bnj = binomial[l2];
        const double *bnl = binomial[l4];
        double sum{0.0};
        if (fabs(cdx) < 1.e-14) {
            if (fabs(abx) < 1.e-14) {
                return G[l12][l34];
            }
            double x12t{1.};
            for (int i=0; i<=l2; ++i) {
                sum += bnj[i] * x12t * G[l12-i][l34];
                x12t *= abx;
            }
            return sum;
        }
        if (fabs(abx) < 1.-14) {
            double x34t{1.0};
            for (int j=0; j<=l4; ++j) {
                sum += bnl[j] * x34t * G[l12][l34-j];
                x34t *= cdx;
            }
            return sum;
        }
        double x12t{1.0};
        for (int i=0; i<=l2; ++i) {
            double x34t = bnj[i]*x12t;
            for (int j=0; j<=l4; ++j) {
                sum += bnl[j] * x34t * G[l12-i][l34-j];
                x34t *= cdx;
            }
            x12t *= abx;
        }
        return sum;
    }

    constexpr void RecurKernel(double **G,int l12,int l34) noexcept {
        G[0][0] = 1.0;
        G[0][1] = Cp;
        G[1][0] = C;
        G[1][1] = B00 + C * Cp;
        if (l12 < 2 && l34 < 2) return;
        for (int j = 1; j < l34; j++) {
            G[0][j + 1] = j * B1p * G[0][j - 1] + Cp * G[0][j];
            G[1][j + 1] = j * B1p * G[1][j - 1] + B00 * G[0][j] + Cp * G[1][j];
        }
        for (int i = 2; i <= l12; i++) {
            G[i][0] = (i - 1) * B1 * G[i - 2][0] + C * G[i - 1][0];
            G[i][1] = i * B00 * G[i - 1][0] + Cp * G[i][0];
            for (int j = 1; j < l34; j++) {
                G[i][j + 1] =j * B1p * G[i][j - 1] +
                             i * B00 * G[i - 1][j] + Cp * G[i][j];
            }
        }
    }
};

}
#endif