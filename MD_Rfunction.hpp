#ifndef UNOMOL_MD_RFUNCTION_HPP
#define UNOMOL_MD_RFUNCTION_HPP
#include <iostream>
#include <cmath>
#include "Util.hpp"
using namespace std;

namespace unomol {

/////////////////////////////////////////////////////////////
//  C++ Implementation of the MacMurchieDavidson Rfunction
//  Ref: JCP 26,(218-231),1978
/////////////////////////////////////////////////////////////
class MD_Rfunction {
private:
    int lmax,asize;
    double*** r;
    double**** rz;
    double* dfact;
public:
    MD_Rfunction(int maxl) {
        lmax=maxl;
        asize=4*maxl+1;
        r=new_tensor3<double>(asize,asize,asize);
        rz=new_tensor4<double>(asize,asize,asize,asize);
        dfact = new double[(asize+1)];
        dfact[0] = 1.0;
        for (int i=1; i<=asize; ++i) dfact[i] = dfact[i-1] * ( i + i - 1);
    }

    ~MD_Rfunction() {
        delete [] dfact;
        delete_tensor4<double>(rz,asize,asize,asize);
        delete_tensor3<double>(r,asize,asize);
        r=0;
    }

    /** returns value of tensor for indices i,j, and k */
    inline double getValue(int i,int j,int k) const noexcept {
        return r[i][j][k];
    }

    /** return a page of the tensor for fast access */
    inline const double* getRow(int i,int j) const noexcept {
        return r[i][j];
    }

    /** find values for tensor */
    constexpr void eval(double sr,double t,double w,
                        double pq[],int ltot) noexcept {
        double * r0 = rz[0][0][0];
        Fgamma (r0, t, ltot);
        double wterm = -(w+w);
        double sterm = sr;
        for (int m = 0; m <= ltot; ++m) {
            r0[m] = sterm * r0[m];
            sterm *= wterm;
        }
        if ( ltot <= 12) {
            unroll_eval(pq,ltot);
        } else {
            loop_eval(pq,ltot);
        }
    }

    /**
     *  this function an array with the value for the
     *   integral from 0 to 1 of  pow(u,2*m)exp(-t*u*u) du.
     */
    constexpr void Fgamma(double fm[],double t,int m) const noexcept {
        const double tcrit=20.0;
        const double sqrtpi=0.88622692545275801365;
        const double eps=1.e-15;
        if (t>tcrit) {
            double oot = 1.0/t;
            fm[0]=sqrtpi*sqrt(oot);
            for (int i=1; i<=m; i++) fm[i]=fm[i-1]*(i-0.5)*oot;
            return;
        }
        if (t<eps) {
            for (int i=0; i<=m; ++i) fm[i] = 0.5/(0.5+i);
            return;
        }
        double mphalf=m+0.5;
        double term=0.5/mphalf;
        double sum=term;
        for (int i=1; i<=1000; i++) {
            term*=(t/(mphalf+i));
            sum+=term;
            if (term<eps) break;
        }
        double twot=2.0*t;
        double expt=exp(-t);
        fm[m]=sum*expt;
        for (int i=m-1; i>=0; --i) fm[i]=(fm[i+1]*twot+expt)/(i+i+1.0);
    }


    constexpr void loop_eval(const double *pq,int ltot) noexcept {
        // slower evaluation using loops and rz
        const double x = pq[0];
        const double y = pq[1];
        const double z = pq[2];
        // lz = 1
        int m_max = ltot-1;
        for (int m=0; m<=m_max; ++m) {
            rz[0][0][1][m] = z * rz[0][0][0][m+1];
        }
        r[0][0][1] = rz[0][0][1][0];
        // lz = 2
        m_max = ltot-2;
        for (int m=0; m<=m_max; ++m) {
            rz[0][0][2][m] = z * rz[0][0][1][m+1] + rz[0][0][0][m+1];
        }
        r[0][0][2] = rz[0][0][2][0];
        // lz > 2
        for (int lz=3; lz<=ltot; ++lz) {
            const int lzm1 = lz - 1;
            const int lzm2 = lz - 2;
            m_max = ltot - lz;
            for (int m=0; m<=m_max; ++m) {
                rz[0][0][lz][m] = z * rz[0][0][lzm1][m+1] +
                                  lzm1 * rz[0][0][lzm2][m+1];
            }
            r[0][0][lz] = rz[0][0][lz][0];
        }
        // ly = 1
        int lzend = ltot - 1;
        for (int lz=0; lz<=lzend; ++lz) {
            m_max = lzend - lz;
            for (int m=0; m<=m_max; ++m) {
                rz[0][1][lz][m] = y * rz[0][0][lz][m+1];
            }
            r[0][1][lz] = rz[0][1][lz][0];
        }
        // ly = 2
        lzend = ltot - 2;
        for (int lz=0; lz<=lzend; ++lz) {
            m_max = ltot - lz;
            for (int m=0; m<=m_max; ++m) {
                rz[0][2][lz][m] = y * rz[0][1][lz][m+1] + rz[0][0][lz][m+1];
            }
            r[0][2][lz] = rz[0][2][lz][0];
        }
        // ly > 2
        for (int ly=3; ly<=ltot; ++ly) {
            lzend = ltot - ly;
            for (int lz=0; lz<=lzend; ++lz) {
                int lym1 = ly - 1;
                int lym2 = ly - 2;
                m_max = lzend - lz;
                for (int m=0; m<=m_max; ++m) {
                    rz[0][ly][lz][m+1] =
                        y * rz[0][lym1][lz][m+1] +
                        lym1 * rz[0][lym2][lz][m+1];
                }
                r[0][ly][lz] = rz[0][ly][lz][0];
            }
        }
        // lx = 1;
        int lyend = ltot - 1;
        for (int ly=0; ly<=lyend; ++ly) {
            lzend = lyend - ly;
            for (int lz=0; lz<=lzend; ++lz) {
                m_max = lzend - lz;
                for (int m=0; m<=m_max; ++m) {
                    rz[1][ly][lz][m+1] =
                        x * rz[0][ly][lz][m+1];
                }
                r[1][ly][lz] = rz[1][ly][lz][0];
            }
        }
        // lx = 2;
        lyend = ltot - 2;
        for (int ly=0; ly<=lyend; ++ly) {
            lzend = lyend - ly;
            for (int lz=0; lz<=lzend; ++lz) {
                m_max = lzend - lz;
                for (int m=0; m<=m_max; ++m) {
                    rz[2][ly][lz][m+1] =
                        x * rz[1][ly][lz][m+1] +
                        rz[0][ly][lz][m+1];
                }
                r[2][ly][lz] = rz[2][ly][lz][0];
            }
        }
        // lx > 2
        for (int lx=3; lx<=ltot; ++lx) {
            lyend = ltot - lx;
            int lxm1 = lx - 1;
            int lxm2 = lx - 2;
            for (int ly=0; ly<=lyend; ++ly) {
                lzend = lyend - ly;
                for (int lz = 0; lz<=lzend; ++lz) {
                    m_max = lzend - lz;
                    for (int m=0; m<=m_max; ++m) {
                        rz[lx][ly][lz][m] =
                            x * rz[lxm1][ly][lz][m+1] +
                            lxm1 * rz[lxm2][ly][lz][m+1];

                    }
                    r[lx][ly][lz] = rz[lx][ly][lz][0];
                }
            }
        }
    }

    constexpr void loop_eval_one_cen(int ltot) noexcept {
        // slower evaluation using loops and rz
        // lz = 1
        int m_max = ltot-1;
        for (int m=0; m<=m_max; ++m) {
            rz[0][0][1][m] = 0.0;
        }
        r[0][0][1] = rz[0][0][1][0];
        // lz = 2
        m_max = ltot-1;
        for (int m=0; m<=m_max; ++m) {
            rz[0][0][2][m] = rz[0][0][0][m+1];
        }
        r[0][0][2] = rz[0][0][2][0];
        // lz > 2
        for (int lz=3; lz<=ltot; ++lz) {
            int lzm1 = lz - 1;
            int lzm2 = lz - 2;
            m_max = ltot - lz;
            for (int m=0; m<=m_max; ++m) {
                rz[0][0][lz][m] = lzm1 * rz[0][0][lzm2][m+1];
            }
            r[0][0][lz] = rz[0][0][lz][0];
        }
        // ly = 1
        int lzend = ltot - 1;
        for (int lz=0; lz<=lzend; ++lz) {
            m_max = lzend - lz;
            for (int m=0; m<=m_max; ++m) {
                rz[0][1][lz][m] = 0.0;
            }
            r[0][1][lz] = rz[0][1][lz][0];
        }
        // ly = 2
        lzend = ltot - 2;
        for (int lz=0; lz<=lzend; ++lz) {
            m_max = ltot - lz;
            for (int m=0; m<=m_max; ++m) {
                rz[0][2][lz][m] = rz[0][0][lz][m+1];
            }
            r[0][2][lz] = rz[0][2][lz][0];
        }
        // ly > 2
        for (int ly=3; ly<=ltot; ++ly) {
            lzend = ltot - ly;
            for (int lz=0; lz<=lzend; ++lz) {
                int lym1 = ly - 1;
                int lym2 = ly - 2;
                m_max = lzend - lz;
                for (int m=0; m<=m_max; ++m) {
                    rz[0][ly][lz][m+1] = lym1 * rz[0][lym2][lz][m+1];
                }
                r[0][ly][lz] = rz[0][ly][lz][0];
            }
        }
        // lx = 1;
        int lyend = ltot - 1;
        for (int ly=0; ly<=lyend; ++ly) {
            lzend = lyend - ly;
            for (int lz=0; lz<=lzend; ++lz) {
                m_max = lzend - lz;
                for (int m=0; m<=m_max; ++m) {
                    rz[1][ly][lz][m+1] = 0.0;
                }
                r[1][ly][lz] = rz[1][ly][lz][0];
            }
        }
        // lx = 2;
        lyend = ltot - 2;
        for (int ly=0; ly<=lyend; ++ly) {
            lzend = lyend - ly;
            for (int lz=0; lz<=lzend; ++lz) {
                m_max = lzend - lz;
                for (int m=0; m<=m_max; ++m) {
                    rz[2][ly][lz][m+1] =
                        rz[0][ly][lz][m+1];
                }
                r[2][ly][lz] = rz[2][ly][lz][0];
            }
        }
        // lx > 2
        for (int lx=3; lx<=ltot; ++lx) {
            lyend = ltot - lx;
            int lxm1 = lx - 1;
            int lxm2 = lx - 2;
            for (int ly=0; ly<=lyend; ++ly) {
                lzend = lyend - ly;
                for (int lz = 0; lz<=lzend; ++lz) {
                    m_max = lzend - lz;
                    for (int m=0; m<=m_max; ++m) {
                        rz[lx][ly][lz][m] =
                            lxm1 * rz[lxm2][ly][lz][m+1];

                    }
                    r[lx][ly][lz] = rz[lx][ly][lz][0];
                }
            }
        }
    }

    void eval_one_cen(double sr,double w,int ltot) noexcept {
        double *r0 = rz[0][0][0];
        for (int l=0; l<=ltot; ++l) r0[l] = 0.5/(0.5 + l);
        double sterm = sr;
        double wterm = -(w+w);
        for (int lx=0; lx<=ltot; ++lx) {
            r0[lx] *= sterm;
            sterm *= wterm;
        }
        for (int lx=0; lx<=ltot; ++lx) {
            int lyend = ltot - lx;
            for (int ly=0; ly<=lyend; ++ly) {
                int lzend = lyend - ly;
                for (int lz=0; lz<=lzend; ++lz) {
                    r[lx][ly][lz] = 0.0;
                }
            }
        }
        for (int lx=0; lx<=ltot; lx+=2) {
            int lyend = ltot - lx;
            int lxh = lx >> 1;
            double dx = dfact[lxh];
            for (int ly=0; ly<=lyend; ly+=2) {
                int lzend = lyend - ly;
                int lyh = ly >> 1;
                double dxy = dfact[lyh] * dx;
                for (int lz=0; lz<=lzend; lz+=2) {
                    int lzh = lz >> 1;
                    int index = lxh + lyh + lzh;
                    r[lx][ly][lz] = dxy * dfact[lzh] * r0[index];
                }
            }
        }
    }

    constexpr void unroll_eval(const double *pq,int ltot) noexcept {
        const double *rm = rz[0][0][0];
        r[0][0][0] = rm[0];
        if ( ltot == 0) return;
        double x1 = pq[0];
        double y1 = pq[1];
        double z1 = pq[2];
        r[0][0][1] = z1*rm[1];
        r[0][1][0] = y1*rm[1];
        r[1][0][0] = x1*rm[1];
        if ( ltot == 1) return;

        double x2= x1*x1;
        double y2= y1*y1;
        double z2= z1*z1;
        r[0][0][2] = z2*rm[2] + rm[1];
        r[0][1][1] = y1*z1*rm[2];
        r[0][2][0] = y2*rm[2] + rm[1];
        r[1][0][1] = x1*z1*rm[2];
        r[1][1][0] = x1*y1*rm[2];
        r[2][0][0] = x2*rm[2] + rm[1];
        if ( ltot == 2) return;

        double x3= x2*x1;
        double y3= y2*y1;
        double z3= z2*z1;
        double z2p = z2*rm[3] + rm[2];
        double y2p = y2*rm[3] + rm[2];
        double x2p = x2*rm[3] + rm[2];
        r[0][0][3] = z3*rm[3] + 3*z1*rm[2];
        r[0][1][2] = y1*z2p;
        r[0][2][1] = z1*y2p;
        r[0][3][0] = y3*rm[3] + 3*y1*rm[2];
        r[1][0][2] = x1*z2p;
        r[1][1][1] = x1*y1*z1*rm[3];
        r[1][2][0] = x1*y2p;
        r[2][0][1] = z1*x2p;
        r[2][1][0] = y1*x2p;
        r[3][0][0] = x3*rm[3] + 3*x1*rm[2];
        if ( ltot == 3) return;

        double x4= x3*x1;
        double y4= y3*y1;
        double z4= z3*z1;
        double x3p = x2*rm[4] + 3. * rm[3];
        double y3p = y2*rm[4] + 3. * rm[3];
        double z3p = z2*rm[4] + 3. * rm[3];
        r[0][0][4] = z2*(z2*rm[4] + 6*rm[3]) + 3*rm[2];
        r[0][1][3] = y1*z1*z3p;
        r[0][2][2] = y2*z2*rm[4] + (z2+y2)*rm[3] + rm[2];
        r[0][3][1] = z1*y1*y3p;
        r[0][4][0] = y4*rm[4] + 6*y2*rm[3] + 3*rm[2];
        r[1][0][3] = x1*z1*z3p;
        r[1][1][2] = x1*y1*(z2*rm[4] + rm[3]);
        r[1][2][1] = x1*z1*(y2*rm[4] + rm[3]);
        r[1][3][0] = x1*y1*y3p;
        r[2][0][2] = x2*z2*rm[4] + (z2+x2)*rm[3] + rm[2];
        r[2][1][1] = y1*z1*(x2*rm[4] + rm[3]);
        r[2][2][0] = x2*y2*rm[4] + (y2+x2)*rm[3] + rm[2];
        r[3][0][1] = x1*z1*x3p;
        r[3][1][0] = x1*y1*x3p;
        r[4][0][0] = x4*rm[4] + 6*x2*rm[3] + 3*rm[2];
        if ( ltot == 4) return;

        double x5= x4*x1;
        double y5= y4*y1;
        double z5= z4*z1;
        r[0][0][5] = z5*rm[5] + 10*z3*rm[4] + 15*z1*rm[3];
        r[0][1][4] = y1*z4*rm[5] + 6*y1*z2*rm[4] + 3*y1*rm[3];
        r[0][2][3] = y2*z3*rm[5] + (z3+3*y2*z1)*rm[4] + 3*z1*rm[3];
        r[0][3][2] = y3*z2*rm[5] + (3*y1*z2+y3)*rm[4] + 3*y1*rm[3];
        r[0][4][1] = y4*z1*rm[5] + 6*y2*z1*rm[4] + 3*z1*rm[3];
        r[0][5][0] = y5*rm[5] + 10*y3*rm[4] + 15*y1*rm[3];
        r[1][0][4] = x1*z4*rm[5] + 6*x1*z2*rm[4] + 3*x1*rm[3];
        r[1][1][3] = x1*y1*z3*rm[5] + 3*x1*y1*z1*rm[4];
        r[1][2][2] = x1*y2*z2*rm[5] + (x1*z2+x1*y2)*rm[4] + x1*rm[3];
        r[1][3][1] = x1*y3*z1*rm[5] + 3*x1*y1*z1*rm[4];
        r[1][4][0] = x1*y4*rm[5] + 6*x1*y2*rm[4] + 3*x1*rm[3];
        r[2][0][3] = x2*z3*rm[5] + (z3+3*x2*z1)*rm[4] + 3*z1*rm[3];
        r[2][1][2] = x2*y1*z2*rm[5] + (y1*z2+x2*y1)*rm[4] + y1*rm[3];
        r[2][2][1] = x2*y2*z1*rm[5] + (y2*z1+x2*z1)*rm[4] + z1*rm[3];
        r[2][3][0] = x2*y3*rm[5] + (y3+3*x2*y1)*rm[4] + 3*y1*rm[3];
        r[3][0][2] = x3*z2*rm[5] + (3*x1*z2+x3)*rm[4] + 3*x1*rm[3];
        r[3][1][1] = x3*y1*z1*rm[5] + 3*x1*y1*z1*rm[4];
        r[3][2][0] = x3*y2*rm[5] + (3*x1*y2+x3)*rm[4] + 3*x1*rm[3];
        r[4][0][1] = x4*z1*rm[5] + 6*x2*z1*rm[4] + 3*z1*rm[3];
        r[4][1][0] = x4*y1*rm[5] + 6*x2*y1*rm[4] + 3*y1*rm[3];
        r[5][0][0] = x5*rm[5] + 10*x3*rm[4] + 15*x1*rm[3];
        if ( ltot == 5) return;

        double x6= x5*x1;
        double y6= y5*y1;
        double z6= z5*z1;
        r[0][0][6] = z6*rm[6] + 15*z4*rm[5] + 45*z2*rm[4] + 15*rm[3];
        r[0][1][5] = y1*z5*rm[6] + 10*y1*z3*rm[5] + 15*y1*z1*rm[4];
        r[0][2][4] = y2*z4*rm[6] + (z4+6*y2*z2)*rm[5] + (6*z2+3*y2)*rm[4] + 3*rm[3];
        r[0][3][3] = y3*z3*rm[6] + (3*y1*z3+3*y3*z1)*rm[5] + 9*y1*z1*rm[4];
        r[0][4][2] = y4*z2*rm[6] + (6*y2*z2+y4)*rm[5] + (3*z2+6*y2)*rm[4] + 3*rm[3];
        r[0][5][1] = y5*z1*rm[6] + 10*y3*z1*rm[5] + 15*y1*z1*rm[4];
        r[0][6][0] = y6*rm[6] + 15*y4*rm[5] + 45*y2*rm[4] + 15*rm[3];
        r[1][0][5] = x1*z5*rm[6] + 10*x1*z3*rm[5] + 15*x1*z1*rm[4];
        r[1][1][4] = x1*y1*z4*rm[6] + 6*x1*y1*z2*rm[5] + 3*x1*y1*rm[4];
        r[1][2][3] = x1*y2*z3*rm[6] + (x1*z3+3*x1*y2*z1)*rm[5] + 3*x1*z1*rm[4];
        r[1][3][2] = x1*y3*z2*rm[6] + (3*x1*y1*z2+x1*y3)*rm[5] + 3*x1*y1*rm[4];
        r[1][4][1] = x1*y4*z1*rm[6] + 6*x1*y2*z1*rm[5] + 3*x1*z1*rm[4];
        r[1][5][0] = x1*y5*rm[6] + 10*x1*y3*rm[5] + 15*x1*y1*rm[4];
        r[2][0][4] = x2*z4*rm[6] + (z4+6*x2*z2)*rm[5] + (6*z2+3*x2)*rm[4] + 3*rm[3];
        r[2][1][3] = x2*y1*z3*rm[6] + (y1*z3+3*x2*y1*z1)*rm[5] + 3*y1*z1*rm[4];
        r[2][2][2] = x2*y2*z2*rm[6] + (y2*z2+x2*z2+x2*y2)*rm[5] + (z2+y2+x2)*rm[4] + rm[3];
        r[2][3][1] = x2*y3*z1*rm[6] + (y3*z1+3*x2*y1*z1)*rm[5] + 3*y1*z1*rm[4];
        r[2][4][0] = x2*y4*rm[6] + (y4+6*x2*y2)*rm[5] + (6*y2+3*x2)*rm[4] + 3*rm[3];
        r[3][0][3] = x3*z3*rm[6] + (3*x1*z3+3*x3*z1)*rm[5] + 9*x1*z1*rm[4];
        r[3][1][2] = x3*y1*z2*rm[6] + (3*x1*y1*z2+x3*y1)*rm[5] + 3*x1*y1*rm[4];
        r[3][2][1] = x3*y2*z1*rm[6] + (3*x1*y2*z1+x3*z1)*rm[5] + 3*x1*z1*rm[4];
        r[3][3][0] = x3*y3*rm[6] + (3*x1*y3+3*x3*y1)*rm[5] + 9*x1*y1*rm[4];
        r[4][0][2] = x4*z2*rm[6] + (6*x2*z2+x4)*rm[5] + (3*z2+6*x2)*rm[4] + 3*rm[3];
        r[4][1][1] = x4*y1*z1*rm[6] + 6*x2*y1*z1*rm[5] + 3*y1*z1*rm[4];
        r[4][2][0] = x4*y2*rm[6] + (6*x2*y2+x4)*rm[5] + (3*y2+6*x2)*rm[4] + 3*rm[3];
        r[5][0][1] = x5*z1*rm[6] + 10*x3*z1*rm[5] + 15*x1*z1*rm[4];
        r[5][1][0] = x5*y1*rm[6] + 10*x3*y1*rm[5] + 15*x1*y1*rm[4];
        r[6][0][0] = x6*rm[6] + 15*x4*rm[5] + 45*x2*rm[4] + 15*rm[3];
        if ( ltot == 6) return;

        double x7= x6*x1;
        double y7= y6*y1;
        double z7= z6*z1;
        r[0][0][7] = z7*rm[7] + 21*z5*rm[6] + 105*z3*rm[5] + 105*z1*rm[4];
        r[0][1][6] = y1*z6*rm[7] + 15*y1*z4*rm[6] + 45*y1*z2*rm[5] + 15*y1*rm[4];
        r[0][2][5] = y2*z5*rm[7] + (z5+10*y2*z3)*rm[6] + (10*z3+15*y2*z1)*rm[5] + 15*z1*rm[4];
        r[0][3][4] = y3*z4*rm[7] + (3*y1*z4+6*y3*z2)*rm[6] + (18*y1*z2+3*y3)*rm[5] + 9*y1*rm[4];
        r[0][4][3] = y4*z3*rm[7] + (6*y2*z3+3*y4*z1)*rm[6] + (3*z3+18*y2*z1)*rm[5] + 9*z1*rm[4];
        r[0][5][2] = y5*z2*rm[7] + (10*y3*z2+y5)*rm[6] + (15*y1*z2+10*y3)*rm[5] + 15*y1*rm[4];
        r[0][6][1] = y6*z1*rm[7] + 15*y4*z1*rm[6] + 45*y2*z1*rm[5] + 15*z1*rm[4];
        r[0][7][0] = y7*rm[7] + 21*y5*rm[6] + 105*y3*rm[5] + 105*y1*rm[4];
        r[1][0][6] = x1*z6*rm[7] + 15*x1*z4*rm[6] + 45*x1*z2*rm[5] + 15*x1*rm[4];
        r[1][1][5] = x1*y1*z5*rm[7] + 10*x1*y1*z3*rm[6] + 15*x1*y1*z1*rm[5];
        r[1][2][4] = x1*y2*z4*rm[7] + (x1*z4+6*x1*y2*z2)*rm[6] + (6*x1*z2+3*x1*y2)*rm[5] + 3*x1*rm[4];
        r[1][3][3] = x1*y3*z3*rm[7] + (3*x1*y1*z3+3*x1*y3*z1)*rm[6] + 9*x1*y1*z1*rm[5];
        r[1][4][2] = x1*y4*z2*rm[7] + (6*x1*y2*z2+x1*y4)*rm[6] + (3*x1*z2+6*x1*y2)*rm[5] + 3*x1*rm[4];
        r[1][5][1] = x1*y5*z1*rm[7] + 10*x1*y3*z1*rm[6] + 15*x1*y1*z1*rm[5];
        r[1][6][0] = x1*y6*rm[7] + 15*x1*y4*rm[6] + 45*x1*y2*rm[5] + 15*x1*rm[4];
        r[2][0][5] = x2*z5*rm[7] + (z5+10*x2*z3)*rm[6] + (10*z3+15*x2*z1)*rm[5] + 15*z1*rm[4];
        r[2][1][4] = x2*y1*z4*rm[7] + (y1*z4+6*x2*y1*z2)*rm[6] + (6*y1*z2+3*x2*y1)*rm[5] + 3*y1*rm[4];
        r[2][2][3] = x2*y2*z3*rm[7] + (y2*z3+x2*z3+3*x2*y2*z1)*rm[6] + (z3+3*y2*z1+3*x2*z1)*rm[5] +
                     3*z1*rm[4];
        r[2][3][2] = x2*y3*z2*rm[7] + (y3*z2+3*x2*y1*z2+x2*y3)*rm[6] + (3*y1*z2+y3+3*x2*y1)*rm[5] +
                     3*y1*rm[4];
        r[2][4][1] = x2*y4*z1*rm[7] + (y4*z1+6*x2*y2*z1)*rm[6] + (6*y2*z1+3*x2*z1)*rm[5] + 3*z1*rm[4];
        r[2][5][0] = x2*y5*rm[7] + (y5+10*x2*y3)*rm[6] + (10*y3+15*x2*y1)*rm[5] + 15*y1*rm[4];
        r[3][0][4] = x3*z4*rm[7] + (3*x1*z4+6*x3*z2)*rm[6] + (18*x1*z2+3*x3)*rm[5] + 9*x1*rm[4];
        r[3][1][3] = x3*y1*z3*rm[7] + (3*x1*y1*z3+3*x3*y1*z1)*rm[6] + 9*x1*y1*z1*rm[5];
        r[3][2][2] = x3*y2*z2*rm[7] + (3*x1*y2*z2+x3*z2+x3*y2)*rm[6] + (3*x1*z2+3*x1*y2+x3)*rm[5] +
                     3*x1*rm[4];
        r[3][3][1] = x3*y3*z1*rm[7] + (3*x1*y3*z1+3*x3*y1*z1)*rm[6] + 9*x1*y1*z1*rm[5];
        r[3][4][0] = x3*y4*rm[7] + (3*x1*y4+6*x3*y2)*rm[6] + (18*x1*y2+3*x3)*rm[5] + 9*x1*rm[4];
        r[4][0][3] = x4*z3*rm[7] + (6*x2*z3+3*x4*z1)*rm[6] + (3*z3+18*x2*z1)*rm[5] + 9*z1*rm[4];
        r[4][1][2] = x4*y1*z2*rm[7] + (6*x2*y1*z2+x4*y1)*rm[6] + (3*y1*z2+6*x2*y1)*rm[5] + 3*y1*rm[4];
        r[4][2][1] = x4*y2*z1*rm[7] + (6*x2*y2*z1+x4*z1)*rm[6] + (3*y2*z1+6*x2*z1)*rm[5] + 3*z1*rm[4];
        r[4][3][0] = x4*y3*rm[7] + (6*x2*y3+3*x4*y1)*rm[6] + (3*y3+18*x2*y1)*rm[5] + 9*y1*rm[4];
        r[5][0][2] = x5*z2*rm[7] + (10*x3*z2+x5)*rm[6] + (15*x1*z2+10*x3)*rm[5] + 15*x1*rm[4];
        r[5][1][1] = x5*y1*z1*rm[7] + 10*x3*y1*z1*rm[6] + 15*x1*y1*z1*rm[5];
        r[5][2][0] = x5*y2*rm[7] + (10*x3*y2+x5)*rm[6] + (15*x1*y2+10*x3)*rm[5] + 15*x1*rm[4];
        r[6][0][1] = x6*z1*rm[7] + 15*x4*z1*rm[6] + 45*x2*z1*rm[5] + 15*z1*rm[4];
        r[6][1][0] = x6*y1*rm[7] + 15*x4*y1*rm[6] + 45*x2*y1*rm[5] + 15*y1*rm[4];
        r[7][0][0] = x7*rm[7] + 21*x5*rm[6] + 105*x3*rm[5] + 105*x1*rm[4];
        if ( ltot == 7) return;

        double x8= x7*x1;
        double y8= y7*y1;
        double z8= z7*z1;
        r[0][0][8] = z8*rm[8] + 28*z6*rm[7] + 210*z4*rm[6] + 420*z2*rm[5] + 105*rm[4];
        r[0][1][7] = y1*z7*rm[8] + 21*y1*z5*rm[7] + 105*y1*z3*rm[6] + 105*y1*z1*rm[5];
        r[0][2][6] = y2*z6*rm[8] + (z6+15*y2*z4)*rm[7] + (15*z4+45*y2*z2)*rm[6] +
                     (45*z2+15*y2)*rm[5] + 15*rm[4];
        r[0][3][5] = y3*z5*rm[8] + (3*y1*z5+10*y3*z3)*rm[7] + (30*y1*z3+15*y3*z1)*rm[6] + 45*y1*z1*rm[5];
        r[0][4][4] = y4*z4*rm[8] + (6*y2*z4+6*y4*z2)*rm[7] + (3*z4+36*y2*z2+3*y4)*rm[6] +
                     (18*z2+18*y2)*rm[5] + 9*rm[4];
        r[0][5][3] = y5*z3*rm[8] + (10*y3*z3+3*y5*z1)*rm[7] + (15*y1*z3+30*y3*z1)*rm[6] + 45*y1*z1*rm[5];
        r[0][6][2] = y6*z2*rm[8] + (15*y4*z2+y6)*rm[7] + (45*y2*z2+15*y4)*rm[6] +
                     (15*z2+45*y2)*rm[5] + 15*rm[4];
        r[0][7][1] = y7*z1*rm[8] + 21*y5*z1*rm[7] + 105*y3*z1*rm[6] + 105*y1*z1*rm[5];
        r[0][8][0] = y8*rm[8] + 28*y6*rm[7] + 210*y4*rm[6] + 420*y2*rm[5] + 105*rm[4];
        r[1][0][7] = x1*z7*rm[8] + 21*x1*z5*rm[7] + 105*x1*z3*rm[6] + 105*x1*z1*rm[5];
        r[1][1][6] = x1*y1*z6*rm[8] + 15*x1*y1*z4*rm[7] + 45*x1*y1*z2*rm[6] + 15*x1*y1*rm[5];
        r[1][2][5] = x1*y2*z5*rm[8] + (x1*z5+10*x1*y2*z3)*rm[7] + (10*x1*z3+15*x1*y2*z1)*rm[6] +
                     15*x1*z1*rm[5];
        r[1][3][4] = x1*y3*z4*rm[8] + (3*x1*y1*z4+6*x1*y3*z2)*rm[7] + (18*x1*y1*z2+3*x1*y3)*rm[6] +
                     9*x1*y1*rm[5];
        r[1][4][3] = x1*y4*z3*rm[8] + (6*x1*y2*z3+3*x1*y4*z1)*rm[7] + (3*x1*z3+18*x1*y2*z1)*rm[6] +
                     9*x1*z1*rm[5];
        r[1][5][2] = x1*y5*z2*rm[8] + (10*x1*y3*z2+x1*y5)*rm[7] + (15*x1*y1*z2+10*x1*y3)*rm[6] +
                     15*x1*y1*rm[5];
        r[1][6][1] = x1*y6*z1*rm[8] + 15*x1*y4*z1*rm[7] + 45*x1*y2*z1*rm[6] + 15*x1*z1*rm[5];
        r[1][7][0] = x1*y7*rm[8] + 21*x1*y5*rm[7] + 105*x1*y3*rm[6] + 105*x1*y1*rm[5];
        r[2][0][6] = x2*z6*rm[8] + (z6+15*x2*z4)*rm[7] + (15*z4+45*x2*z2)*rm[6] +
                     (45*z2+15*x2)*rm[5] + 15*rm[4];
        r[2][1][5] = x2*y1*z5*rm[8] + (y1*z5+10*x2*y1*z3)*rm[7] + (10*y1*z3+15*x2*y1*z1)*rm[6] +
                     15*y1*z1*rm[5];
        r[2][2][4] = x2*y2*z4*rm[8] + (y2*z4+x2*z4+6*x2*y2*z2)*rm[7] + (z4+6*y2*z2+6*x2*z2+3*x2*y2)*rm[6]
                     + (6*z2+3*y2+3*x2)*rm[5] + 3*rm[4];
        r[2][3][3] = x2*y3*z3*rm[8] + (y3*z3+3*x2*y1*z3+3*x2*y3*z1)*rm[7] + (3*y1*z3+3*y3*z1+9*x2*y1*z1)
                     *rm[6] + 9*y1*z1*rm[5];
        r[2][4][2] = x2*y4*z2*rm[8] + (y4*z2+6*x2*y2*z2+x2*y4)*rm[7] + (6*y2*z2+3*x2*z2+y4+6*x2*y2)*rm[6]
                     + (3*z2+6*y2+3*x2)*rm[5] + 3*rm[4];
        r[2][5][1] = x2*y5*z1*rm[8] + (y5*z1+10*x2*y3*z1)*rm[7] + (10*y3*z1+15*x2*y1*z1)*rm[6] +
                     15*y1*z1*rm[5];
        r[2][6][0] = x2*y6*rm[8] + (y6+15*x2*y4)*rm[7] + (15*y4+45*x2*y2)*rm[6] +
                     (45*y2+15*x2)*rm[5] + 15*rm[4];
        r[3][0][5] = x3*z5*rm[8] + (3*x1*z5+10*x3*z3)*rm[7] + (30*x1*z3+15*x3*z1)*rm[6] + 45*x1*z1*rm[5];
        r[3][1][4] = x3*y1*z4*rm[8] + (3*x1*y1*z4+6*x3*y1*z2)*rm[7] + (18*x1*y1*z2+3*x3*y1)*rm[6] +
                     9*x1*y1*rm[5];
        r[3][2][3] = x3*y2*z3*rm[8] + (3*x1*y2*z3+x3*z3+3*x3*y2*z1)*rm[7] + (3*x1*z3+9*x1*y2*z1+3*x3*z1)
                     *rm[6] + 9*x1*z1*rm[5];
        r[3][3][2] = x3*y3*z2*rm[8] + (3*x1*y3*z2+3*x3*y1*z2+x3*y3)*rm[7] + (9*x1*y1*z2+3*x1*y3+3*x3*y1)
                     *rm[6] + 9*x1*y1*rm[5];
        r[3][4][1] = x3*y4*z1*rm[8] + (3*x1*y4*z1+6*x3*y2*z1)*rm[7] + (18*x1*y2*z1+3*x3*z1)*rm[6] +
                     9*x1*z1*rm[5];
        r[3][5][0] = x3*y5*rm[8] + (3*x1*y5+10*x3*y3)*rm[7] + (30*x1*y3+15*x3*y1)*rm[6] + 45*x1*y1*rm[5];
        r[4][0][4] = x4*z4*rm[8] + (6*x2*z4+6*x4*z2)*rm[7] + (3*z4+36*x2*z2+3*x4)*rm[6] +
                     (18*z2+18*x2)*rm[5] + 9*rm[4];
        r[4][1][3] = x4*y1*z3*rm[8] + (6*x2*y1*z3+3*x4*y1*z1)*rm[7] + (3*y1*z3+18*x2*y1*z1)*rm[6] +
                     9*y1*z1*rm[5];
        r[4][2][2] = x4*y2*z2*rm[8] + (6*x2*y2*z2+x4*z2+x4*y2)*rm[7] + (3*y2*z2+6*x2*z2+6*x2*y2+x4)*rm[6]
                     + (3*z2+3*y2+6*x2)*rm[5] + 3*rm[4];
        r[4][3][1] = x4*y3*z1*rm[8] + (6*x2*y3*z1+3*x4*y1*z1)*rm[7] + (3*y3*z1+18*x2*y1*z1)*rm[6] +
                     9*y1*z1*rm[5];
        r[4][4][0] = x4*y4*rm[8] + (6*x2*y4+6*x4*y2)*rm[7] + (3*y4+36*x2*y2+3*x4)*rm[6] +
                     (18*y2+18*x2)*rm[5] + 9*rm[4];
        r[5][0][3] = x5*z3*rm[8] + (10*x3*z3+3*x5*z1)*rm[7] + (15*x1*z3+30*x3*z1)*rm[6] + 45*x1*z1*rm[5];
        r[5][1][2] = x5*y1*z2*rm[8] + (10*x3*y1*z2+x5*y1)*rm[7] + (15*x1*y1*z2+10*x3*y1)*rm[6] +
                     15*x1*y1*rm[5];
        r[5][2][1] = x5*y2*z1*rm[8] + (10*x3*y2*z1+x5*z1)*rm[7] + (15*x1*y2*z1+10*x3*z1)*rm[6] +
                     15*x1*z1*rm[5];
        r[5][3][0] = x5*y3*rm[8] + (10*x3*y3+3*x5*y1)*rm[7] + (15*x1*y3+30*x3*y1)*rm[6] + 45*x1*y1*rm[5];
        r[6][0][2] = x6*z2*rm[8] + (15*x4*z2+x6)*rm[7] + (45*x2*z2+15*x4)*rm[6] +
                     (15*z2+45*x2)*rm[5] + 15*rm[4];
        r[6][1][1] = x6*y1*z1*rm[8] + 15*x4*y1*z1*rm[7] + 45*x2*y1*z1*rm[6] + 15*y1*z1*rm[5];
        r[6][2][0] = x6*y2*rm[8] + (15*x4*y2+x6)*rm[7] + (45*x2*y2+15*x4)*rm[6] +
                     (15*y2+45*x2)*rm[5] + 15*rm[4];
        r[7][0][1] = x7*z1*rm[8] + 21*x5*z1*rm[7] + 105*x3*z1*rm[6] + 105*x1*z1*rm[5];
        r[7][1][0] = x7*y1*rm[8] + 21*x5*y1*rm[7] + 105*x3*y1*rm[6] + 105*x1*y1*rm[5];
        r[8][0][0] = x8*rm[8] + 28*x6*rm[7] + 210*x4*rm[6] + 420*x2*rm[5] + 105*rm[4];
        if ( ltot == 8) return;

        double x9= x8*x1;
        double y9= y8*y1;
        double z9= z8*z1;
        r[0][0][9] = z9*rm[9] + 36*z7*rm[8] + 378*z5*rm[7] + 1260*z3*rm[6] + 945*z1*rm[5];
        r[0][1][8] = y1*z8*rm[9] + 28*y1*z6*rm[8] + 210*y1*z4*rm[7] + 420*y1*z2*rm[6] + 105*y1*rm[5];
        r[0][2][7] = y2*z7*rm[9] + (z7+21*y2*z5)*rm[8] + (21*z5+105*y2*z3)*rm[7] +
                     (105*z3+105*y2*z1)*rm[6] + 105*z1*rm[5];
        r[0][3][6] = y3*z6*rm[9] + (3*y1*z6+15*y3*z4)*rm[8] + (45*y1*z4+45*y3*z2)*rm[7] +
                     (135*y1*z2+15*y3)*rm[6] + 45*y1*rm[5];
        r[0][4][5] = y4*z5*rm[9] + (6*y2*z5+10*y4*z3)*rm[8] + (3*z5+60*y2*z3+15*y4*z1)*rm[7] +
                     (30*z3+90*y2*z1)*rm[6] + 45*z1*rm[5];
        r[0][5][4] = y5*z4*rm[9] + (10*y3*z4+6*y5*z2)*rm[8] + (15*y1*z4+60*y3*z2+3*y5)*rm[7] +
                     (90*y1*z2+30*y3)*rm[6] + 45*y1*rm[5];
        r[0][6][3] = y6*z3*rm[9] + (15*y4*z3+3*y6*z1)*rm[8] + (45*y2*z3+45*y4*z1)*rm[7] +
                     (15*z3+135*y2*z1)*rm[6] + 45*z1*rm[5];
        r[0][7][2] = y7*z2*rm[9] + (21*y5*z2+y7)*rm[8] + (105*y3*z2+21*y5)*rm[7] +
                     (105*y1*z2+105*y3)*rm[6] + 105*y1*rm[5];
        r[0][8][1] = y8*z1*rm[9] + 28*y6*z1*rm[8] + 210*y4*z1*rm[7] + 420*y2*z1*rm[6] + 105*z1*rm[5];
        r[0][9][0] = y9*rm[9] + 36*y7*rm[8] + 378*y5*rm[7] + 1260*y3*rm[6] + 945*y1*rm[5];
        r[1][0][8] = x1*z8*rm[9] + 28*x1*z6*rm[8] + 210*x1*z4*rm[7] + 420*x1*z2*rm[6] + 105*x1*rm[5];
        r[1][1][7] = x1*y1*z7*rm[9] + 21*x1*y1*z5*rm[8] + 105*x1*y1*z3*rm[7] + 105*x1*y1*z1*rm[6];
        r[1][2][6] = x1*y2*z6*rm[9] + (x1*z6+15*x1*y2*z4)*rm[8] + (15*x1*z4+45*x1*y2*z2)*rm[7] +
                     (45*x1*z2+15*x1*y2)*rm[6] + 15*x1*rm[5];
        r[1][3][5] = x1*y3*z5*rm[9] + (3*x1*y1*z5+10*x1*y3*z3)*rm[8] + (30*x1*y1*z3+15*x1*y3*z1)*rm[7] +
                     45*x1*y1*z1*rm[6];
        r[1][4][4] = x1*y4*z4*rm[9] + (6*x1*y2*z4+6*x1*y4*z2)*rm[8] + (3*x1*z4+36*x1*y2*z2+3*x1*y4)*rm[7] +
                     (18*x1*z2+18*x1*y2)*rm[6] + 9*x1*rm[5];
        r[1][5][3] = x1*y5*z3*rm[9] + (10*x1*y3*z3+3*x1*y5*z1)*rm[8] + (15*x1*y1*z3+30*x1*y3*z1)*rm[7] +
                     45*x1*y1*z1*rm[6];
        r[1][6][2] = x1*y6*z2*rm[9] + (15*x1*y4*z2+x1*y6)*rm[8] + (45*x1*y2*z2+15*x1*y4)*rm[7] +
                     (15*x1*z2+45*x1*y2)*rm[6] + 15*x1*rm[5];
        r[1][7][1] = x1*y7*z1*rm[9] + 21*x1*y5*z1*rm[8] + 105*x1*y3*z1*rm[7] + 105*x1*y1*z1*rm[6];
        r[1][8][0] = x1*y8*rm[9] + 28*x1*y6*rm[8] + 210*x1*y4*rm[7] + 420*x1*y2*rm[6] + 105*x1*rm[5];
        r[2][0][7] = x2*z7*rm[9] + (z7+21*x2*z5)*rm[8] + (21*z5+105*x2*z3)*rm[7] +
                     (105*z3+105*x2*z1)*rm[6] + 105*z1*rm[5];
        r[2][1][6] = x2*y1*z6*rm[9] + (y1*z6+15*x2*y1*z4)*rm[8] + (15*y1*z4+45*x2*y1*z2)*rm[7] +
                     (45*y1*z2+15*x2*y1)*rm[6] + 15*y1*rm[5];
        r[2][2][5] = x2*y2*z5*rm[9] + (y2*z5+x2*z5+10*x2*y2*z3)*rm[8] + (z5+10*y2*z3+10*x2*z3+15*x2*y2*z1)
                     *rm[7] +
                     (10*z3+15*y2*z1+15*x2*z1)*rm[6] + 15*z1*rm[5];
        r[2][3][4] = x2*y3*z4*rm[9] + (y3*z4+3*x2*y1*z4+6*x2*y3*z2)*rm[8] + (3*y1*z4+6*y3*z2+18*x2*y1*z2
                     +3*x2*y3)*rm[7] +
                     (18*y1*z2+3*y3+9*x2*y1)*rm[6] + 9*y1*rm[5];
        r[2][4][3] = x2*y4*z3*rm[9] + (y4*z3+6*x2*y2*z3+3*x2*y4*z1)*rm[8] + (6*y2*z3+3*x2*z3+3*y4*z1
                     +18*x2*y2*z1)*rm[7] +
                     (3*z3+18*y2*z1+9*x2*z1)*rm[6] + 9*z1*rm[5];
        r[2][5][2] = x2*y5*z2*rm[9] + (y5*z2+10*x2*y3*z2+x2*y5)*rm[8] + (10*y3*z2+15*x2*y1*z2+y5+10*x2*y3)
                     *rm[7] +
                     (15*y1*z2+10*y3+15*x2*y1)*rm[6] + 15*y1*rm[5];
        r[2][6][1] = x2*y6*z1*rm[9] + (y6*z1+15*x2*y4*z1)*rm[8] + (15*y4*z1+45*x2*y2*z1)*rm[7] +
                     (45*y2*z1+15*x2*z1)*rm[6] + 15*z1*rm[5];
        r[2][7][0] = x2*y7*rm[9] + (y7+21*x2*y5)*rm[8] + (21*y5+105*x2*y3)*rm[7] +
                     (105*y3+105*x2*y1)*rm[6] + 105*y1*rm[5];
        r[3][0][6] = x3*z6*rm[9] + (3*x1*z6+15*x3*z4)*rm[8] + (45*x1*z4+45*x3*z2)*rm[7] +
                     (135*x1*z2+15*x3)*rm[6] + 45*x1*rm[5];
        r[3][1][5] = x3*y1*z5*rm[9] + (3*x1*y1*z5+10*x3*y1*z3)*rm[8] + (30*x1*y1*z3+15*x3*y1*z1)*rm[7] +
                     45*x1*y1*z1*rm[6];
        r[3][2][4] = x3*y2*z4*rm[9] + (3*x1*y2*z4+x3*z4+6*x3*y2*z2)*rm[8] + (3*x1*z4+18*x1*y2*z2+6*x3*z2
                     +3*x3*y2)*rm[7] +
                     (18*x1*z2+9*x1*y2+3*x3)*rm[6] + 9*x1*rm[5];
        r[3][3][3] = x3*y3*z3*rm[9] + (3*x1*y3*z3+3*x3*y1*z3+3*x3*y3*z1)*rm[8] +
                     (9*x1*y1*z3+9*x1*y3*z1+9*x3*y1*z1)*rm[7] +
                     27*x1*y1*z1*rm[6];
        r[3][4][2] = x3*y4*z2*rm[9] + (3*x1*y4*z2+6*x3*y2*z2+x3*y4)*rm[8] + (18*x1*y2*z2+3*x3*z2+3*x1*y4
                     +6*x3*y2)*rm[7] +
                     (9*x1*z2+18*x1*y2+3*x3)*rm[6] + 9*x1*rm[5];
        r[3][5][1] = x3*y5*z1*rm[9] + (3*x1*y5*z1+10*x3*y3*z1)*rm[8] + (30*x1*y3*z1+15*x3*y1*z1)*rm[7] +
                     45*x1*y1*z1*rm[6];
        r[3][6][0] = x3*y6*rm[9] + (3*x1*y6+15*x3*y4)*rm[8] + (45*x1*y4+45*x3*y2)*rm[7] +
                     (135*x1*y2+15*x3)*rm[6] + 45*x1*rm[5];
        r[4][0][5] = x4*z5*rm[9] + (6*x2*z5+10*x4*z3)*rm[8] + (3*z5+60*x2*z3+15*x4*z1)*rm[7] +
                     (30*z3+90*x2*z1)*rm[6] + 45*z1*rm[5];
        r[4][1][4] = x4*y1*z4*rm[9] + (6*x2*y1*z4+6*x4*y1*z2)*rm[8] + (3*y1*z4+36*x2*y1*z2+3*x4*y1)*rm[7] +
                     (18*y1*z2+18*x2*y1)*rm[6] + 9*y1*rm[5];
        r[4][2][3] = x4*y2*z3*rm[9] + (6*x2*y2*z3+x4*z3+3*x4*y2*z1)*rm[8] + (3*y2*z3+6*x2*z3+18*x2*y2*z1
                     +3*x4*z1)*rm[7] +
                     (3*z3+9*y2*z1+18*x2*z1)*rm[6] + 9*z1*rm[5];
        r[4][3][2] = x4*y3*z2*rm[9] + (6*x2*y3*z2+3*x4*y1*z2+x4*y3)*rm[8] + (3*y3*z2+18*x2*y1*z2+6*x2*y3
                     +3*x4*y1)*rm[7] +
                     (9*y1*z2+3*y3+18*x2*y1)*rm[6] + 9*y1*rm[5];
        r[4][4][1] = x4*y4*z1*rm[9] + (6*x2*y4*z1+6*x4*y2*z1)*rm[8] + (3*y4*z1+36*x2*y2*z1+3*x4*z1)*rm[7] +
                     (18*y2*z1+18*x2*z1)*rm[6] + 9*z1*rm[5];
        r[4][5][0] = x4*y5*rm[9] + (6*x2*y5+10*x4*y3)*rm[8] + (3*y5+60*x2*y3+15*x4*y1)*rm[7] +
                     (30*y3+90*x2*y1)*rm[6] + 45*y1*rm[5];
        r[5][0][4] = x5*z4*rm[9] + (10*x3*z4+6*x5*z2)*rm[8] + (15*x1*z4+60*x3*z2+3*x5)*rm[7] +
                     (90*x1*z2+30*x3)*rm[6] + 45*x1*rm[5];
        r[5][1][3] = x5*y1*z3*rm[9] + (10*x3*y1*z3+3*x5*y1*z1)*rm[8] + (15*x1*y1*z3+30*x3*y1*z1)*rm[7] +
                     45*x1*y1*z1*rm[6];
        r[5][2][2] = x5*y2*z2*rm[9] + (10*x3*y2*z2+x5*z2+x5*y2)*rm[8] + (15*x1*y2*z2+10*x3*z2+10*x3*y2+x5)
                     *rm[7] +
                     (15*x1*z2+15*x1*y2+10*x3)*rm[6] + 15*x1*rm[5];
        r[5][3][1] = x5*y3*z1*rm[9] + (10*x3*y3*z1+3*x5*y1*z1)*rm[8] + (15*x1*y3*z1+30*x3*y1*z1)*rm[7] +
                     45*x1*y1*z1*rm[6];
        r[5][4][0] = x5*y4*rm[9] + (10*x3*y4+6*x5*y2)*rm[8] + (15*x1*y4+60*x3*y2+3*x5)*rm[7] +
                     (90*x1*y2+30*x3)*rm[6] + 45*x1*rm[5];
        r[6][0][3] = x6*z3*rm[9] + (15*x4*z3+3*x6*z1)*rm[8] + (45*x2*z3+45*x4*z1)*rm[7] +
                     (15*z3+135*x2*z1)*rm[6] + 45*z1*rm[5];
        r[6][1][2] = x6*y1*z2*rm[9] + (15*x4*y1*z2+x6*y1)*rm[8] + (45*x2*y1*z2+15*x4*y1)*rm[7] +
                     (15*y1*z2+45*x2*y1)*rm[6] + 15*y1*rm[5];
        r[6][2][1] = x6*y2*z1*rm[9] + (15*x4*y2*z1+x6*z1)*rm[8] + (45*x2*y2*z1+15*x4*z1)*rm[7] +
                     (15*y2*z1+45*x2*z1)*rm[6] + 15*z1*rm[5];
        r[6][3][0] = x6*y3*rm[9] + (15*x4*y3+3*x6*y1)*rm[8] + (45*x2*y3+45*x4*y1)*rm[7] +
                     (15*y3+135*x2*y1)*rm[6] + 45*y1*rm[5];
        r[7][0][2] = x7*z2*rm[9] + (21*x5*z2+x7)*rm[8] + (105*x3*z2+21*x5)*rm[7] +
                     (105*x1*z2+105*x3)*rm[6] + 105*x1*rm[5];
        r[7][1][1] = x7*y1*z1*rm[9] + 21*x5*y1*z1*rm[8] + 105*x3*y1*z1*rm[7] + 105*x1*y1*z1*rm[6];
        r[7][2][0] = x7*y2*rm[9] + (21*x5*y2+x7)*rm[8] + (105*x3*y2+21*x5)*rm[7] +
                     (105*x1*y2+105*x3)*rm[6] + 105*x1*rm[5];
        r[8][0][1] = x8*z1*rm[9] + 28*x6*z1*rm[8] + 210*x4*z1*rm[7] + 420*x2*z1*rm[6] + 105*z1*rm[5];
        r[8][1][0] = x8*y1*rm[9] + 28*x6*y1*rm[8] + 210*x4*y1*rm[7] + 420*x2*y1*rm[6] + 105*y1*rm[5];
        r[9][0][0] = x9*rm[9] + 36*x7*rm[8] + 378*x5*rm[7] + 1260*x3*rm[6] + 945*x1*rm[5];
        if ( ltot == 9) return;

        double x10= x9*x1;
        double y10= y9*y1;
        double z10= z9*z1;
        r[0][0][10] = z10*rm[10] + 45*z8*rm[9] + 630*z6*rm[8] + 3150*z4*rm[7] + 4725*z2*rm[6] + 945*rm[5];
        r[0][1][9] = y1*z9*rm[10] + 36*y1*z7*rm[9] + 378*y1*z5*rm[8] + 1260*y1*z3*rm[7] + 945*y1*z1*rm[6];
        r[0][2][8] = y2*z8*rm[10] + (z8+28*y2*z6)*rm[9] + (28*z6+210*y2*z4)*rm[8] +
                     (210*z4+420*y2*z2)*rm[7] +
                     (420*z2+105*y2)*rm[6] + 105*rm[5];
        r[0][3][7] = y3*z7*rm[10] + (3*y1*z7+21*y3*z5)*rm[9] + (63*y1*z5+105*y3*z3)*rm[8] +
                     (315*y1*z3+105*y3*z1)*rm[7] + 315*y1*z1*rm[6];
        r[0][4][6] = y4*z6*rm[10] + (6*y2*z6+15*y4*z4)*rm[9] + (3*z6+90*y2*z4+45*y4*z2)*rm[8] +
                     (45*z4+270*y2*z2+15*y4)*rm[7] + (135*z2+90*y2)*rm[6] + 45*rm[5];
        r[0][5][5] = y5*z5*rm[10] + (10*y3*z5+10*y5*z3)*rm[9] + (15*y1*z5+100*y3*z3+15*y5*z1)*rm[8] +
                     (150*y1*z3+150*y3*z1)*rm[7] + 225*y1*z1*rm[6];
        r[0][6][4] = y6*z4*rm[10] + (15*y4*z4+6*y6*z2)*rm[9] + (45*y2*z4+90*y4*z2+3*y6)*rm[8] +
                     (15*z4+270*y2*z2+45*y4)*rm[7] + (90*z2+135*y2)*rm[6] + 45*rm[5];
        r[0][7][3] = y7*z3*rm[10] + (21*y5*z3+3*y7*z1)*rm[9] + (105*y3*z3+63*y5*z1)*rm[8] +
                     (105*y1*z3+315*y3*z1)*rm[7] + 315*y1*z1*rm[6];
        r[0][8][2] = y8*z2*rm[10] + (28*y6*z2+y8)*rm[9] + (210*y4*z2+28*y6)*rm[8] +
                     (420*y2*z2+210*y4)*rm[7] +
                     (105*z2+420*y2)*rm[6] + 105*rm[5];
        r[0][9][1] = y9*z1*rm[10] + 36*y7*z1*rm[9] + 378*y5*z1*rm[8] + 1260*y3*z1*rm[7] + 945*y1*z1*rm[6];
        r[0][10][0] = y10*rm[10] + 45*y8*rm[9] + 630*y6*rm[8] + 3150*y4*rm[7] + 4725*y2*rm[6] + 945*rm[5];
        r[1][0][9] = x1*z9*rm[10] + 36*x1*z7*rm[9] + 378*x1*z5*rm[8] + 1260*x1*z3*rm[7] + 945*x1*z1*rm[6];
        r[1][1][8] = x1*y1*z8*rm[10] + 28*x1*y1*z6*rm[9] + 210*x1*y1*z4*rm[8] + 420*x1*y1*z2*rm[7] +
                     105*x1*y1*rm[6];
        r[1][2][7] = x1*y2*z7*rm[10] + (x1*z7+21*x1*y2*z5)*rm[9] + (21*x1*z5+105*x1*y2*z3)*rm[8] +
                     (105*x1*z3+105*x1*y2*z1)*rm[7] + 105*x1*z1*rm[6];
        r[1][3][6] = x1*y3*z6*rm[10] + (3*x1*y1*z6+15*x1*y3*z4)*rm[9] + (45*x1*y1*z4+45*x1*y3*z2)*rm[8] +
                     (135*x1*y1*z2+15*x1*y3)*rm[7] + 45*x1*y1*rm[6];
        r[1][4][5] = x1*y4*z5*rm[10] + (6*x1*y2*z5+10*x1*y4*z3)*rm[9] + (3*x1*z5+60*x1*y2*z3+15*x1*y4*z1)
                     *rm[8] +
                     (30*x1*z3+90*x1*y2*z1)*rm[7] + 45*x1*z1*rm[6];
        r[1][5][4] = x1*y5*z4*rm[10] + (10*x1*y3*z4+6*x1*y5*z2)*rm[9] + (15*x1*y1*z4+60*x1*y3*z2+3*x1*y5)
                     *rm[8] +
                     (90*x1*y1*z2+30*x1*y3)*rm[7] + 45*x1*y1*rm[6];
        r[1][6][3] = x1*y6*z3*rm[10] + (15*x1*y4*z3+3*x1*y6*z1)*rm[9] + (45*x1*y2*z3+45*x1*y4*z1)*rm[8] +
                     (15*x1*z3+135*x1*y2*z1)*rm[7] + 45*x1*z1*rm[6];
        r[1][7][2] = x1*y7*z2*rm[10] + (21*x1*y5*z2+x1*y7)*rm[9] + (105*x1*y3*z2+21*x1*y5)*rm[8] +
                     (105*x1*y1*z2+105*x1*y3)*rm[7] + 105*x1*y1*rm[6];
        r[1][8][1] = x1*y8*z1*rm[10] + 28*x1*y6*z1*rm[9] + 210*x1*y4*z1*rm[8] + 420*x1*y2*z1*rm[7] +
                     105*x1*z1*rm[6];
        r[1][9][0] = x1*y9*rm[10] + 36*x1*y7*rm[9] + 378*x1*y5*rm[8] + 1260*x1*y3*rm[7] + 945*x1*y1*rm[6];
        r[2][0][8] = x2*z8*rm[10] + (z8+28*x2*z6)*rm[9] + (28*z6+210*x2*z4)*rm[8] +
                     (210*z4+420*x2*z2)*rm[7] +
                     (420*z2+105*x2)*rm[6] + 105*rm[5];
        r[2][1][7] = x2*y1*z7*rm[10] + (y1*z7+21*x2*y1*z5)*rm[9] + (21*y1*z5+105*x2*y1*z3)*rm[8] +
                     (105*y1*z3+105*x2*y1*z1)*rm[7] + 105*y1*z1*rm[6];
        r[2][2][6] = x2*y2*z6*rm[10] + (y2*z6+x2*z6+15*x2*y2*z4)*rm[9] + (z6+15*y2*z4+15*x2*z4+45*x2*y2*z2)
                     *rm[8] +
                     (15*z4+45*y2*z2+45*x2*z2+15*x2*y2)*rm[7] + (45*z2+15*y2+15*x2)*rm[6] + 15*rm[5];
        r[2][3][5] = x2*y3*z5*rm[10] + (y3*z5+3*x2*y1*z5+10*x2*y3*z3)*rm[9] +
                     (3*y1*z5+10*y3*z3+30*x2*y1*z3+15*x2*y3*z1)*rm[8]
                     + (30*y1*z3+15*y3*z1+45*x2*y1*z1)*rm[7] + 45*y1*z1*rm[6];
        r[2][4][4] = x2*y4*z4*rm[10] + (y4*z4+6*x2*y2*z4+6*x2*y4*z2)*rm[9] + (6*y2*z4+3*x2*z4+6*y4*z2
                     +36*x2*y2*z2+3*x2*y4)*rm[8]
                     + (3*z4+36*y2*z2+18*x2*z2+3*y4+18*x2*y2)*rm[7] + (18*z2+18*y2+9*x2)*rm[6] + 9*rm[5];
        r[2][5][3] = x2*y5*z3*rm[10] + (y5*z3+10*x2*y3*z3+3*x2*y5*z1)*rm[9] +
                     (10*y3*z3+15*x2*y1*z3+3*y5*z1+30*x2*y3*z1)*rm[8]
                     + (15*y1*z3+30*y3*z1+45*x2*y1*z1)*rm[7] + 45*y1*z1*rm[6];
        r[2][6][2] = x2*y6*z2*rm[10] + (y6*z2+15*x2*y4*z2+x2*y6)*rm[9] + (15*y4*z2+45*x2*y2*z2+y6+15*x2*y4)
                     *rm[8] +
                     (45*y2*z2+15*x2*z2+15*y4+45*x2*y2)*rm[7] + (15*z2+45*y2+15*x2)*rm[6] + 15*rm[5];
        r[2][7][1] = x2*y7*z1*rm[10] + (y7*z1+21*x2*y5*z1)*rm[9] + (21*y5*z1+105*x2*y3*z1)*rm[8] +
                     (105*y3*z1+105*x2*y1*z1)*rm[7] + 105*y1*z1*rm[6];
        r[2][8][0] = x2*y8*rm[10] + (y8+28*x2*y6)*rm[9] + (28*y6+210*x2*y4)*rm[8] +
                     (210*y4+420*x2*y2)*rm[7] +
                     (420*y2+105*x2)*rm[6] + 105*rm[5];
        r[3][0][7] = x3*z7*rm[10] + (3*x1*z7+21*x3*z5)*rm[9] + (63*x1*z5+105*x3*z3)*rm[8] +
                     (315*x1*z3+105*x3*z1)*rm[7] + 315*x1*z1*rm[6];
        r[3][1][6] = x3*y1*z6*rm[10] + (3*x1*y1*z6+15*x3*y1*z4)*rm[9] + (45*x1*y1*z4+45*x3*y1*z2)*rm[8] +
                     (135*x1*y1*z2+15*x3*y1)*rm[7] + 45*x1*y1*rm[6];
        r[3][2][5] = x3*y2*z5*rm[10] + (3*x1*y2*z5+x3*z5+10*x3*y2*z3)*rm[9] +
                     (3*x1*z5+30*x1*y2*z3+10*x3*z3+15*x3*y2*z1)*rm[8]
                     + (30*x1*z3+45*x1*y2*z1+15*x3*z1)*rm[7] + 45*x1*z1*rm[6];
        r[3][3][4] = x3*y3*z4*rm[10] + (3*x1*y3*z4+3*x3*y1*z4+6*x3*y3*z2)*rm[9] +
                     (9*x1*y1*z4+18*x1*y3*z2+18*x3*y1*z2+3*x3*y3)
                     *rm[8] + (54*x1*y1*z2+9*x1*y3+9*x3*y1)*rm[7] + 27*x1*y1*rm[6];
        r[3][4][3] = x3*y4*z3*rm[10] + (3*x1*y4*z3+6*x3*y2*z3+3*x3*y4*z1)*rm[9] +
                     (18*x1*y2*z3+3*x3*z3+9*x1*y4*z1+18*x3*y2*z1)
                     *rm[8] + (9*x1*z3+54*x1*y2*z1+9*x3*z1)*rm[7] + 27*x1*z1*rm[6];
        r[3][5][2] = x3*y5*z2*rm[10] + (3*x1*y5*z2+10*x3*y3*z2+x3*y5)*rm[9] +
                     (30*x1*y3*z2+15*x3*y1*z2+3*x1*y5+10*x3*y3)*rm[8]
                     + (45*x1*y1*z2+30*x1*y3+15*x3*y1)*rm[7] + 45*x1*y1*rm[6];
        r[3][6][1] = x3*y6*z1*rm[10] + (3*x1*y6*z1+15*x3*y4*z1)*rm[9] + (45*x1*y4*z1+45*x3*y2*z1)*rm[8] +
                     (135*x1*y2*z1+15*x3*z1)*rm[7] + 45*x1*z1*rm[6];
        r[3][7][0] = x3*y7*rm[10] + (3*x1*y7+21*x3*y5)*rm[9] + (63*x1*y5+105*x3*y3)*rm[8] +
                     (315*x1*y3+105*x3*y1)*rm[7] + 315*x1*y1*rm[6];
        r[4][0][6] = x4*z6*rm[10] + (6*x2*z6+15*x4*z4)*rm[9] + (3*z6+90*x2*z4+45*x4*z2)*rm[8] +
                     (45*z4+270*x2*z2+15*x4)*rm[7] + (135*z2+90*x2)*rm[6] + 45*rm[5];
        r[4][1][5] = x4*y1*z5*rm[10] + (6*x2*y1*z5+10*x4*y1*z3)*rm[9] + (3*y1*z5+60*x2*y1*z3+15*x4*y1*z1)
                     *rm[8] +
                     (30*y1*z3+90*x2*y1*z1)*rm[7] + 45*y1*z1*rm[6];
        r[4][2][4] = x4*y2*z4*rm[10] + (6*x2*y2*z4+x4*z4+6*x4*y2*z2)*rm[9] + (3*y2*z4+6*x2*z4+36*x2*y2*z2
                     +6*x4*z2+3*x4*y2)*rm[8]
                     + (3*z4+18*y2*z2+36*x2*z2+18*x2*y2+3*x4)*rm[7] + (18*z2+9*y2+18*x2)*rm[6] + 9*rm[5];
        r[4][3][3] = x4*y3*z3*rm[10] + (6*x2*y3*z3+3*x4*y1*z3+3*x4*y3*z1)*rm[9] +
                     (3*y3*z3+18*x2*y1*z3+18*x2*y3*z1+9*x4*y1*z1)
                     *rm[8] + (9*y1*z3+9*y3*z1+54*x2*y1*z1)*rm[7] + 27*y1*z1*rm[6];
        r[4][4][2] = x4*y4*z2*rm[10] + (6*x2*y4*z2+6*x4*y2*z2+x4*y4)*rm[9] + (3*y4*z2+36*x2*y2*z2+3*x4*z2
                     +6*x2*y4+6*x4*y2)*rm[8]
                     + (18*y2*z2+18*x2*z2+3*y4+36*x2*y2+3*x4)*rm[7] + (9*z2+18*y2+18*x2)*rm[6] + 9*rm[5];
        r[4][5][1] = x4*y5*z1*rm[10] + (6*x2*y5*z1+10*x4*y3*z1)*rm[9] + (3*y5*z1+60*x2*y3*z1+15*x4*y1*z1)
                     *rm[8] +
                     (30*y3*z1+90*x2*y1*z1)*rm[7] + 45*y1*z1*rm[6];
        r[4][6][0] = x4*y6*rm[10] + (6*x2*y6+15*x4*y4)*rm[9] + (3*y6+90*x2*y4+45*x4*y2)*rm[8] +
                     (45*y4+270*x2*y2+15*x4)*rm[7] + (135*y2+90*x2)*rm[6] + 45*rm[5];
        r[5][0][5] = x5*z5*rm[10] + (10*x3*z5+10*x5*z3)*rm[9] + (15*x1*z5+100*x3*z3+15*x5*z1)*rm[8] +
                     (150*x1*z3+150*x3*z1)*rm[7] + 225*x1*z1*rm[6];
        r[5][1][4] = x5*y1*z4*rm[10] + (10*x3*y1*z4+6*x5*y1*z2)*rm[9] + (15*x1*y1*z4+60*x3*y1*z2+3*x5*y1)
                     *rm[8] +
                     (90*x1*y1*z2+30*x3*y1)*rm[7] + 45*x1*y1*rm[6];
        r[5][2][3] = x5*y2*z3*rm[10] + (10*x3*y2*z3+x5*z3+3*x5*y2*z1)*rm[9] +
                     (15*x1*y2*z3+10*x3*z3+30*x3*y2*z1+3*x5*z1)*rm[8]
                     + (15*x1*z3+45*x1*y2*z1+30*x3*z1)*rm[7] + 45*x1*z1*rm[6];
        r[5][3][2] = x5*y3*z2*rm[10] + (10*x3*y3*z2+3*x5*y1*z2+x5*y3)*rm[9] +
                     (15*x1*y3*z2+30*x3*y1*z2+10*x3*y3+3*x5*y1)*rm[8]
                     + (45*x1*y1*z2+15*x1*y3+30*x3*y1)*rm[7] + 45*x1*y1*rm[6];
        r[5][4][1] = x5*y4*z1*rm[10] + (10*x3*y4*z1+6*x5*y2*z1)*rm[9] + (15*x1*y4*z1+60*x3*y2*z1+3*x5*z1)
                     *rm[8] +
                     (90*x1*y2*z1+30*x3*z1)*rm[7] + 45*x1*z1*rm[6];
        r[5][5][0] = x5*y5*rm[10] + (10*x3*y5+10*x5*y3)*rm[9] + (15*x1*y5+100*x3*y3+15*x5*y1)*rm[8] +
                     (150*x1*y3+150*x3*y1)*rm[7] + 225*x1*y1*rm[6];
        r[6][0][4] = x6*z4*rm[10] + (15*x4*z4+6*x6*z2)*rm[9] + (45*x2*z4+90*x4*z2+3*x6)*rm[8] +
                     (15*z4+270*x2*z2+45*x4)*rm[7] + (90*z2+135*x2)*rm[6] + 45*rm[5];
        r[6][1][3] = x6*y1*z3*rm[10] + (15*x4*y1*z3+3*x6*y1*z1)*rm[9] + (45*x2*y1*z3+45*x4*y1*z1)*rm[8] +
                     (15*y1*z3+135*x2*y1*z1)*rm[7] + 45*y1*z1*rm[6];
        r[6][2][2] = x6*y2*z2*rm[10] + (15*x4*y2*z2+x6*z2+x6*y2)*rm[9] + (45*x2*y2*z2+15*x4*z2+15*x4*y2+x6)
                     *rm[8] +
                     (15*y2*z2+45*x2*z2+45*x2*y2+15*x4)*rm[7] + (15*z2+15*y2+45*x2)*rm[6] + 15*rm[5];
        r[6][3][1] = x6*y3*z1*rm[10] + (15*x4*y3*z1+3*x6*y1*z1)*rm[9] + (45*x2*y3*z1+45*x4*y1*z1)*rm[8] +
                     (15*y3*z1+135*x2*y1*z1)*rm[7] + 45*y1*z1*rm[6];
        r[6][4][0] = x6*y4*rm[10] + (15*x4*y4+6*x6*y2)*rm[9] + (45*x2*y4+90*x4*y2+3*x6)*rm[8] +
                     (15*y4+270*x2*y2+45*x4)*rm[7] + (90*y2+135*x2)*rm[6] + 45*rm[5];
        r[7][0][3] = x7*z3*rm[10] + (21*x5*z3+3*x7*z1)*rm[9] + (105*x3*z3+63*x5*z1)*rm[8] +
                     (105*x1*z3+315*x3*z1)*rm[7] + 315*x1*z1*rm[6];
        r[7][1][2] = x7*y1*z2*rm[10] + (21*x5*y1*z2+x7*y1)*rm[9] + (105*x3*y1*z2+21*x5*y1)*rm[8] +
                     (105*x1*y1*z2+105*x3*y1)*rm[7] + 105*x1*y1*rm[6];
        r[7][2][1] = x7*y2*z1*rm[10] + (21*x5*y2*z1+x7*z1)*rm[9] + (105*x3*y2*z1+21*x5*z1)*rm[8] +
                     (105*x1*y2*z1+105*x3*z1)*rm[7] + 105*x1*z1*rm[6];
        r[7][3][0] = x7*y3*rm[10] + (21*x5*y3+3*x7*y1)*rm[9] + (105*x3*y3+63*x5*y1)*rm[8] +
                     (105*x1*y3+315*x3*y1)*rm[7] + 315*x1*y1*rm[6];
        r[8][0][2] = x8*z2*rm[10] + (28*x6*z2+x8)*rm[9] + (210*x4*z2+28*x6)*rm[8] +
                     (420*x2*z2+210*x4)*rm[7] +
                     (105*z2+420*x2)*rm[6] + 105*rm[5];
        r[8][1][1] = x8*y1*z1*rm[10] + 28*x6*y1*z1*rm[9] + 210*x4*y1*z1*rm[8] + 420*x2*y1*z1*rm[7] +
                     105*y1*z1*rm[6];
        r[8][2][0] = x8*y2*rm[10] + (28*x6*y2+x8)*rm[9] + (210*x4*y2+28*x6)*rm[8] +
                     (420*x2*y2+210*x4)*rm[7] +
                     (105*y2+420*x2)*rm[6] + 105*rm[5];
        r[9][0][1] = x9*z1*rm[10] + 36*x7*z1*rm[9] + 378*x5*z1*rm[8] + 1260*x3*z1*rm[7] + 945*x1*z1*rm[6];
        r[9][1][0] = x9*y1*rm[10] + 36*x7*y1*rm[9] + 378*x5*y1*rm[8] + 1260*x3*y1*rm[7] + 945*x1*y1*rm[6];
        r[10][0][0] = x10*rm[10] + 45*x8*rm[9] + 630*x6*rm[8] + 3150*x4*rm[7] + 4725*x2*rm[6] + 945*rm[5];
        if ( ltot == 10) return;

        double x11= x10*x1;
        double y11= y10*y1;
        double z11= z10*z1;
        r[0][0][11] = z11*rm[11] + 55*z9*rm[10] + 990*z7*rm[9] + 6930*z5*rm[8] + 17325*z3*rm[7] +
                      10395*z1*rm[6];
        r[0][1][10] = y1*z10*rm[11] + 45*y1*z8*rm[10] + 630*y1*z6*rm[9] + 3150*y1*z4*rm[8] +
                      4725*y1*z2*rm[7] + 945*y1*rm[6];
        r[0][2][9] = y2*z9*rm[11] + (z9+36*y2*z7)*rm[10] + (36*z7+378*y2*z5)*rm[9] +
                     (378*z5+1260*y2*z3)*rm[8] +
                     (1260*z3+945*y2*z1)*rm[7] + 945*z1*rm[6];
        r[0][3][8] = y3*z8*rm[11] + (3*y1*z8+28*y3*z6)*rm[10] + (84*y1*z6+210*y3*z4)*rm[9] +
                     (630*y1*z4+420*y3*z2)*rm[8] +
                     (1260*y1*z2+105*y3)*rm[7] + 315*y1*rm[6];
        r[0][4][7] = y4*z7*rm[11] + (6*y2*z7+21*y4*z5)*rm[10] + (3*z7+126*y2*z5+105*y4*z3)*rm[9] +
                     (63*z5+630*y2*z3+105*y4*z1)*rm[8] + (315*z3+630*y2*z1)*rm[7] + 315*z1*rm[6];
        r[0][5][6] = y5*z6*rm[11] + (10*y3*z6+15*y5*z4)*rm[10] + (15*y1*z6+150*y3*z4+45*y5*z2)*rm[9] +
                     (225*y1*z4+450*y3*z2+15*y5)*rm[8] + (675*y1*z2+150*y3)*rm[7] + 225*y1*rm[6];
        r[0][6][5] = y6*z5*rm[11] + (15*y4*z5+10*y6*z3)*rm[10] + (45*y2*z5+150*y4*z3+15*y6*z1)*rm[9] +
                     (15*z5+450*y2*z3+225*y4*z1)*rm[8] + (150*z3+675*y2*z1)*rm[7] + 225*z1*rm[6];
        r[0][7][4] = y7*z4*rm[11] + (21*y5*z4+6*y7*z2)*rm[10] + (105*y3*z4+126*y5*z2+3*y7)*rm[9] +
                     (105*y1*z4+630*y3*z2+63*y5)*rm[8] + (630*y1*z2+315*y3)*rm[7] + 315*y1*rm[6];
        r[0][8][3] = y8*z3*rm[11] + (28*y6*z3+3*y8*z1)*rm[10] + (210*y4*z3+84*y6*z1)*rm[9] +
                     (420*y2*z3+630*y4*z1)*rm[8] +
                     (105*z3+1260*y2*z1)*rm[7] + 315*z1*rm[6];
        r[0][9][2] = y9*z2*rm[11] + (36*y7*z2+y9)*rm[10] + (378*y5*z2+36*y7)*rm[9] +
                     (1260*y3*z2+378*y5)*rm[8] +
                     (945*y1*z2+1260*y3)*rm[7] + 945*y1*rm[6];
        r[0][10][1] = y10*z1*rm[11] + 45*y8*z1*rm[10] + 630*y6*z1*rm[9] + 3150*y4*z1*rm[8] +
                      4725*y2*z1*rm[7] + 945*z1*rm[6];
        r[0][11][0] = y11*rm[11] + 55*y9*rm[10] + 990*y7*rm[9] + 6930*y5*rm[8] + 17325*y3*rm[7] +
                      10395*y1*rm[6];
        r[1][0][10] = x1*z10*rm[11] + 45*x1*z8*rm[10] + 630*x1*z6*rm[9] + 3150*x1*z4*rm[8] +
                      4725*x1*z2*rm[7] + 945*x1*rm[6];
        r[1][1][9] = x1*y1*z9*rm[11] + 36*x1*y1*z7*rm[10] + 378*x1*y1*z5*rm[9] + 1260*x1*y1*z3*rm[8] +
                     945*x1*y1*z1*rm[7];
        r[1][2][8] = x1*y2*z8*rm[11] + (x1*z8+28*x1*y2*z6)*rm[10] + (28*x1*z6+210*x1*y2*z4)*rm[9] +
                     (210*x1*z4+420*x1*y2*z2)*rm[8] + (420*x1*z2+105*x1*y2)*rm[7] + 105*x1*rm[6];
        r[1][3][7] = x1*y3*z7*rm[11] + (3*x1*y1*z7+21*x1*y3*z5)*rm[10] + (63*x1*y1*z5+105*x1*y3*z3)*rm[9] +
                     (315*x1*y1*z3+105*x1*y3*z1)*rm[8] + 315*x1*y1*z1*rm[7];
        r[1][4][6] = x1*y4*z6*rm[11] + (6*x1*y2*z6+15*x1*y4*z4)*rm[10] + (3*x1*z6+90*x1*y2*z4+45*x1*y4*z2)
                     *rm[9] +
                     (45*x1*z4+270*x1*y2*z2+15*x1*y4)*rm[8] + (135*x1*z2+90*x1*y2)*rm[7] + 45*x1*rm[6];
        r[1][5][5] = x1*y5*z5*rm[11] + (10*x1*y3*z5+10*x1*y5*z3)*rm[10] + (15*x1*y1*z5+100*x1*y3*z3
                     +15*x1*y5*z1)*rm[9] +
                     (150*x1*y1*z3+150*x1*y3*z1)*rm[8] + 225*x1*y1*z1*rm[7];
        r[1][6][4] = x1*y6*z4*rm[11] + (15*x1*y4*z4+6*x1*y6*z2)*rm[10] + (45*x1*y2*z4+90*x1*y4*z2+3*x1*y6)
                     *rm[9] +
                     (15*x1*z4+270*x1*y2*z2+45*x1*y4)*rm[8] + (90*x1*z2+135*x1*y2)*rm[7] + 45*x1*rm[6];
        r[1][7][3] = x1*y7*z3*rm[11] + (21*x1*y5*z3+3*x1*y7*z1)*rm[10] + (105*x1*y3*z3+63*x1*y5*z1)*rm[9] +
                     (105*x1*y1*z3+315*x1*y3*z1)*rm[8] + 315*x1*y1*z1*rm[7];
        r[1][8][2] = x1*y8*z2*rm[11] + (28*x1*y6*z2+x1*y8)*rm[10] + (210*x1*y4*z2+28*x1*y6)*rm[9] +
                     (420*x1*y2*z2+210*x1*y4)*rm[8] + (105*x1*z2+420*x1*y2)*rm[7] + 105*x1*rm[6];
        r[1][9][1] = x1*y9*z1*rm[11] + 36*x1*y7*z1*rm[10] + 378*x1*y5*z1*rm[9] + 1260*x1*y3*z1*rm[8] +
                     945*x1*y1*z1*rm[7];
        r[1][10][0] = x1*y10*rm[11] + 45*x1*y8*rm[10] + 630*x1*y6*rm[9] + 3150*x1*y4*rm[8] +
                      4725*x1*y2*rm[7] + 945*x1*rm[6];
        r[2][0][9] = x2*z9*rm[11] + (z9+36*x2*z7)*rm[10] + (36*z7+378*x2*z5)*rm[9] +
                     (378*z5+1260*x2*z3)*rm[8] +
                     (1260*z3+945*x2*z1)*rm[7] + 945*z1*rm[6];
        r[2][1][8] = x2*y1*z8*rm[11] + (y1*z8+28*x2*y1*z6)*rm[10] + (28*y1*z6+210*x2*y1*z4)*rm[9] +
                     (210*y1*z4+420*x2*y1*z2)*rm[8] + (420*y1*z2+105*x2*y1)*rm[7] + 105*y1*rm[6];
        r[2][2][7] = x2*y2*z7*rm[11] + (y2*z7+x2*z7+21*x2*y2*z5)*rm[10] + (z7+21*y2*z5+21*x2*z5
                     +105*x2*y2*z3)*rm[9] +
                     (21*z5+105*y2*z3+105*x2*z3+105*x2*y2*z1)*rm[8] + (105*z3+105*y2*z1+105*x2*z1)*rm[7] + 105*z1*rm[6];
        r[2][3][6] = x2*y3*z6*rm[11] + (y3*z6+3*x2*y1*z6+15*x2*y3*z4)*rm[10] +
                     (3*y1*z6+15*y3*z4+45*x2*y1*z4+45*x2*y3*z2)*rm[9]
                     + (45*y1*z4+45*y3*z2+135*x2*y1*z2+15*x2*y3)*rm[8] + (135*y1*z2+15*y3+45*x2*y1)*rm[7] + 45*y1*rm[6];
        r[2][4][5] = x2*y4*z5*rm[11] + (y4*z5+6*x2*y2*z5+10*x2*y4*z3)*rm[10] +
                     (6*y2*z5+3*x2*z5+10*y4*z3+60*x2*y2*z3
                      +15*x2*y4*z1)*rm[9] + (3*z5+60*y2*z3+30*x2*z3+15*y4*z1+90*x2*y2*z1)*rm[8] +
                     (30*z3+90*y2*z1+45*x2*z1)*rm[7] +
                     45*z1*rm[6];
        r[2][5][4] = x2*y5*z4*rm[11] + (y5*z4+10*x2*y3*z4+6*x2*y5*z2)*rm[10] +
                     (10*y3*z4+15*x2*y1*z4+6*y5*z2+60*x2*y3*z2
                      +3*x2*y5)*rm[9] + (15*y1*z4+60*y3*z2+90*x2*y1*z2+3*y5+30*x2*y3)*rm[8] +
                     (90*y1*z2+30*y3+45*x2*y1)*rm[7] + 45*y1*rm[6];
        r[2][6][3] = x2*y6*z3*rm[11] + (y6*z3+15*x2*y4*z3+3*x2*y6*z1)*rm[10] +
                     (15*y4*z3+45*x2*y2*z3+3*y6*z1+45*x2*y4*z1)*rm[9]
                     + (45*y2*z3+15*x2*z3+45*y4*z1+135*x2*y2*z1)*rm[8] + (15*z3+135*y2*z1+45*x2*z1)*rm[7] + 45*z1*rm[6];
        r[2][7][2] = x2*y7*z2*rm[11] + (y7*z2+21*x2*y5*z2+x2*y7)*rm[10] + (21*y5*z2+105*x2*y3*z2+y7
                     +21*x2*y5)*rm[9] +
                     (105*y3*z2+105*x2*y1*z2+21*y5+105*x2*y3)*rm[8] + (105*y1*z2+105*y3+105*x2*y1)*rm[7] + 105*y1*rm[6];
        r[2][8][1] = x2*y8*z1*rm[11] + (y8*z1+28*x2*y6*z1)*rm[10] + (28*y6*z1+210*x2*y4*z1)*rm[9] +
                     (210*y4*z1+420*x2*y2*z1)*rm[8] + (420*y2*z1+105*x2*z1)*rm[7] + 105*z1*rm[6];
        r[2][9][0] = x2*y9*rm[11] + (y9+36*x2*y7)*rm[10] + (36*y7+378*x2*y5)*rm[9] +
                     (378*y5+1260*x2*y3)*rm[8] +
                     (1260*y3+945*x2*y1)*rm[7] + 945*y1*rm[6];
        r[3][0][8] = x3*z8*rm[11] + (3*x1*z8+28*x3*z6)*rm[10] + (84*x1*z6+210*x3*z4)*rm[9] +
                     (630*x1*z4+420*x3*z2)*rm[8] +
                     (1260*x1*z2+105*x3)*rm[7] + 315*x1*rm[6];
        r[3][1][7] = x3*y1*z7*rm[11] + (3*x1*y1*z7+21*x3*y1*z5)*rm[10] + (63*x1*y1*z5+105*x3*y1*z3)*rm[9] +
                     (315*x1*y1*z3+105*x3*y1*z1)*rm[8] + 315*x1*y1*z1*rm[7];
        r[3][2][6] = x3*y2*z6*rm[11] + (3*x1*y2*z6+x3*z6+15*x3*y2*z4)*rm[10] +
                     (3*x1*z6+45*x1*y2*z4+15*x3*z4+45*x3*y2*z2)*rm[9]
                     + (45*x1*z4+135*x1*y2*z2+45*x3*z2+15*x3*y2)*rm[8] + (135*x1*z2+45*x1*y2+15*x3)*rm[7] + 45*x1*rm[6];
        r[3][3][5] = x3*y3*z5*rm[11] + (3*x1*y3*z5+3*x3*y1*z5+10*x3*y3*z3)*rm[10] +
                     (9*x1*y1*z5+30*x1*y3*z3+30*x3*y1*z3
                      +15*x3*y3*z1)*rm[9] + (90*x1*y1*z3+45*x1*y3*z1+45*x3*y1*z1)*rm[8] + 135*x1*y1*z1*rm[7];
        r[3][4][4] = x3*y4*z4*rm[11] + (3*x1*y4*z4+6*x3*y2*z4+6*x3*y4*z2)*rm[10] +
                     (18*x1*y2*z4+3*x3*z4+18*x1*y4*z2+36*x3*y2*z2
                      +3*x3*y4)*rm[9] + (9*x1*z4+108*x1*y2*z2+18*x3*z2+9*x1*y4+18*x3*y2)*rm[8] +
                     (54*x1*z2+54*x1*y2+9*x3)*rm[7] + 27*x1*rm[6];
        r[3][5][3] = x3*y5*z3*rm[11] + (3*x1*y5*z3+10*x3*y3*z3+3*x3*y5*z1)*rm[10] +
                     (30*x1*y3*z3+15*x3*y1*z3+9*x1*y5*z1
                      +30*x3*y3*z1)*rm[9] + (45*x1*y1*z3+90*x1*y3*z1+45*x3*y1*z1)*rm[8] + 135*x1*y1*z1*rm[7];
        r[3][6][2] = x3*y6*z2*rm[11] + (3*x1*y6*z2+15*x3*y4*z2+x3*y6)*rm[10] +
                     (45*x1*y4*z2+45*x3*y2*z2+3*x1*y6+15*x3*y4)*rm[9]
                     + (135*x1*y2*z2+15*x3*z2+45*x1*y4+45*x3*y2)*rm[8] + (45*x1*z2+135*x1*y2+15*x3)*rm[7] + 45*x1*rm[6];
        r[3][7][1] = x3*y7*z1*rm[11] + (3*x1*y7*z1+21*x3*y5*z1)*rm[10] + (63*x1*y5*z1+105*x3*y3*z1)*rm[9] +
                     (315*x1*y3*z1+105*x3*y1*z1)*rm[8] + 315*x1*y1*z1*rm[7];
        r[3][8][0] = x3*y8*rm[11] + (3*x1*y8+28*x3*y6)*rm[10] + (84*x1*y6+210*x3*y4)*rm[9] +
                     (630*x1*y4+420*x3*y2)*rm[8] +
                     (1260*x1*y2+105*x3)*rm[7] + 315*x1*rm[6];
        r[4][0][7] = x4*z7*rm[11] + (6*x2*z7+21*x4*z5)*rm[10] + (3*z7+126*x2*z5+105*x4*z3)*rm[9] +
                     (63*z5+630*x2*z3+105*x4*z1)*rm[8] + (315*z3+630*x2*z1)*rm[7] + 315*z1*rm[6];
        r[4][1][6] = x4*y1*z6*rm[11] + (6*x2*y1*z6+15*x4*y1*z4)*rm[10] + (3*y1*z6+90*x2*y1*z4+45*x4*y1*z2)
                     *rm[9] +
                     (45*y1*z4+270*x2*y1*z2+15*x4*y1)*rm[8] + (135*y1*z2+90*x2*y1)*rm[7] + 45*y1*rm[6];
        r[4][2][5] = x4*y2*z5*rm[11] + (6*x2*y2*z5+x4*z5+10*x4*y2*z3)*rm[10] +
                     (3*y2*z5+6*x2*z5+60*x2*y2*z3+10*x4*z3
                      +15*x4*y2*z1)*rm[9] + (3*z5+30*y2*z3+60*x2*z3+90*x2*y2*z1+15*x4*z1)*rm[8] +
                     (30*z3+45*y2*z1+90*x2*z1)*rm[7] +
                     45*z1*rm[6];
        r[4][3][4] = x4*y3*z4*rm[11] + (6*x2*y3*z4+3*x4*y1*z4+6*x4*y3*z2)*rm[10] +
                     (3*y3*z4+18*x2*y1*z4+36*x2*y3*z2+18*x4*y1*z2
                      +3*x4*y3)*rm[9] + (9*y1*z4+18*y3*z2+108*x2*y1*z2+18*x2*y3+9*x4*y1)*rm[8] +
                     (54*y1*z2+9*y3+54*x2*y1)*rm[7] + 27*y1*rm[6];
        r[4][4][3] = x4*y4*z3*rm[11] + (6*x2*y4*z3+6*x4*y2*z3+3*x4*y4*z1)*rm[10] +
                     (3*y4*z3+36*x2*y2*z3+3*x4*z3+18*x2*y4*z1
                      +18*x4*y2*z1)*rm[9] + (18*y2*z3+18*x2*z3+9*y4*z1+108*x2*y2*z1+9*x4*z1)*rm[8] +
                     (9*z3+54*y2*z1+54*x2*z1)*rm[7] +
                     27*z1*rm[6];
        r[4][5][2] = x4*y5*z2*rm[11] + (6*x2*y5*z2+10*x4*y3*z2+x4*y5)*rm[10] +
                     (3*y5*z2+60*x2*y3*z2+15*x4*y1*z2+6*x2*y5
                      +10*x4*y3)*rm[9] + (30*y3*z2+90*x2*y1*z2+3*y5+60*x2*y3+15*x4*y1)*rm[8] +
                     (45*y1*z2+30*y3+90*x2*y1)*rm[7] + 45*y1*rm[6];
        r[4][6][1] = x4*y6*z1*rm[11] + (6*x2*y6*z1+15*x4*y4*z1)*rm[10] + (3*y6*z1+90*x2*y4*z1+45*x4*y2*z1)
                     *rm[9] +
                     (45*y4*z1+270*x2*y2*z1+15*x4*z1)*rm[8] + (135*y2*z1+90*x2*z1)*rm[7] + 45*z1*rm[6];
        r[4][7][0] = x4*y7*rm[11] + (6*x2*y7+21*x4*y5)*rm[10] + (3*y7+126*x2*y5+105*x4*y3)*rm[9] +
                     (63*y5+630*x2*y3+105*x4*y1)*rm[8] + (315*y3+630*x2*y1)*rm[7] + 315*y1*rm[6];
        r[5][0][6] = x5*z6*rm[11] + (10*x3*z6+15*x5*z4)*rm[10] + (15*x1*z6+150*x3*z4+45*x5*z2)*rm[9] +
                     (225*x1*z4+450*x3*z2+15*x5)*rm[8] + (675*x1*z2+150*x3)*rm[7] + 225*x1*rm[6];
        r[5][1][5] = x5*y1*z5*rm[11] + (10*x3*y1*z5+10*x5*y1*z3)*rm[10] + (15*x1*y1*z5+100*x3*y1*z3
                     +15*x5*y1*z1)*rm[9] +
                     (150*x1*y1*z3+150*x3*y1*z1)*rm[8] + 225*x1*y1*z1*rm[7];
        r[5][2][4] = x5*y2*z4*rm[11] + (10*x3*y2*z4+x5*z4+6*x5*y2*z2)*rm[10] +
                     (15*x1*y2*z4+10*x3*z4+60*x3*y2*z2+6*x5*z2
                      +3*x5*y2)*rm[9] + (15*x1*z4+90*x1*y2*z2+60*x3*z2+30*x3*y2+3*x5)*rm[8] +
                     (90*x1*z2+45*x1*y2+30*x3)*rm[7] + 45*x1*rm[6];
        r[5][3][3] = x5*y3*z3*rm[11] + (10*x3*y3*z3+3*x5*y1*z3+3*x5*y3*z1)*rm[10] +
                     (15*x1*y3*z3+30*x3*y1*z3+30*x3*y3*z1
                      +9*x5*y1*z1)*rm[9] + (45*x1*y1*z3+45*x1*y3*z1+90*x3*y1*z1)*rm[8] + 135*x1*y1*z1*rm[7];
        r[5][4][2] = x5*y4*z2*rm[11] + (10*x3*y4*z2+6*x5*y2*z2+x5*y4)*rm[10] +
                     (15*x1*y4*z2+60*x3*y2*z2+3*x5*z2+10*x3*y4
                      +6*x5*y2)*rm[9] + (90*x1*y2*z2+30*x3*z2+15*x1*y4+60*x3*y2+3*x5)*rm[8] +
                     (45*x1*z2+90*x1*y2+30*x3)*rm[7] + 45*x1*rm[6];
        r[5][5][1] = x5*y5*z1*rm[11] + (10*x3*y5*z1+10*x5*y3*z1)*rm[10] + (15*x1*y5*z1+100*x3*y3*z1
                     +15*x5*y1*z1)*rm[9] +
                     (150*x1*y3*z1+150*x3*y1*z1)*rm[8] + 225*x1*y1*z1*rm[7];
        r[5][6][0] = x5*y6*rm[11] + (10*x3*y6+15*x5*y4)*rm[10] + (15*x1*y6+150*x3*y4+45*x5*y2)*rm[9] +
                     (225*x1*y4+450*x3*y2+15*x5)*rm[8] + (675*x1*y2+150*x3)*rm[7] + 225*x1*rm[6];
        r[6][0][5] = x6*z5*rm[11] + (15*x4*z5+10*x6*z3)*rm[10] + (45*x2*z5+150*x4*z3+15*x6*z1)*rm[9] +
                     (15*z5+450*x2*z3+225*x4*z1)*rm[8] + (150*z3+675*x2*z1)*rm[7] + 225*z1*rm[6];
        r[6][1][4] = x6*y1*z4*rm[11] + (15*x4*y1*z4+6*x6*y1*z2)*rm[10] + (45*x2*y1*z4+90*x4*y1*z2+3*x6*y1)
                     *rm[9] +
                     (15*y1*z4+270*x2*y1*z2+45*x4*y1)*rm[8] + (90*y1*z2+135*x2*y1)*rm[7] + 45*y1*rm[6];
        r[6][2][3] = x6*y2*z3*rm[11] + (15*x4*y2*z3+x6*z3+3*x6*y2*z1)*rm[10] +
                     (45*x2*y2*z3+15*x4*z3+45*x4*y2*z1+3*x6*z1)*rm[9]
                     + (15*y2*z3+45*x2*z3+135*x2*y2*z1+45*x4*z1)*rm[8] + (15*z3+45*y2*z1+135*x2*z1)*rm[7] + 45*z1*rm[6];
        r[6][3][2] = x6*y3*z2*rm[11] + (15*x4*y3*z2+3*x6*y1*z2+x6*y3)*rm[10] +
                     (45*x2*y3*z2+45*x4*y1*z2+15*x4*y3+3*x6*y1)*rm[9]
                     + (15*y3*z2+135*x2*y1*z2+45*x2*y3+45*x4*y1)*rm[8] + (45*y1*z2+15*y3+135*x2*y1)*rm[7] + 45*y1*rm[6];
        r[6][4][1] = x6*y4*z1*rm[11] + (15*x4*y4*z1+6*x6*y2*z1)*rm[10] + (45*x2*y4*z1+90*x4*y2*z1+3*x6*z1)
                     *rm[9] +
                     (15*y4*z1+270*x2*y2*z1+45*x4*z1)*rm[8] + (90*y2*z1+135*x2*z1)*rm[7] + 45*z1*rm[6];
        r[6][5][0] = x6*y5*rm[11] + (15*x4*y5+10*x6*y3)*rm[10] + (45*x2*y5+150*x4*y3+15*x6*y1)*rm[9] +
                     (15*y5+450*x2*y3+225*x4*y1)*rm[8] + (150*y3+675*x2*y1)*rm[7] + 225*y1*rm[6];
        r[7][0][4] = x7*z4*rm[11] + (21*x5*z4+6*x7*z2)*rm[10] + (105*x3*z4+126*x5*z2+3*x7)*rm[9] +
                     (105*x1*z4+630*x3*z2+63*x5)*rm[8] + (630*x1*z2+315*x3)*rm[7] + 315*x1*rm[6];
        r[7][1][3] = x7*y1*z3*rm[11] + (21*x5*y1*z3+3*x7*y1*z1)*rm[10] + (105*x3*y1*z3+63*x5*y1*z1)*rm[9] +
                     (105*x1*y1*z3+315*x3*y1*z1)*rm[8] + 315*x1*y1*z1*rm[7];
        r[7][2][2] = x7*y2*z2*rm[11] + (21*x5*y2*z2+x7*z2+x7*y2)*rm[10] + (105*x3*y2*z2+21*x5*z2+21*x5*y2
                     +x7)*rm[9] +
                     (105*x1*y2*z2+105*x3*z2+105*x3*y2+21*x5)*rm[8] + (105*x1*z2+105*x1*y2+105*x3)*rm[7] + 105*x1*rm[6];
        r[7][3][1] = x7*y3*z1*rm[11] + (21*x5*y3*z1+3*x7*y1*z1)*rm[10] + (105*x3*y3*z1+63*x5*y1*z1)*rm[9] +
                     (105*x1*y3*z1+315*x3*y1*z1)*rm[8] + 315*x1*y1*z1*rm[7];
        r[7][4][0] = x7*y4*rm[11] + (21*x5*y4+6*x7*y2)*rm[10] + (105*x3*y4+126*x5*y2+3*x7)*rm[9] +
                     (105*x1*y4+630*x3*y2+63*x5)*rm[8] + (630*x1*y2+315*x3)*rm[7] + 315*x1*rm[6];
        r[8][0][3] = x8*z3*rm[11] + (28*x6*z3+3*x8*z1)*rm[10] + (210*x4*z3+84*x6*z1)*rm[9] +
                     (420*x2*z3+630*x4*z1)*rm[8] +
                     (105*z3+1260*x2*z1)*rm[7] + 315*z1*rm[6];
        r[8][1][2] = x8*y1*z2*rm[11] + (28*x6*y1*z2+x8*y1)*rm[10] + (210*x4*y1*z2+28*x6*y1)*rm[9] +
                     (420*x2*y1*z2+210*x4*y1)*rm[8] + (105*y1*z2+420*x2*y1)*rm[7] + 105*y1*rm[6];
        r[8][2][1] = x8*y2*z1*rm[11] + (28*x6*y2*z1+x8*z1)*rm[10] + (210*x4*y2*z1+28*x6*z1)*rm[9] +
                     (420*x2*y2*z1+210*x4*z1)*rm[8] + (105*y2*z1+420*x2*z1)*rm[7] + 105*z1*rm[6];
        r[8][3][0] = x8*y3*rm[11] + (28*x6*y3+3*x8*y1)*rm[10] + (210*x4*y3+84*x6*y1)*rm[9] +
                     (420*x2*y3+630*x4*y1)*rm[8] +
                     (105*y3+1260*x2*y1)*rm[7] + 315*y1*rm[6];
        r[9][0][2] = x9*z2*rm[11] + (36*x7*z2+x9)*rm[10] + (378*x5*z2+36*x7)*rm[9] +
                     (1260*x3*z2+378*x5)*rm[8] +
                     (945*x1*z2+1260*x3)*rm[7] + 945*x1*rm[6];
        r[9][1][1] = x9*y1*z1*rm[11] + 36*x7*y1*z1*rm[10] + 378*x5*y1*z1*rm[9] + 1260*x3*y1*z1*rm[8] +
                     945*x1*y1*z1*rm[7];
        r[9][2][0] = x9*y2*rm[11] + (36*x7*y2+x9)*rm[10] + (378*x5*y2+36*x7)*rm[9] +
                     (1260*x3*y2+378*x5)*rm[8] +
                     (945*x1*y2+1260*x3)*rm[7] + 945*x1*rm[6];
        r[10][0][1] = x10*z1*rm[11] + 45*x8*z1*rm[10] + 630*x6*z1*rm[9] + 3150*x4*z1*rm[8] +
                      4725*x2*z1*rm[7] + 945*z1*rm[6];
        r[10][1][0] = x10*y1*rm[11] + 45*x8*y1*rm[10] + 630*x6*y1*rm[9] + 3150*x4*y1*rm[8] +
                      4725*x2*y1*rm[7] + 945*y1*rm[6];
        r[11][0][0] = x11*rm[11] + 55*x9*rm[10] + 990*x7*rm[9] + 6930*x5*rm[8] + 17325*x3*rm[7] +
                      10395*x1*rm[6];
        if ( ltot == 11) return;

        double x12= x11*x1;
        double y12= y11*y1;
        double z12= z11*z1;
        r[0][0][12] = z12*rm[12] + 66*z10*rm[11] + 1485*z8*rm[10] + 13860*z6*rm[9] + 51975*z4*rm[8] +
                      62370*z2*rm[7] +
                      10395*rm[6];
        r[0][1][11] = y1*z11*rm[12] + 55*y1*z9*rm[11] + 990*y1*z7*rm[10] + 6930*y1*z5*rm[9] +
                      17325*y1*z3*rm[8] +
                      10395*y1*z1*rm[7];
        r[0][2][10] = y2*z10*rm[12] + (z10+45*y2*z8)*rm[11] + (45*z8+630*y2*z6)*rm[10] +
                      (630*z6+3150*y2*z4)*rm[9] +
                      (3150*z4+4725*y2*z2)*rm[8] + (4725*z2+945*y2)*rm[7] + 945*rm[6];
        r[0][3][9] = y3*z9*rm[12] + (3*y1*z9+36*y3*z7)*rm[11] + (108*y1*z7+378*y3*z5)*rm[10] +
                     (1134*y1*z5+1260*y3*z3)*rm[9] +
                     (3780*y1*z3+945*y3*z1)*rm[8] + 2835*y1*z1*rm[7];
        r[0][4][8] = y4*z8*rm[12] + (6*y2*z8+28*y4*z6)*rm[11] + (3*z8+168*y2*z6+210*y4*z4)*rm[10] +
                     (84*z6+1260*y2*z4+420*y4*z2)*rm[9] + (630*z4+2520*y2*z2+105*y4)*rm[8] +
                     (1260*z2+630*y2)*rm[7] + 315*rm[6];
        r[0][5][7] = y5*z7*rm[12] + (10*y3*z7+21*y5*z5)*rm[11] + (15*y1*z7+210*y3*z5+105*y5*z3)*rm[10] +
                     (315*y1*z5+1050*y3*z3+105*y5*z1)*rm[9] + (1575*y1*z3+1050*y3*z1)*rm[8] + 1575*y1*z1*rm[7];
        r[0][6][6] = y6*z6*rm[12] + (15*y4*z6+15*y6*z4)*rm[11] + (45*y2*z6+225*y4*z4+45*y6*z2)*rm[10] +
                     (15*z6+675*y2*z4+675*y4*z2+15*y6)*rm[9] + (225*z4+2025*y2*z2+225*y4)*rm[8] +
                     (675*z2+675*y2)*rm[7] + 225*rm[6];
        r[0][7][5] = y7*z5*rm[12] + (21*y5*z5+10*y7*z3)*rm[11] + (105*y3*z5+210*y5*z3+15*y7*z1)*rm[10] +
                     (105*y1*z5+1050*y3*z3+315*y5*z1)*rm[9] + (1050*y1*z3+1575*y3*z1)*rm[8] + 1575*y1*z1*rm[7];
        r[0][8][4] = y8*z4*rm[12] + (28*y6*z4+6*y8*z2)*rm[11] + (210*y4*z4+168*y6*z2+3*y8)*rm[10] +
                     (420*y2*z4+1260*y4*z2+84*y6)*rm[9] + (105*z4+2520*y2*z2+630*y4)*rm[8] +
                     (630*z2+1260*y2)*rm[7] + 315*rm[6];
        r[0][9][3] = y9*z3*rm[12] + (36*y7*z3+3*y9*z1)*rm[11] + (378*y5*z3+108*y7*z1)*rm[10] +
                     (1260*y3*z3+1134*y5*z1)*rm[9] +
                     (945*y1*z3+3780*y3*z1)*rm[8] + 2835*y1*z1*rm[7];
        r[0][10][2] = y10*z2*rm[12] + (45*y8*z2+y10)*rm[11] + (630*y6*z2+45*y8)*rm[10] +
                      (3150*y4*z2+630*y6)*rm[9] +
                      (4725*y2*z2+3150*y4)*rm[8] + (945*z2+4725*y2)*rm[7] + 945*rm[6];
        r[0][11][1] = y11*z1*rm[12] + 55*y9*z1*rm[11] + 990*y7*z1*rm[10] + 6930*y5*z1*rm[9] +
                      17325*y3*z1*rm[8] +
                      10395*y1*z1*rm[7];
        r[0][12][0] = y12*rm[12] + 66*y10*rm[11] + 1485*y8*rm[10] + 13860*y6*rm[9] + 51975*y4*rm[8] +
                      62370*y2*rm[7] +
                      10395*rm[6];
        r[1][0][11] = x1*z11*rm[12] + 55*x1*z9*rm[11] + 990*x1*z7*rm[10] + 6930*x1*z5*rm[9] +
                      17325*x1*z3*rm[8] +
                      10395*x1*z1*rm[7];
        r[1][1][10] = x1*y1*z10*rm[12] + 45*x1*y1*z8*rm[11] + 630*x1*y1*z6*rm[10] + 3150*x1*y1*z4*rm[9] +
                      4725*x1*y1*z2*rm[8] +
                      945*x1*y1*rm[7];
        r[1][2][9] = x1*y2*z9*rm[12] + (x1*z9+36*x1*y2*z7)*rm[11] + (36*x1*z7+378*x1*y2*z5)*rm[10] +
                     (378*x1*z5+1260*x1*y2*z3)*rm[9] + (1260*x1*z3+945*x1*y2*z1)*rm[8] + 945*x1*z1*rm[7];
        r[1][3][8] = x1*y3*z8*rm[12] + (3*x1*y1*z8+28*x1*y3*z6)*rm[11] + (84*x1*y1*z6+210*x1*y3*z4)*rm[10] +
                     (630*x1*y1*z4+420*x1*y3*z2)*rm[9] + (1260*x1*y1*z2+105*x1*y3)*rm[8] + 315*x1*y1*rm[7];
        r[1][4][7] = x1*y4*z7*rm[12] + (6*x1*y2*z7+21*x1*y4*z5)*rm[11] + (3*x1*z7+126*x1*y2*z5+105*x1*y4*z3)
                     *rm[10] +
                     (63*x1*z5+630*x1*y2*z3+105*x1*y4*z1)*rm[9] + (315*x1*z3+630*x1*y2*z1)*rm[8] + 315*x1*z1*rm[7];
        r[1][5][6] = x1*y5*z6*rm[12] + (10*x1*y3*z6+15*x1*y5*z4)*rm[11] + (15*x1*y1*z6+150*x1*y3*z4
                     +45*x1*y5*z2)*rm[10] +
                     (225*x1*y1*z4+450*x1*y3*z2+15*x1*y5)*rm[9] + (675*x1*y1*z2+150*x1*y3)*rm[8] + 225*x1*y1*rm[7];
        r[1][6][5] = x1*y6*z5*rm[12] + (15*x1*y4*z5+10*x1*y6*z3)*rm[11] + (45*x1*y2*z5+150*x1*y4*z3
                     +15*x1*y6*z1)*rm[10] +
                     (15*x1*z5+450*x1*y2*z3+225*x1*y4*z1)*rm[9] + (150*x1*z3+675*x1*y2*z1)*rm[8] + 225*x1*z1*rm[7];
        r[1][7][4] = x1*y7*z4*rm[12] + (21*x1*y5*z4+6*x1*y7*z2)*rm[11] + (105*x1*y3*z4+126*x1*y5*z2+3*x1*y7)
                     *rm[10] +
                     (105*x1*y1*z4+630*x1*y3*z2+63*x1*y5)*rm[9] + (630*x1*y1*z2+315*x1*y3)*rm[8] + 315*x1*y1*rm[7];
        r[1][8][3] = x1*y8*z3*rm[12] + (28*x1*y6*z3+3*x1*y8*z1)*rm[11] + (210*x1*y4*z3+84*x1*y6*z1)*rm[10] +
                     (420*x1*y2*z3+630*x1*y4*z1)*rm[9] + (105*x1*z3+1260*x1*y2*z1)*rm[8] + 315*x1*z1*rm[7];
        r[1][9][2] = x1*y9*z2*rm[12] + (36*x1*y7*z2+x1*y9)*rm[11] + (378*x1*y5*z2+36*x1*y7)*rm[10] +
                     (1260*x1*y3*z2+378*x1*y5)*rm[9] + (945*x1*y1*z2+1260*x1*y3)*rm[8] + 945*x1*y1*rm[7];
        r[1][10][1] = x1*y10*z1*rm[12] + 45*x1*y8*z1*rm[11] + 630*x1*y6*z1*rm[10] + 3150*x1*y4*z1*rm[9] +
                      4725*x1*y2*z1*rm[8] +
                      945*x1*z1*rm[7];
        r[1][11][0] = x1*y11*rm[12] + 55*x1*y9*rm[11] + 990*x1*y7*rm[10] + 6930*x1*y5*rm[9] +
                      17325*x1*y3*rm[8] +
                      10395*x1*y1*rm[7];
        r[2][0][10] = x2*z10*rm[12] + (z10+45*x2*z8)*rm[11] + (45*z8+630*x2*z6)*rm[10] +
                      (630*z6+3150*x2*z4)*rm[9] +
                      (3150*z4+4725*x2*z2)*rm[8] + (4725*z2+945*x2)*rm[7] + 945*rm[6];
        r[2][1][9] = x2*y1*z9*rm[12] + (y1*z9+36*x2*y1*z7)*rm[11] + (36*y1*z7+378*x2*y1*z5)*rm[10] +
                     (378*y1*z5+1260*x2*y1*z3)*rm[9] + (1260*y1*z3+945*x2*y1*z1)*rm[8] + 945*y1*z1*rm[7];
        r[2][2][8] = x2*y2*z8*rm[12] + (y2*z8+x2*z8+28*x2*y2*z6)*rm[11] + (z8+28*y2*z6+28*x2*z6
                     +210*x2*y2*z4)*rm[10] +
                     (28*z6+210*y2*z4+210*x2*z4+420*x2*y2*z2)*rm[9] + (210*z4+420*y2*z2+420*x2*z2+105*x2*y2)*rm[8] +
                     (420*z2+105*y2+105*x2)*rm[7] + 105*rm[6];
        r[2][3][7] = x2*y3*z7*rm[12] + (y3*z7+3*x2*y1*z7+21*x2*y3*z5)*rm[11] +
                     (3*y1*z7+21*y3*z5+63*x2*y1*z5+105*x2*y3*z3)
                     *rm[10] + (63*y1*z5+105*y3*z3+315*x2*y1*z3+105*x2*y3*z1)*rm[9] + (315*y1*z3+105*y3*z1+315*x2*y1*z1)
                     *rm[8] +
                     315*y1*z1*rm[7];
        r[2][4][6] = x2*y4*z6*rm[12] + (y4*z6+6*x2*y2*z6+15*x2*y4*z4)*rm[11] +
                     (6*y2*z6+3*x2*z6+15*y4*z4+90*x2*y2*z4
                      +45*x2*y4*z2)*rm[10] + (3*z6+90*y2*z4+45*x2*z4+45*y4*z2+270*x2*y2*z2+15*x2*y4)*rm[9] +
                     (45*z4+270*y2*z2+135*x2*z2+15*y4+90*x2*y2)*rm[8] + (135*z2+90*y2+45*x2)*rm[7] + 45*rm[6];
        r[2][5][5] = x2*y5*z5*rm[12] + (y5*z5+10*x2*y3*z5+10*x2*y5*z3)*rm[11] +
                     (10*y3*z5+15*x2*y1*z5+10*y5*z3+100*x2*y3*z3
                      +15*x2*y5*z1)*rm[10] + (15*y1*z5+100*y3*z3+150*x2*y1*z3+15*y5*z1+150*x2*y3*z1)*rm[9] +
                     (150*y1*z3+150*y3*z1+225*x2*y1*z1)*rm[8] + 225*y1*z1*rm[7];
        r[2][6][4] = x2*y6*z4*rm[12] + (y6*z4+15*x2*y4*z4+6*x2*y6*z2)*rm[11] +
                     (15*y4*z4+45*x2*y2*z4+6*y6*z2+90*x2*y4*z2
                      +3*x2*y6)*rm[10] + (45*y2*z4+15*x2*z4+90*y4*z2+270*x2*y2*z2+3*y6+45*x2*y4)*rm[9] +
                     (15*z4+270*y2*z2+90*x2*z2+45*y4
                      +135*x2*y2)*rm[8] + (90*z2+135*y2+45*x2)*rm[7] + 45*rm[6];
        r[2][7][3] = x2*y7*z3*rm[12] + (y7*z3+21*x2*y5*z3+3*x2*y7*z1)*rm[11] +
                     (21*y5*z3+105*x2*y3*z3+3*y7*z1+63*x2*y5*z1)
                     *rm[10] + (105*y3*z3+105*x2*y1*z3+63*y5*z1+315*x2*y3*z1)*rm[9] + (105*y1*z3+315*y3*z1+315*x2*y1*z1)
                     *rm[8] +
                     315*y1*z1*rm[7];
        r[2][8][2] = x2*y8*z2*rm[12] + (y8*z2+28*x2*y6*z2+x2*y8)*rm[11] + (28*y6*z2+210*x2*y4*z2+y8
                     +28*x2*y6)*rm[10] +
                     (210*y4*z2+420*x2*y2*z2+28*y6+210*x2*y4)*rm[9] + (420*y2*z2+105*x2*z2+210*y4+420*x2*y2)*rm[8] +
                     (105*z2+420*y2+105*x2)*rm[7] + 105*rm[6];
        r[2][9][1] = x2*y9*z1*rm[12] + (y9*z1+36*x2*y7*z1)*rm[11] + (36*y7*z1+378*x2*y5*z1)*rm[10] +
                     (378*y5*z1+1260*x2*y3*z1)*rm[9] + (1260*y3*z1+945*x2*y1*z1)*rm[8] + 945*y1*z1*rm[7];
        r[2][10][0] = x2*y10*rm[12] + (y10+45*x2*y8)*rm[11] + (45*y8+630*x2*y6)*rm[10] +
                      (630*y6+3150*x2*y4)*rm[9] +
                      (3150*y4+4725*x2*y2)*rm[8] + (4725*y2+945*x2)*rm[7] + 945*rm[6];
        r[3][0][9] = x3*z9*rm[12] + (3*x1*z9+36*x3*z7)*rm[11] + (108*x1*z7+378*x3*z5)*rm[10] +
                     (1134*x1*z5+1260*x3*z3)*rm[9] +
                     (3780*x1*z3+945*x3*z1)*rm[8] + 2835*x1*z1*rm[7];
        r[3][1][8] = x3*y1*z8*rm[12] + (3*x1*y1*z8+28*x3*y1*z6)*rm[11] + (84*x1*y1*z6+210*x3*y1*z4)*rm[10] +
                     (630*x1*y1*z4+420*x3*y1*z2)*rm[9] + (1260*x1*y1*z2+105*x3*y1)*rm[8] + 315*x1*y1*rm[7];
        r[3][2][7] = x3*y2*z7*rm[12] + (3*x1*y2*z7+x3*z7+21*x3*y2*z5)*rm[11] +
                     (3*x1*z7+63*x1*y2*z5+21*x3*z5+105*x3*y2*z3)
                     *rm[10] + (63*x1*z5+315*x1*y2*z3+105*x3*z3+105*x3*y2*z1)*rm[9] + (315*x1*z3+315*x1*y2*z1+105*x3*z1)
                     *rm[8] +
                     315*x1*z1*rm[7];
        r[3][3][6] = x3*y3*z6*rm[12] + (3*x1*y3*z6+3*x3*y1*z6+15*x3*y3*z4)*rm[11] +
                     (9*x1*y1*z6+45*x1*y3*z4+45*x3*y1*z4
                      +45*x3*y3*z2)*rm[10] + (135*x1*y1*z4+135*x1*y3*z2+135*x3*y1*z2+15*x3*y3)*rm[9] +
                     (405*x1*y1*z2+45*x1*y3+45*x3*y1)*rm[8]
                     + 135*x1*y1*rm[7];
        r[3][4][5] = x3*y4*z5*rm[12] + (3*x1*y4*z5+6*x3*y2*z5+10*x3*y4*z3)*rm[11] +
                     (18*x1*y2*z5+3*x3*z5+30*x1*y4*z3+60*x3*y2*z3
                      +15*x3*y4*z1)*rm[10] + (9*x1*z5+180*x1*y2*z3+30*x3*z3+45*x1*y4*z1+90*x3*y2*z1)*rm[9] +
                     (90*x1*z3+270*x1*y2*z1+45*x3*z1)*rm[8] + 135*x1*z1*rm[7];
        r[3][5][4] = x3*y5*z4*rm[12] + (3*x1*y5*z4+10*x3*y3*z4+6*x3*y5*z2)*rm[11] +
                     (30*x1*y3*z4+15*x3*y1*z4+18*x1*y5*z2
                      +60*x3*y3*z2+3*x3*y5)*rm[10] + (45*x1*y1*z4+180*x1*y3*z2+90*x3*y1*z2+9*x1*y5+30*x3*y3)*rm[9] +
                     (270*x1*y1*z2+90*x1*y3+45*x3*y1)*rm[8] + 135*x1*y1*rm[7];
        r[3][6][3] = x3*y6*z3*rm[12] + (3*x1*y6*z3+15*x3*y4*z3+3*x3*y6*z1)*rm[11] +
                     (45*x1*y4*z3+45*x3*y2*z3+9*x1*y6*z1
                      +45*x3*y4*z1)*rm[10] + (135*x1*y2*z3+15*x3*z3+135*x1*y4*z1+135*x3*y2*z1)*rm[9] +
                     (45*x1*z3+405*x1*y2*z1+45*x3*z1)*rm[8]
                     + 135*x1*z1*rm[7];
        r[3][7][2] = x3*y7*z2*rm[12] + (3*x1*y7*z2+21*x3*y5*z2+x3*y7)*rm[11] +
                     (63*x1*y5*z2+105*x3*y3*z2+3*x1*y7+21*x3*y5)
                     *rm[10] + (315*x1*y3*z2+105*x3*y1*z2+63*x1*y5+105*x3*y3)*rm[9] + (315*x1*y1*z2+315*x1*y3+105*x3*y1)
                     *rm[8] +
                     315*x1*y1*rm[7];
        r[3][8][1] = x3*y8*z1*rm[12] + (3*x1*y8*z1+28*x3*y6*z1)*rm[11] + (84*x1*y6*z1+210*x3*y4*z1)*rm[10] +
                     (630*x1*y4*z1+420*x3*y2*z1)*rm[9] + (1260*x1*y2*z1+105*x3*z1)*rm[8] + 315*x1*z1*rm[7];
        r[3][9][0] = x3*y9*rm[12] + (3*x1*y9+36*x3*y7)*rm[11] + (108*x1*y7+378*x3*y5)*rm[10] +
                     (1134*x1*y5+1260*x3*y3)*rm[9] +
                     (3780*x1*y3+945*x3*y1)*rm[8] + 2835*x1*y1*rm[7];
        r[4][0][8] = x4*z8*rm[12] + (6*x2*z8+28*x4*z6)*rm[11] + (3*z8+168*x2*z6+210*x4*z4)*rm[10] +
                     (84*z6+1260*x2*z4+420*x4*z2)*rm[9] + (630*z4+2520*x2*z2+105*x4)*rm[8] +
                     (1260*z2+630*x2)*rm[7] + 315*rm[6];
        r[4][1][7] = x4*y1*z7*rm[12] + (6*x2*y1*z7+21*x4*y1*z5)*rm[11] + (3*y1*z7+126*x2*y1*z5+105*x4*y1*z3)
                     *rm[10] +
                     (63*y1*z5+630*x2*y1*z3+105*x4*y1*z1)*rm[9] + (315*y1*z3+630*x2*y1*z1)*rm[8] + 315*y1*z1*rm[7];
        r[4][2][6] = x4*y2*z6*rm[12] + (6*x2*y2*z6+x4*z6+15*x4*y2*z4)*rm[11] +
                     (3*y2*z6+6*x2*z6+90*x2*y2*z4+15*x4*z4
                      +45*x4*y2*z2)*rm[10] + (3*z6+45*y2*z4+90*x2*z4+270*x2*y2*z2+45*x4*z2+15*x4*y2)*rm[9] +
                     (45*z4+135*y2*z2+270*x2*z2+90*x2*y2+15*x4)*rm[8] + (135*z2+45*y2+90*x2)*rm[7] + 45*rm[6];
        r[4][3][5] = x4*y3*z5*rm[12] + (6*x2*y3*z5+3*x4*y1*z5+10*x4*y3*z3)*rm[11] +
                     (3*y3*z5+18*x2*y1*z5+60*x2*y3*z3+30*x4*y1*z3
                      +15*x4*y3*z1)*rm[10] + (9*y1*z5+30*y3*z3+180*x2*y1*z3+90*x2*y3*z1+45*x4*y1*z1)*rm[9] +
                     (90*y1*z3+45*y3*z1+270*x2*y1*z1)*rm[8] + 135*y1*z1*rm[7];
        r[4][4][4] = x4*y4*z4*rm[12] + (6*x2*y4*z4+6*x4*y2*z4+6*x4*y4*z2)*rm[11] +
                     (3*y4*z4+36*x2*y2*z4+3*x4*z4+36*x2*y4*z2
                      +36*x4*y2*z2+3*x4*y4)*rm[10] + (18*y2*z4+18*x2*z4+18*y4*z2+216*x2*y2*z2+18*x4*z2+18*x2*y4+18*x4*y2)
                     *rm[9] +
                     (9*z4+108*y2*z2+108*x2*z2+9*y4+108*x2*y2+9*x4)*rm[8] + (54*z2+54*y2+54*x2)*rm[7] + 27*rm[6];
        r[4][5][3] = x4*y5*z3*rm[12] + (6*x2*y5*z3+10*x4*y3*z3+3*x4*y5*z1)*rm[11] +
                     (3*y5*z3+60*x2*y3*z3+15*x4*y1*z3+18*x2*y5*z1
                      +30*x4*y3*z1)*rm[10] + (30*y3*z3+90*x2*y1*z3+9*y5*z1+180*x2*y3*z1+45*x4*y1*z1)*rm[9] +
                     (45*y1*z3+90*y3*z1+270*x2*y1*z1)*rm[8] + 135*y1*z1*rm[7];
        r[4][6][2] = x4*y6*z2*rm[12] + (6*x2*y6*z2+15*x4*y4*z2+x4*y6)*rm[11] +
                     (3*y6*z2+90*x2*y4*z2+45*x4*y2*z2+6*x2*y6
                      +15*x4*y4)*rm[10] + (45*y4*z2+270*x2*y2*z2+15*x4*z2+3*y6+90*x2*y4+45*x4*y2)*rm[9] +
                     (135*y2*z2+90*x2*z2+45*y4+270*x2*y2+15*x4)*rm[8] + (45*z2+135*y2+90*x2)*rm[7] + 45*rm[6];
        r[4][7][1] = x4*y7*z1*rm[12] + (6*x2*y7*z1+21*x4*y5*z1)*rm[11] + (3*y7*z1+126*x2*y5*z1+105*x4*y3*z1)
                     *rm[10] +
                     (63*y5*z1+630*x2*y3*z1+105*x4*y1*z1)*rm[9] + (315*y3*z1+630*x2*y1*z1)*rm[8] + 315*y1*z1*rm[7];
        r[4][8][0] = x4*y8*rm[12] + (6*x2*y8+28*x4*y6)*rm[11] + (3*y8+168*x2*y6+210*x4*y4)*rm[10] +
                     (84*y6+1260*x2*y4+420*x4*y2)*rm[9] + (630*y4+2520*x2*y2+105*x4)*rm[8] +
                     (1260*y2+630*x2)*rm[7] + 315*rm[6];
        r[5][0][7] = x5*z7*rm[12] + (10*x3*z7+21*x5*z5)*rm[11] + (15*x1*z7+210*x3*z5+105*x5*z3)*rm[10] +
                     (315*x1*z5+1050*x3*z3+105*x5*z1)*rm[9] + (1575*x1*z3+1050*x3*z1)*rm[8] + 1575*x1*z1*rm[7];
        r[5][1][6] = x5*y1*z6*rm[12] + (10*x3*y1*z6+15*x5*y1*z4)*rm[11] + (15*x1*y1*z6+150*x3*y1*z4
                     +45*x5*y1*z2)*rm[10] +
                     (225*x1*y1*z4+450*x3*y1*z2+15*x5*y1)*rm[9] + (675*x1*y1*z2+150*x3*y1)*rm[8] + 225*x1*y1*rm[7];
        r[5][2][5] = x5*y2*z5*rm[12] + (10*x3*y2*z5+x5*z5+10*x5*y2*z3)*rm[11] +
                     (15*x1*y2*z5+10*x3*z5+100*x3*y2*z3+10*x5*z3
                      +15*x5*y2*z1)*rm[10] + (15*x1*z5+150*x1*y2*z3+100*x3*z3+150*x3*y2*z1+15*x5*z1)*rm[9] +
                     (150*x1*z3+225*x1*y2*z1+150*x3*z1)*rm[8] + 225*x1*z1*rm[7];
        r[5][3][4] = x5*y3*z4*rm[12] + (10*x3*y3*z4+3*x5*y1*z4+6*x5*y3*z2)*rm[11] +
                     (15*x1*y3*z4+30*x3*y1*z4+60*x3*y3*z2
                      +18*x5*y1*z2+3*x5*y3)*rm[10] + (45*x1*y1*z4+90*x1*y3*z2+180*x3*y1*z2+30*x3*y3+9*x5*y1)*rm[9] +
                     (270*x1*y1*z2+45*x1*y3+90*x3*y1)*rm[8] + 135*x1*y1*rm[7];
        r[5][4][3] = x5*y4*z3*rm[12] + (10*x3*y4*z3+6*x5*y2*z3+3*x5*y4*z1)*rm[11] +
                     (15*x1*y4*z3+60*x3*y2*z3+3*x5*z3+30*x3*y4*z1
                      +18*x5*y2*z1)*rm[10] + (90*x1*y2*z3+30*x3*z3+45*x1*y4*z1+180*x3*y2*z1+9*x5*z1)*rm[9] +
                     (45*x1*z3+270*x1*y2*z1+90*x3*z1)*rm[8] + 135*x1*z1*rm[7];
        r[5][5][2] = x5*y5*z2*rm[12] + (10*x3*y5*z2+10*x5*y3*z2+x5*y5)*rm[11] +
                     (15*x1*y5*z2+100*x3*y3*z2+15*x5*y1*z2+10*x3*y5
                      +10*x5*y3)*rm[10] + (150*x1*y3*z2+150*x3*y1*z2+15*x1*y5+100*x3*y3+15*x5*y1)*rm[9] +
                     (225*x1*y1*z2+150*x1*y3+150*x3*y1)*rm[8] + 225*x1*y1*rm[7];
        r[5][6][1] = x5*y6*z1*rm[12] + (10*x3*y6*z1+15*x5*y4*z1)*rm[11] + (15*x1*y6*z1+150*x3*y4*z1
                     +45*x5*y2*z1)*rm[10] +
                     (225*x1*y4*z1+450*x3*y2*z1+15*x5*z1)*rm[9] + (675*x1*y2*z1+150*x3*z1)*rm[8] + 225*x1*z1*rm[7];
        r[5][7][0] = x5*y7*rm[12] + (10*x3*y7+21*x5*y5)*rm[11] + (15*x1*y7+210*x3*y5+105*x5*y3)*rm[10] +
                     (315*x1*y5+1050*x3*y3+105*x5*y1)*rm[9] + (1575*x1*y3+1050*x3*y1)*rm[8] + 1575*x1*y1*rm[7];
        r[6][0][6] = x6*z6*rm[12] + (15*x4*z6+15*x6*z4)*rm[11] + (45*x2*z6+225*x4*z4+45*x6*z2)*rm[10] +
                     (15*z6+675*x2*z4+675*x4*z2+15*x6)*rm[9] + (225*z4+2025*x2*z2+225*x4)*rm[8] +
                     (675*z2+675*x2)*rm[7] + 225*rm[6];
        r[6][1][5] = x6*y1*z5*rm[12] + (15*x4*y1*z5+10*x6*y1*z3)*rm[11] + (45*x2*y1*z5+150*x4*y1*z3
                     +15*x6*y1*z1)*rm[10] +
                     (15*y1*z5+450*x2*y1*z3+225*x4*y1*z1)*rm[9] + (150*y1*z3+675*x2*y1*z1)*rm[8] + 225*y1*z1*rm[7];
        r[6][2][4] = x6*y2*z4*rm[12] + (15*x4*y2*z4+x6*z4+6*x6*y2*z2)*rm[11] +
                     (45*x2*y2*z4+15*x4*z4+90*x4*y2*z2+6*x6*z2
                      +3*x6*y2)*rm[10] + (15*y2*z4+45*x2*z4+270*x2*y2*z2+90*x4*z2+45*x4*y2+3*x6)*rm[9] +
                     (15*z4+90*y2*z2+270*x2*z2+135*x2*y2
                      +45*x4)*rm[8] + (90*z2+45*y2+135*x2)*rm[7] + 45*rm[6];
        r[6][3][3] = x6*y3*z3*rm[12] + (15*x4*y3*z3+3*x6*y1*z3+3*x6*y3*z1)*rm[11] +
                     (45*x2*y3*z3+45*x4*y1*z3+45*x4*y3*z1
                      +9*x6*y1*z1)*rm[10] + (15*y3*z3+135*x2*y1*z3+135*x2*y3*z1+135*x4*y1*z1)*rm[9] +
                     (45*y1*z3+45*y3*z1+405*x2*y1*z1)*rm[8] +
                     135*y1*z1*rm[7];
        r[6][4][2] = x6*y4*z2*rm[12] + (15*x4*y4*z2+6*x6*y2*z2+x6*y4)*rm[11] +
                     (45*x2*y4*z2+90*x4*y2*z2+3*x6*z2+15*x4*y4
                      +6*x6*y2)*rm[10] + (15*y4*z2+270*x2*y2*z2+45*x4*z2+45*x2*y4+90*x4*y2+3*x6)*rm[9] +
                     (90*y2*z2+135*x2*z2+15*y4+270*x2*y2
                      +45*x4)*rm[8] + (45*z2+90*y2+135*x2)*rm[7] + 45*rm[6];
        r[6][5][1] = x6*y5*z1*rm[12] + (15*x4*y5*z1+10*x6*y3*z1)*rm[11] + (45*x2*y5*z1+150*x4*y3*z1
                     +15*x6*y1*z1)*rm[10] +
                     (15*y5*z1+450*x2*y3*z1+225*x4*y1*z1)*rm[9] + (150*y3*z1+675*x2*y1*z1)*rm[8] + 225*y1*z1*rm[7];
        r[6][6][0] = x6*y6*rm[12] + (15*x4*y6+15*x6*y4)*rm[11] + (45*x2*y6+225*x4*y4+45*x6*y2)*rm[10] +
                     (15*y6+675*x2*y4+675*x4*y2+15*x6)*rm[9] + (225*y4+2025*x2*y2+225*x4)*rm[8] +
                     (675*y2+675*x2)*rm[7] + 225*rm[6];
        r[7][0][5] = x7*z5*rm[12] + (21*x5*z5+10*x7*z3)*rm[11] + (105*x3*z5+210*x5*z3+15*x7*z1)*rm[10] +
                     (105*x1*z5+1050*x3*z3+315*x5*z1)*rm[9] + (1050*x1*z3+1575*x3*z1)*rm[8] + 1575*x1*z1*rm[7];
        r[7][1][4] = x7*y1*z4*rm[12] + (21*x5*y1*z4+6*x7*y1*z2)*rm[11] + (105*x3*y1*z4+126*x5*y1*z2+3*x7*y1)
                     *rm[10] +
                     (105*x1*y1*z4+630*x3*y1*z2+63*x5*y1)*rm[9] + (630*x1*y1*z2+315*x3*y1)*rm[8] + 315*x1*y1*rm[7];
        r[7][2][3] = x7*y2*z3*rm[12] + (21*x5*y2*z3+x7*z3+3*x7*y2*z1)*rm[11] +
                     (105*x3*y2*z3+21*x5*z3+63*x5*y2*z1+3*x7*z1)
                     *rm[10] + (105*x1*y2*z3+105*x3*z3+315*x3*y2*z1+63*x5*z1)*rm[9] + (105*x1*z3+315*x1*y2*z1+315*x3*z1)
                     *rm[8] +
                     315*x1*z1*rm[7];
        r[7][3][2] = x7*y3*z2*rm[12] + (21*x5*y3*z2+3*x7*y1*z2+x7*y3)*rm[11] +
                     (105*x3*y3*z2+63*x5*y1*z2+21*x5*y3+3*x7*y1)
                     *rm[10] + (105*x1*y3*z2+315*x3*y1*z2+105*x3*y3+63*x5*y1)*rm[9] + (315*x1*y1*z2+105*x1*y3+315*x3*y1)
                     *rm[8] +
                     315*x1*y1*rm[7];
        r[7][4][1] = x7*y4*z1*rm[12] + (21*x5*y4*z1+6*x7*y2*z1)*rm[11] + (105*x3*y4*z1+126*x5*y2*z1+3*x7*z1)
                     *rm[10] +
                     (105*x1*y4*z1+630*x3*y2*z1+63*x5*z1)*rm[9] + (630*x1*y2*z1+315*x3*z1)*rm[8] + 315*x1*z1*rm[7];
        r[7][5][0] = x7*y5*rm[12] + (21*x5*y5+10*x7*y3)*rm[11] + (105*x3*y5+210*x5*y3+15*x7*y1)*rm[10] +
                     (105*x1*y5+1050*x3*y3+315*x5*y1)*rm[9] + (1050*x1*y3+1575*x3*y1)*rm[8] + 1575*x1*y1*rm[7];
        r[8][0][4] = x8*z4*rm[12] + (28*x6*z4+6*x8*z2)*rm[11] + (210*x4*z4+168*x6*z2+3*x8)*rm[10] +
                     (420*x2*z4+1260*x4*z2+84*x6)*rm[9] + (105*z4+2520*x2*z2+630*x4)*rm[8] +
                     (630*z2+1260*x2)*rm[7] + 315*rm[6];
        r[8][1][3] = x8*y1*z3*rm[12] + (28*x6*y1*z3+3*x8*y1*z1)*rm[11] + (210*x4*y1*z3+84*x6*y1*z1)*rm[10] +
                     (420*x2*y1*z3+630*x4*y1*z1)*rm[9] + (105*y1*z3+1260*x2*y1*z1)*rm[8] + 315*y1*z1*rm[7];
        r[8][2][2] = x8*y2*z2*rm[12] + (28*x6*y2*z2+x8*z2+x8*y2)*rm[11] + (210*x4*y2*z2+28*x6*z2+28*x6*y2
                     +x8)*rm[10] +
                     (420*x2*y2*z2+210*x4*z2+210*x4*y2+28*x6)*rm[9] + (105*y2*z2+420*x2*z2+420*x2*y2+210*x4)*rm[8] +
                     (105*z2+105*y2+420*x2)*rm[7] + 105*rm[6];
        r[8][3][1] = x8*y3*z1*rm[12] + (28*x6*y3*z1+3*x8*y1*z1)*rm[11] + (210*x4*y3*z1+84*x6*y1*z1)*rm[10] +
                     (420*x2*y3*z1+630*x4*y1*z1)*rm[9] + (105*y3*z1+1260*x2*y1*z1)*rm[8] + 315*y1*z1*rm[7];
        r[8][4][0] = x8*y4*rm[12] + (28*x6*y4+6*x8*y2)*rm[11] + (210*x4*y4+168*x6*y2+3*x8)*rm[10] +
                     (420*x2*y4+1260*x4*y2+84*x6)*rm[9] + (105*y4+2520*x2*y2+630*x4)*rm[8] +
                     (630*y2+1260*x2)*rm[7] + 315*rm[6];
        r[9][0][3] = x9*z3*rm[12] + (36*x7*z3+3*x9*z1)*rm[11] + (378*x5*z3+108*x7*z1)*rm[10] +
                     (1260*x3*z3+1134*x5*z1)*rm[9] +
                     (945*x1*z3+3780*x3*z1)*rm[8] + 2835*x1*z1*rm[7];
        r[9][1][2] = x9*y1*z2*rm[12] + (36*x7*y1*z2+x9*y1)*rm[11] + (378*x5*y1*z2+36*x7*y1)*rm[10] +
                     (1260*x3*y1*z2+378*x5*y1)*rm[9] + (945*x1*y1*z2+1260*x3*y1)*rm[8] + 945*x1*y1*rm[7];
        r[9][2][1] = x9*y2*z1*rm[12] + (36*x7*y2*z1+x9*z1)*rm[11] + (378*x5*y2*z1+36*x7*z1)*rm[10] +
                     (1260*x3*y2*z1+378*x5*z1)*rm[9] + (945*x1*y2*z1+1260*x3*z1)*rm[8] + 945*x1*z1*rm[7];
        r[9][3][0] = x9*y3*rm[12] + (36*x7*y3+3*x9*y1)*rm[11] + (378*x5*y3+108*x7*y1)*rm[10] +
                     (1260*x3*y3+1134*x5*y1)*rm[9] +
                     (945*x1*y3+3780*x3*y1)*rm[8] + 2835*x1*y1*rm[7];
        r[10][0][2] = x10*z2*rm[12] + (45*x8*z2+x10)*rm[11] + (630*x6*z2+45*x8)*rm[10] +
                      (3150*x4*z2+630*x6)*rm[9] +
                      (4725*x2*z2+3150*x4)*rm[8] + (945*z2+4725*x2)*rm[7] + 945*rm[6];
        r[10][1][1] = x10*y1*z1*rm[12] + 45*x8*y1*z1*rm[11] + 630*x6*y1*z1*rm[10] + 3150*x4*y1*z1*rm[9] +
                      4725*x2*y1*z1*rm[8] +
                      945*y1*z1*rm[7];
        r[10][2][0] = x10*y2*rm[12] + (45*x8*y2+x10)*rm[11] + (630*x6*y2+45*x8)*rm[10] +
                      (3150*x4*y2+630*x6)*rm[9] +
                      (4725*x2*y2+3150*x4)*rm[8] + (945*y2+4725*x2)*rm[7] + 945*rm[6];
        r[11][0][1] = x11*z1*rm[12] + 55*x9*z1*rm[11] + 990*x7*z1*rm[10] + 6930*x5*z1*rm[9] +
                      17325*x3*z1*rm[8] +
                      10395*x1*z1*rm[7];
        r[11][1][0] = x11*y1*rm[12] + 55*x9*y1*rm[11] + 990*x7*y1*rm[10] + 6930*x5*y1*rm[9] +
                      17325*x3*y1*rm[8] +
                      10395*x1*y1*rm[7];
        r[12][0][0] = x12*rm[12] + 66*x10*rm[11] + 1485*x8*rm[10] + 13860*x6*rm[9] + 51975*x4*rm[8] +
                      62370*x2*rm[7] +
                      10395*rm[6];
        if ( ltot == 12) return;
        std::cerr << "Lmax exceeded in " << __FUNCTION__ <<"\n";
        exit(-1);
    }
};
}
#endif
