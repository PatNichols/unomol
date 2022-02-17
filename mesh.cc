#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
using namespace std;


inline void Cheb2ndKind(int npts,
                        double *rpts,double *wpts) {
    register int i;
    double cs,sn,t;
    const double arg0=M_PI/(npts+1.);
    const double sn0=sin(arg0);
    t=sin(0.5*arg0);
    t=t*t;
    t=t+t;
    const double cs0=t;
    cs=1.-cs0;
    sn=sn0;
    for (i=0; i<npts; ++i) {
        rpts[i]=cs;
        wpts[i]=arg0*sn;
        t=sn-sn*cs0+cs*sn0;
        cs=cs-cs*cs0-sn*sn0;
        sn=t;
    }
}

inline void Cheb2ndKindGen(int npts,
                           double *rpts,double *wpts) {
    register int i;
    double cs,sn,t;
    const double arg0=M_PI/(npts+1.);
    const double c1=2.0/M_PI;
    const double c2=2.0*c1/3.0;
    const double cw=16./3.0/(npts+1.);
    const double sn0=sin(arg0);
    t=sin(0.5*arg0);
    t=t*t;
    t=t+t;
    const double cs0=t;
    cs=1.-cs0;
    sn=sn0;
    for (i=0; i<npts; ++i) {
        double sn2=sn*sn;
        double sn4=sn2*sn2;
        rpts[i]=1.-2.*arg0*(i+1)+sn*cs*(c1+c2*sn2);
        wpts[i]=cw*sn4;
        t=sn-sn*cs0+cs*sn0;
        cs=cs-cs*cs0-sn*sn0;
        sn=t;
    }
}

inline void ClenshawCurtis(
