#ifndef UNOMOL_GDPMINTS_HPP_
#define UNOMOL_GDPMINTS_HPP_
#include <iostream>
#include <cmath>
#include "Basis.hpp"
#include "AuxFunctions.hpp"
#include "Rys.hpp"
#include "Util.hpp"
#include "Structs.hpp"
using namespace std;

inline  void calc_gdpm_ints(
    const ShellPairData& sq,
    double *svals,
    const AuxFunctions& aux,Rys& rys,
    double ***Gx,double ***Gy,double ***Gz,const double* p)
{
    double q[3],pq[3],pq2;
    const double SRterm= 34.9868366552497250;
    const int npalf=6;
    const double palf[] = { 6.8505018000, 4.0491646000, 3.5941062989,
                            1.2478274000, 0.7927690989, 0.3377107978
                          };
    const double pcof[] = { 0.0766926784, 0.1483475912, 0.0462334641,
                            0.0717376427, 0.0447149792, 0.0069678529
                          };
    const double abx=sq.a[0]-sq.b[0];
    const double aby=sq.a[1]-sq.b[1];
    const double abz=sq.a[2]-sq.b[2];
    const double cdx=0.0;
    const int lz=0;
    const int lvt12=sq.lv1+sq.lv2;
    const int nroots=lvt12/2+1;
    for (int i=0;i<npalf;++i)
    {
        double pxp=palf[i];
        double c12=pcof[i];
        double abi=1.0/pxp;
        for (int k=0;k<sq.npr1;++k)
        {
            double cxp=sq.al1[k];
            double c3=sq.co1[k];
            double f34=1.0;
            int lend=sq.npr2;
            if (sq.al1==sq.al2)
            {
                f34=2.0;
                lend=k+1;
            }
            for (int l=0;l<lend;++l)
            {
                if (k==l) f34=1.0;
                double c34=c3*f34*sq.co2[l];
                double dxp=sq.al2[l];
                double qxp=cxp+dxp;
                double cdi=1.0/qxp;
                double s34=exp(-cxp*dxp*sq.ab2*cdi);
                double txp=pxp+qxp;
                double sr=SRterm*s34*cdi*abi/sqrt(txp);
                double w=pxp*qxp/txp;
                q[0]=(cxp*sq.a[0]+dxp*sq.b[0])*cdi;
                q[1]=(cxp*sq.a[1]+dxp*sq.b[1])*cdi;
                q[2]=(cxp*sq.a[2]+dxp*sq.b[2])*cdi;
                pq[0]=p[0]-q[0];
                pq[1]=p[1]-q[1];
                pq[2]=p[2]-q[2];
                pq2=pq[0]*pq[0]+pq[1]*pq[1]+pq[2]*pq[2];
                double t=pq2*w;
                rys.calculateRoots(nroots,t);
                for (register int iroot=0;iroot<nroots;++iroot)
                {
                    double rr=rys.Root(iroot);
                    double drt=rr/(1.0+rr);
                    double fff=drt/txp;
                    rys.B00 = 0.5 * fff;
                    rys.B1 = (0.5 - rys.B00 * qxp) / pxp;
                    rys.B1p = (0.5 - rys.B00 * pxp) / qxp;
                    rys.C = - qxp * (pq[0]) * fff;
                    rys.Cp = (q[0] - sq.a[0]) + pxp * (pq[0]) * fff;
                    rys.Recur (Gx[iroot], lz, lvt12);
                    rys.C =  - qxp * (pq[1]) * fff;
                    rys.Cp = (q[1] - sq.a[1]) + pxp * (pq[1]) * fff;
                    rys.Recur (Gy[iroot], lz, lvt12);
                    rys.C = - qxp * (pq[2]) * fff;
                    rys.Cp = (q[2] - sq.a[2]) + pxp * (pq[2]) * fff;
                    rys.Recur (Gz[iroot], lz, lvt12);
                }
                for (register int kc=0;kc<sq.len;++kc)
                {
                    unsigned short key=sq.lstates[kc];
                    int jls=key&0xF;
                    key>>=4;
                    int ils=key;
                    double nfact=c12*c34*sr*
                                 aux.normalization_factor(sq.lv1,ils)*
                                 aux.normalization_factor(sq.lv2,jls);
                    int l3=aux.Lxyz(sq.lv1,ils,0);
                    int m3=aux.Lxyz(sq.lv1,ils,1);
                    int n3=aux.Lxyz(sq.lv1,ils,2);
                    int l4=aux.Lxyz(sq.lv2,jls,0);
                    int m4=aux.Lxyz(sq.lv2,jls,1);
                    int n4=aux.Lxyz(sq.lv2,jls,2);
                    int l34=l3+l4;
                    int m34=m3+m4;
                    int n34=n3+n4;
                    register double sum=0.0;
                    for (register int iroot=0;iroot<nroots;++iroot)
                    {
                        sum+=rys.Shift(cdx,abx,Gx[iroot],lz,lz,l34,l4)*
                             rys.Shift(cdx,aby,Gy[iroot],lz,lz,m34,m4)*
                             rys.Shift(cdx,abz,Gz[iroot],lz,lz,n34,n4)*
                             rys.Weight(iroot);
                    }
                    svals[kc]+=sum*nfact;
                }
            }
        }
    }
}

inline
void GDPMInts(const Basis& bas,double* Hmat)
{
    int maxl=bas.maxLvalue();
    Rys rys(maxl);
    int nshell=bas.number_of_shells();
    const Shell* shell=bas.shell_ptr();
    const Center* center=bas.center_ptr();
    const AuxFunctions& aux(*bas.auxfun_ptr());
    int ml2=aux.number_of_lstates( maxl );
    ml2=ml2*ml2;
    const double *q=(center+bas.skip_center())->r_vec();
    int maxr=2*maxl+1;
    double ***Gx=new_tensor3<double>(maxr,maxr,maxr);
    double ***Gy=new_tensor3<double>(maxr,maxr,maxr);
    double ***Gz=new_tensor3<double>(maxr,maxr,maxr);
    double *svals=new double[ml2];
    ShellPairData sp;
    sp.lstates=new unsigned short[ml2];
    int * ostates=new int[ml2];
    int ir0=0;
    int nls1=0;
    for (int ishell=0;ishell<nshell;ir0+=nls1,ishell++)
    {
        sp.npr1=(shell+ishell)->number_of_prims();
        sp.lv1=(shell+ishell)->Lvalue();
        int cn1=(shell+ishell)->center();
        sp.al1=(shell+ishell)->alf_ptr();
        sp.co1=(shell+ishell)->cof_ptr();
        sp.a=(center+cn1)->r_vec();
        nls1=aux.number_of_lstates(sp.lv1);
        int jr0=0;
        int nls2=0;
        for (int jshell=0;jshell<=ishell;jr0+=nls2,jshell++)
        {
            sp.npr2=(shell+jshell)->number_of_prims();
            sp.lv2=(shell+jshell)->Lvalue();
            int cn2=(shell+jshell)->center();
            sp.al2=(shell+jshell)->alf_ptr();
            sp.co2=(shell+jshell)->cof_ptr();
            sp.b=(center+cn2)->r_vec();
            sp.ab2=dist_sqr(sp.a,sp.b);
            nls2=aux.number_of_lstates(sp.lv2);
            int knt=0;
            int ir=ir0;
            for (int ils=0;ils<nls1;++ir,ils++)
            {
                int ijr=(ir*(ir+1)/2)+jr0;
                int jend=nls2;
                if (ishell==jshell) jend=ils+1;
                for (int jls=0;jls<jend;++ijr,jls++)
                {
                    sp.lstates[knt]=(unsigned short)((ils<<4)+jls);
                    ostates[knt]=ijr;
                    svals[knt]=0.0;
                    ++knt;
                }
            }
            if (!knt) continue;
            sp.len=knt;
            calc_gdpm_ints(sp,svals,aux,rys,Gx,Gy,Gz,q);
            for (int kc=0;kc<knt;kc++)
            {
                int ijr=ostates[kc];
                Hmat[ijr]-=svals[kc];
            }
        }
    }
    delete [] ostates;
    delete [] sp.lstates;
    delete [] svals;
    delete_tensor3<double>(Gz,maxr,maxr);
    delete_tensor3<double>(Gy,maxr,maxr);
    delete_tensor3<double>(Gx,maxr,maxr);
};

#endif
