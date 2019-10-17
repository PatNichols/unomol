#ifndef _ONE_ELECTRON_INTS_CC_
#define _ONE_ELECTRON_INTS_CC_
#include <iostream>
#include <cmath>
#include "Basis.cc"
#include "AuxFunctions.cc"
#include "MD_Dfunction.cc"
#include "MD_Rfunction.cc"
#include "Util.cc"
#include "Structs.cc"
using namespace std;

inline  void calc_one_electron_ints(
    const ShellPairData& sp,
    double svals[],double tvals[],double vvals[],
    const Center* center,int ncen,int skip,
    const AuxFunctions& aux,
    MD_Dfunction& dx,MD_Dfunction& dy,MD_Dfunction& dz,
    MD_Rfunction& r,double*** rsum)
{
    register int ix,iy,iz;
    const double piterm=5.5683279968317079;
    double p[3];
    double pc[3];
    int lvt=sp.lv1+sp.lv2;
    for (int ip=0;ip<sp.npr1;ip++)
    {
        double axp=sp.al1[ip];
        double c1=sp.co1[ip];
        for (int jp=0;jp<sp.npr2;jp++)
        {
            double c12=c1*sp.co2[jp];
            double bxp=sp.al2[jp];
            double pxp=axp+bxp;
            double abi=1.0/pxp;
            double s12=piterm*exp(-axp*bxp*abi*sp.ab2)/(pxp*sqrt(pxp));
            p[0]=(axp*sp.a[0]+bxp*sp.b[0])*abi;
            p[1]=(axp*sp.a[1]+bxp*sp.b[1])*abi;
            p[2]=(axp*sp.a[2]+bxp*sp.b[2])*abi;
            abi=abi*0.5;
            dx.eval(abi,p[0]-sp.a[0],p[0]-sp.b[0],sp.lv1+1,sp.lv2+1);
            dy.eval(abi,p[1]-sp.a[1],p[1]-sp.b[1],sp.lv1+1,sp.lv2+1);
            dz.eval(abi,p[2]-sp.a[2],p[2]-sp.b[2],sp.lv1+1,sp.lv2+1);
            double sr=2.0*s12*sqrt(pxp/M_PI);
            for (ix=0;ix<=lvt;ix++)
            {
                for (iy=0;iy<=lvt;iy++)
                {
                    for (iz=0;iz<=lvt;iz++)
                    {
                        rsum[ix][iy][iz]=0.0;
                    }
                }
            }
            for (int ic=0;ic<ncen;ic++)
            {
                if (ic==skip) continue;
                double qc=(center+ic)->charge();
                const double* rc=(center+ic)->r_vec();
                pc[0]=p[0]-rc[0];
                pc[1]=p[1]-rc[1];
                pc[2]=p[2]-rc[2];
                double pc2=pc[0]*pc[0]+pc[1]*pc[1]+pc[2]*pc[2];
                r.eval(sr,(pxp*pc2),pxp,pc,lvt);
                for (ix=0;ix<=lvt;ix++)
                {
                    for (iy=0;iy<=lvt;iy++)
                    {
                        for (iz=0;iz<=lvt;iz++)
                        {
                            rsum[ix][iy][iz]-=qc*r.getValue(ix,iy,iz);
                        }
                    }
                }
            }
            for (int kc=0;kc<sp.len;kc++)
            {
                unsigned short key=sp.lstates[kc];
                int jls=key&0xF;
                int ils=key>>4;
                int l1=aux.Lxyz(sp.lv1,ils,0);
                int m1=aux.Lxyz(sp.lv1,ils,1);
                int n1=aux.Lxyz(sp.lv1,ils,2);
                int l2=aux.Lxyz(sp.lv2,jls,0);
                int m2=aux.Lxyz(sp.lv2,jls,1);
                int n2=aux.Lxyz(sp.lv2,jls,2);
                int l12=l1+l2;
                int m12=m1+m2;
                int n12=n1+n2;
                double nfct=aux.normalization_factor(sp.lv1,ils)*
                            aux.normalization_factor(sp.lv2,jls)*c12;
                svals[kc]+=nfct*s12*dx.getValue(l1,l2,0)*dy.getValue(m1,m2,0)*
                           dz.getValue(n1,n2,0);
                double tx=2.0*axp*bxp*dx.getValue(l1+1,l2+1,0);
                if (l1!=0)
                {
                    tx-=l1*bxp*dx.getValue(l1-1,l2+1,0);
                }
                if (l2!=0)
                {
                    tx-=l2*axp*dx.getValue(l1+1,l2-1,0);
                    if (l1!=0)
                    {
                        tx+=0.5*l1*l2*dx.getValue(l1-1,l2-1,0);
                    }
                }
                tx*=dy.getValue(m1,m2,0)*dz.getValue(n1,n2,0);
                double ty=2.0*axp*bxp*dy.getValue(m1+1,m2+1,0);
                if (m1!=0)
                {
                    ty-=m1*bxp*dy.getValue(m1-1,m2+1,0);
                }
                if (m2!=0)
                {
                    ty-=m2*axp*dy.getValue(m1+1,m2-1,0);
                    if (m1!=0)
                    {
                        ty+=0.5*m1*m2*dy.getValue(m1-1,m2-1,0);
                    }
                }
                ty*=dx.getValue(l1,l2,0)*dz.getValue(n1,n2,0);
                double tz=2.0*axp*bxp*dz.getValue(n1+1,n2+1,0);
                if (n1!=0)
                {
                    tz-=n1*bxp*dz.getValue(n1-1,n2+1,0);
                }
                if (n2!=0)
                {
                    tz-=n2*axp*dz.getValue(n1+1,n2-1,0);
                    if (n1!=0)
                    {
                        tz+=0.5*n1*n2*dz.getValue(n1-1,n2-1,0);
                    }
                }
                tz*=dx.getValue(l1,l2,0)*dy.getValue(m1,m2,0);
                tvals[kc]+=nfct*s12*(tx+ty+tz);
                double sum=0.0;
                for (ix=0;ix<=l12;ix++)
                {
                    double fx=dx.getValue(l1,l2,ix);
                    for (iy=0;iy<=m12;iy++)
                    {
                        double fxy=fx*dy.getValue(m1,m2,iy);
                        for (iz=0;iz<=n12;iz++)
                        {
                            sum+=dz.getValue(n1,n2,iz)*fxy*rsum[ix][iy][iz];
                        }
                    }
                }
                vvals[kc]+=nfct*sum;
            }
        }
    }
};

inline
void  OneElectronInts(const Basis& bas,double* Smat,
                      double* Tmat, double* Hmat)
{
    int maxl=bas.maxLvalue();
    MD_Dfunction dx(maxl+1);
    MD_Dfunction dy(maxl+1);
    MD_Dfunction dz(maxl+1);
    MD_Rfunction r(maxl);
    int nshell=bas.number_of_shells();
    const Shell* shell=bas.shell_ptr();
    const Center* center=bas.center_ptr();
    const AuxFunctions& aux(*bas.auxfun_ptr());
    int ncen=bas.number_of_centers();
    int skip=bas.skip_center();
    int ml2=aux.number_of_lstates( maxl );
    ml2=ml2*ml2;
    double *svals=new double[ml2];
    double *tvals=new double[ml2];
    double *vvals=new double[ml2];
    ShellPairData sp;
    sp.lstates=new unsigned short[ml2];
    int * ostates=new int[ml2];
    int rsize=4*maxl+1;
    double*** rsum=new_tensor3< double >(rsize,rsize,rsize);
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
            for (int ils=0;ils<nls1;ils++)
            {
                int ijr=(ir*(ir+1)/2)+jr0;
                int jr=jr0;
                int jend=nls2;
                if (ishell==jshell) jend=ils+1;
                for (int jls=0;jls<jend;jls++)
                {
                    if (jr>ir) break;
                    sp.lstates[knt]=(unsigned short)((ils<<4)+jls);
                    ostates[knt]=ijr;
                    svals[knt]=0.0;
                    tvals[knt]=0.0;
                    vvals[knt]=0.0;
                    ++knt;
                    ++jr;
                    ++ijr;
                }
                ++ir;
            }
            if (!knt) continue;
            sp.len=knt;
            calc_one_electron_ints(sp,svals,tvals,vvals,center,ncen,skip,
                                   aux,dx,dy,dz,r,rsum);
            for (int kc=0;kc<knt;kc++)
            {
                int ijr=ostates[kc];
                Smat[ijr]=svals[kc];
                Tmat[ijr]=tvals[kc];
                Hmat[ijr]=tvals[kc]+vvals[kc];
            }
        }
    }
    delete_tensor3< double >(rsum,rsize,rsize);
    delete [] ostates;
    delete [] sp.lstates;
    delete [] vvals;
    delete [] tvals;
    delete [] svals;
};
#endif
