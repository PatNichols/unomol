#include "GDPMInts.hpp"

namespace unomol {

double calc_dpm_nucrep(
    const Center* center,int ncen,int skip)
{
    int ix,iy,iz;

    const double piterm=5.5683279968317079;
    double pc[3];
    const int npalf=6;
    const double palf[] = { 6.8505018000, 4.0491646000, 3.5941062989,
                            1.2478274000, 0.7927690989, 0.3377107978
                          };
    const double pcof[] = { 0.0766926784, 0.1483475912, 0.0462334641,
                            0.0717376427, 0.0447149792, 0.0069678529
                          };

    MD_Rfunction r(0);
    double rsum;
    double vval = 0.0;
    const double * p = center[skip].r_vec();
    for (int ip=0; ip<npalf; ++ip)
    {
        double pxp = palf[ip];
        double c12 = pcof[ip];
        double abi=1.0/pxp;
        double s12=piterm/(pxp*sqrt(pxp));
        double sr=2.0*s12*sqrt(pxp/M_PI);
        rsum = 0.0;
        for (int ic=0; ic<ncen; ic++) {
            if (ic==skip) continue;
            double qc=(center+ic)->charge();
            const double* rc=(center+ic)->r_vec();
            pc[0]=p[0]-rc[0];
            pc[1]=p[1]-rc[1];
            pc[2]=p[2]-rc[2];
            double pc2=pc[0]*pc[0]+pc[1]*pc[1]+pc[2]*pc[2];
            r.eval(sr,(pxp*pc2),pxp,pc,0);
            rsum += qc * r.getValue(0,0,0);
        }
        vval += rsum * c12;
    }
    return vval;
}


void calc_gdpm_ints(
    const ShellPairData& sp,
    double * vints,
    const AuxFunctions& aux,
    MD_Dfunction& dx,MD_Dfunction& dy,MD_Dfunction& dz,
    MD_Rfunction& rfun,
    const double *q) 
{
    const double SRterm= 34.9868366552497250;
    const int npalf=6;
    const double palf[] = { 6.8505018000, 4.0491646000, 3.5941062989,
                            1.2478274000, 0.7927690989, 0.3377107978
                          };
    const double pcof[] = { 0.0766926784, 0.1483475912, 0.0462334641,
                            0.0717376427, 0.0447149792, 0.0069678529
                          };

    const int lz=0;
    const int lvt12=sp.lv1+sp.lv2;
    int ix,iy,iz;

    double p[3];
    double pq[3];
    int lvt=sp.lv1+sp.lv2;
    for (int ip=0; ip<sp.npr1; ip++) {
        double axp=sp.al1[ip];
        double c1=sp.co1[ip];
        for (int jp=0; jp<sp.npr2; jp++) {
            double c12=c1*sp.co2[jp];
            double bxp=sp.al2[jp];
            double pxp=axp+bxp;
            double abi=1.0/pxp;
            double s12=exp(-axp*bxp*abi*sp.ab2);
            p[0]=(axp*sp.a[0]+bxp*sp.b[0])*abi;
            p[1]=(axp*sp.a[1]+bxp*sp.b[1])*abi;
            p[2]=(axp*sp.a[2]+bxp*sp.b[2])*abi;
            abi=abi*0.5;
            dx.eval(abi,p[0]-sp.a[0],p[0]-sp.b[0],sp.lv1,sp.lv2);
            dy.eval(abi,p[1]-sp.a[1],p[1]-sp.b[1],sp.lv1,sp.lv2);
            dz.eval(abi,p[2]-sp.a[2],p[2]-sp.b[2],sp.lv1,sp.lv2);
            for (int k=0;k<npalf;++k)
            {
                double qxp = palf[k];
                double txp = pxp + qxp;
                double w = pxp * qxp / txp;
                double cdi = 1./qxp;
                double c34 = pcof[k];
                pq[0] = p[0] - q[0];
                pq[1] = p[1] - q[1];
                pq[2] = p[2] - q[2];
                double pq2 = pq[0] * pq[0] + pq[1] * pq[1] + pq[2] * pq[2];
                double t = w * pq2;
                double sr = 2.*  c12 * c34 * SRterm * s12 * abi * cdi / sqrt (txp);
                rfun.eval(sr,t,w,pq,lvt12);
                for (int kc=0; kc<sp.len; ++kc) {
                    auto key=sp.lstates[kc];
                    const int jls=key&0xF;
                    key>>=4;
                    const int ils=key;
                    const int l1 = aux.Lxyz(sp.lv1,ils,0);
                    const int m1 = aux.Lxyz(sp.lv1,ils,1);
                    const int n1 = aux.Lxyz(sp.lv1,ils,2);
                    const int l2 = aux.Lxyz(sp.lv2,jls,0);
                    const int m2 = aux.Lxyz(sp.lv2,jls,1);
                    const int n2 = aux.Lxyz(sp.lv2,jls,2);
                    double nf = aux.normalization_factor(sp.lv1,ils) *
                            aux.normalization_factor(sp.lv2,jls);
                    const int l12 = l1 + l2;
                    const int m12 = m1 + m2;
                    const int n12 = n1 + n2; 
                    double sum = 0.0;
                    for (int ix=0;ix<=l12;++ix)
                    {  
                       const double cx = dx(l1,l2,ix);
                       for (int iy=0;iy<=m12;++iy)
                       {
                          const double cy=dy(m1,m2,iy)*cx;
                          const double * dzp= dz(n1,n2);
                          const double * rz = rfun.getRow(ix,iy);
                          for (int iz=0;iz<=n12;++iz)
                          {
                             sum+=dzp[iz]*rz[iz]*cy;
                          }
                       }
                    }
                    vints[kc] += sum * nf;
                }
            }
        }
    }
}

void GDPMInts(const Basis& bas,double* Hmat) {
    int maxl=bas.maxLvalue();
    MD_Dfunction dx(maxl);
    MD_Dfunction dy(maxl);
    MD_Dfunction dz(maxl);
    MD_Rfunction rfun(maxl);
    int nshell=bas.number_of_shells();
    const Shell* shell=bas.shell_ptr();
    const Center* center=bas.center_ptr();
    const AuxFunctions& aux(*bas.auxfun_ptr());
    int ml2=aux.number_of_lstates( maxl );
    ml2=ml2*ml2;
    const double *q=(center+bas.skip_center())->r_vec();
    int maxr=2*maxl+1;
    double *svals=new double[ml2];
    ShellPairData sp;
    sp.lstates=new unsigned short[ml2];
    int * ostates=new int[ml2];
    for (int ishell=0; ishell<nshell; ++ishell) {
        sp.npr1=(shell+ishell)->number_of_prims();
        sp.lv1=(shell+ishell)->Lvalue();
        int cn1=(shell+ishell)->center();
        sp.al1=(shell+ishell)->alf_ptr();
        sp.co1=(shell+ishell)->cof_ptr();
        sp.a=(center+cn1)->r_vec();
        int ir0 = bas.offset(ishell);
        int nls1=aux.number_of_lstates(sp.lv1);
        for (int jshell=0; jshell<=ishell; ++jshell) {
            sp.npr2=(shell+jshell)->number_of_prims();
            sp.lv2=(shell+jshell)->Lvalue();
            int cn2=(shell+jshell)->center();
            sp.al2=(shell+jshell)->alf_ptr();
            sp.co2=(shell+jshell)->cof_ptr();
            sp.b=(center+cn2)->r_vec();
            sp.ab2=dist_sqr(sp.a,sp.b);
            int jr0 = bas.offset(jshell);
            int nls2=aux.number_of_lstates(sp.lv2);
            int knt=0;
            for (int ils=0; ils<nls1; ++ils) {
                int ir = ir0 + ils;
                int iir = ir * (ir + 1 ) / 2;
                for (int jls=0; jls<nls2; ++jls) {
                    int jr = jr0 + jls;
                    if ( jr > ir ) break;
                    sp.lstates[knt]=((ils<<4)+jls);
                    ostates[knt]= iir + jr;
                    svals[knt]=0.0;
                    ++knt;
                }
            }
            if (!knt) continue;
            sp.len=knt;
            calc_gdpm_ints(sp,svals,
                aux,dx,dy,dz,rfun,
                q);
            for (int kc=0; kc<knt; kc++) {
                int ijr=ostates[kc];
                Hmat[ijr]-=svals[kc];
            }
        }
    }
    delete [] ostates;
    delete [] sp.lstates;
    delete [] svals;
};
}
