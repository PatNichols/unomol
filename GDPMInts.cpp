
#include "GDPMInts.hpp"

namespace unomol {

void calc_gdpm_ints(
    const ShellPairData& sq,
    double *svals,
    const AuxFunctions& aux,Rys& rys,const double *p) {
    double q[3];
    const double SRterm= 34.9868366552497250;
    const int npalf=6;
    const double palf[] = { 6.8505018000, 4.0491646000, 3.5941062989,
                            1.2478274000, 0.7927690989, 0.3377107978
                          };
    const double pcof[] = { 0.0766926784, 0.1483475912, 0.0462334641,
                            0.0717376427, 0.0447149792, 0.0069678529
                          };
    double ab[3],pzero[3],qc[3];
    int lzero[] = { 0, 0, 0 };
    lzero[0] = 0;
    lzero[1] = 0;
    lzero[2] = 0;
    pzero[0] = 0.;
    pzero[1] = 0.;
    pzero[2] = 0.;
    ab[0]=sq.a[0]-sq.b[0];
    ab[1]=sq.a[1]-sq.b[1];
    ab[2]=sq.a[2]-sq.b[2];
    const int lz=0;
    const int lvt12=sq.lv1+sq.lv2;
    const int nroots=lvt12/2+1;
    for (int i=0; i<npalf; ++i) {
        double pxp=palf[i];
        double c12=pcof[i];
        double abi=1.0/pxp;
        for (int k=0; k<sq.npr1; ++k) {
            double cxp=sq.al1[k];
            double c3=sq.co1[k];
            double f34=1.0;
            int lend=sq.npr2;
            if (sq.al1==sq.al2) {
                f34=2.0;
                lend=k+1;
            }
            for (int l=0; l<lend; ++l) {
                if (k==l) f34=1.0;
                double c34=c3*f34*sq.co2[l];
                double dxp=sq.al2[l];
                double qxp=cxp+dxp;
                double cdi=1.0/qxp;
                double s34=exp(-cxp*dxp*sq.ab2*cdi);
                double txp=pxp+qxp;
                double sr=SRterm*s34*cdi*abi/sqrt(txp);
                q[0]=(cxp*sq.a[0]+dxp*sq.b[0])*cdi;
                q[1]=(cxp*sq.a[1]+dxp*sq.b[1])*cdi;
                q[2]=(cxp*sq.a[2]+dxp*sq.b[2])*cdi;
                qc[0] = q[0] - sq.a[0];
                qc[1] = q[1] - sq.a[1];
                qc[2] = q[2] - sq.a[2];
//                rys.Recur(p,q,pa,qc,pxp,qxp,txp,lvt12,lvt34,nroots);
                rys.Recur(p,q,p,qc,pxp,qxp,txp,0,lvt12,nroots);
                for (int kc=0; kc<sq.len; ++kc) {
                    auto key=sq.lstates[kc];
                    int jls=key&0xF;
                    key>>=4;
                    int ils=key;
                    double nfact=c12*c34*sr*
                                 aux.normalization_factor(sq.lv1,ils)*
                                 aux.normalization_factor(sq.lv2,jls);
                    const int * lv1 = lzero;
                    const int * lv2 = lzero;
                    const int * lv3 = aux.l_vector(sq.lv1,ils);
                    const int * lv4 = aux.l_vector(sq.lv2,jls);
                    double sum= rys.Shift(
                                    pzero,ab,lv1,lv2,lv3,lv4,nroots);
                    svals[kc]+=sum*nfact;
                }
            }
        }
    }
}


void GDPMInts(const Basis& bas,double* Hmat) {
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
            calc_gdpm_ints(sp,svals,aux,rys,q);
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
