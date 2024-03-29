#include "TwoElectronInts.hpp"

namespace unomol {

#define UNO_MASK 0xF
#define UNO_SHIFT 4U
#define UNO_SHIFT2 8U

//#define UNOMOL_MD_INTS

#ifdef UNOMOL_MD_INTS

void calc_two_electron_ints_one_cen(const ShellQuartet& sq,
                            const AuxFunctions& aux,
                            MDInts& mds,
                            TwoInts* sints) 
{
    const double SRterm= 34.9868366552497250;
    const double threshold=1.e-12;
    const int lvt12=sq.lv1+sq.lv2;
    const int lvt34=sq.lv3+sq.lv4;
    const int lvt=lvt12+lvt34;
    MD_Dfunction& dx12 = mds.dx12;
    MD_Dfunction& dx34 = mds.dx34;
    MD_Rfunction& rfun = mds.rfun;
//    double pq[3];
//    pq[0] = pq[1] = pq[2] = 0.0;
    for (int i=0; i<sq.npr1; ++i) {
        double axp=sq.al1[i];
        double c1=sq.co1[i];
        double f12=1.0;
        int jend=sq.npr2;
        if (sq.al1==sq.al2) {
            f12=2.0;
            jend=i+1;
        }
        for (int j=0; j<jend; ++j) {
            if (i==j) f12=1.0;
            double s12=c1*f12*sq.co2[j];
            double bxp=sq.al2[j];
            double pxp=axp+bxp;
            double abi=1.0/pxp;
            abi *= 0.5;
            dx12.eval_one_cen(abi,sq.lv1,sq.lv2);
            for (int k=0; k<sq.npr3; ++k) {
                double cxp=sq.al3[k];
                double c3=sq.co3[k];
                double f34=1.0;
                int lend=sq.npr4;
                if (sq.al3==sq.al4) {
                    f34=2.0;
                    lend=k+1;
                }
                for (int l=0; l<lend; ++l) {
                    if (k==l) f34=1.0;
                    double s34=c3*f34*sq.co4[l];
                    double dxp=sq.al4[l];
                    double qxp=cxp+dxp;
                    double cdi=1.0/qxp;
                    double txp=pxp+qxp;
                    double sr= SRterm*s12*s34*2.*abi*cdi/sqrt(txp);
                    cdi *= 0.5;
                    dx34.eval_one_cen(cdi,sq.lv3,sq.lv4);
                    double w = pxp * qxp / txp;
                    rfun.eval_one_cen(sr,w,lvt);
                    for (int kc=0; kc<sq.len; ++kc) {
                        unsigned int key=sq.lstates[kc];
                        unsigned int lls=key&UNO_MASK;
                        key>>=UNO_SHIFT;
                        unsigned int kls=key&UNO_MASK;
                        key>>=UNO_SHIFT;
                        unsigned int jls=key&UNO_MASK;
                        key>>=UNO_SHIFT;
                        unsigned int ils=key&UNO_MASK;
                        const int *lvc1 = aux.l_vector(sq.lv1,ils);
                        const int *lvc2 = aux.l_vector(sq.lv2,jls);
                        const int *lvc3 = aux.l_vector(sq.lv3,kls);
                        const int *lvc4 = aux.l_vector(sq.lv4,lls);

                        int l1 = lvc1[0];
                        int m1 = lvc1[1];
                        int n1 = lvc1[2];

                        int l2 = lvc2[0];
                        int m2 = lvc2[1];
                        int n2 = lvc2[2];

                        int l12 = l1 + l2;
                        int m12 = m1 + m2;
                        int n12 = n1 + n2;

                        int l3 = lvc3[0];
                        int m3 = lvc3[1];
                        int n3 = lvc3[2];

                        int l4 = lvc4[0];
                        int m4 = lvc4[1];
                        int n4 = lvc4[2];

                        int l34 = l3 + l4;
                        int m34 = m3 + m4;
                        int n34 = n3 + n4;

                        double sum = 0;

                        for (int ix12=0; ix12<=l12; ++ix12) {
                            for (int iy12=0; iy12<=m12; ++iy12) {
                                for (int iz12=0; iz12<=n12; ++iz12) {
                                    double v12 =  dx12.getValue(l1,l2,ix12) *
                                                  dx12.getValue(m1,m2,iy12) *
                                                  dx12.getValue(n1,n2,iz12);
                                    for (int ix34=0; ix34<=l34; ++ix34) {
                                        for (int iy34=0; iy34<=m34; ++iy34) {
                                            const double v34 = v12 * dx34.getValue(l3,l4,ix34) *
                                                               dx34.getValue(m3,m4,iy34);
                                            const double *rzp = rfun.getRow((ix12+ix34),(iy12+iy34)) + iz12;
                                            const double *dzp = dx34.getRow(n3,n4);
                                            int sx = ((ix34+iy34)%2) ? -1:1;
                                            for (int iz34=0; iz34<=n34; ++iz34) {
                                                sum += sx * v34 *
                                                       dzp[iz34] *
                                                       rzp[iz34];
                                                sx = -sx;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        (sints+kc)->val += sum * sq.norms[kc];
                    }
                }
            }
        }
    }
}

void calc_two_electron_ints_two_cen(const ShellQuartet& sq,
                            const AuxFunctions& aux,
                            MDInts& mds,
                            TwoInts* sints) {
    const double SRterm= 34.9868366552497250;
    const double threshold=1.e-12;
    const int lvt12=sq.lv1+sq.lv2;
    const int lvt34=sq.lv3+sq.lv4;
    const int lvt=lvt12+lvt34;
    MD_Dfunction& dx12 = mds.dx12;
    MD_Dfunction& dx34 = mds.dx34;
    MD_Rfunction& rfun = mds.rfun;
    double p[3];
    p[0] = sq.a[0];
    p[1] = sq.a[1];
    p[2] = sq.a[2];
    double q[3];
    q[0] = sq.c[0];
    q[1] = sq.c[1];
    q[2] = sq.c[2];
    double pq[3];
    pq[0] = p[0]-q[0];
    pq[1] = p[1]-q[1];
    pq[2] = p[2]-q[2];
    double pq2 = pq[0]*pq[0] + pq[1]*pq[1] + pq[2]*pq[2];
    for (int i=0; i<sq.npr1; ++i) {
        double axp=sq.al1[i];
        double c1=sq.co1[i];
        double f12=1.0;
        int jend=sq.npr2;
        if (sq.al1==sq.al2) {
            f12=2.0;
            jend=i+1;
        }
        for (int j=0; j<jend; ++j) {
            if (i==j) f12=1.0;
            double s12=c1*f12*sq.co2[j];
            double bxp=sq.al2[j];
            double pxp=axp+bxp;
            double abi=1.0/pxp;
            abi *= 0.5;
            dx12.eval_one_cen(abi,sq.lv1,sq.lv2);
            for (int k=0; k<sq.npr3; ++k) {
                double cxp=sq.al3[k];
                double c3=sq.co3[k];
                double f34=1.0;
                int lend=sq.npr4;
                if (sq.al3==sq.al4) {
                    f34=2.0;
                    lend=k+1;
                }
                for (int l=0; l<lend; ++l) {
                    if (k==l) f34=1.0;
                    double s34=c3*f34*sq.co4[l];
                    double dxp=sq.al4[l];
                    double qxp=cxp+dxp;
                    double cdi=1.0/qxp;
                    double txp=pxp+qxp;
                    double sr= SRterm*s12*s34*2.*abi*cdi/sqrt(txp);
                    cdi *= 0.5;
                    dx34.eval_one_cen(cdi,sq.lv3,sq.lv4);
                    double w = pxp * qxp / txp;
                    double t = w * pq2;
                    rfun.eval(sr,t,w,pq,lvt);
                    for (int kc=0; kc<sq.len; ++kc) {
                        unsigned int key=sq.lstates[kc];
                        unsigned int lls=key&UNO_MASK;
                        key>>=UNO_SHIFT;
                        unsigned int kls=key&UNO_MASK;
                        key>>=UNO_SHIFT;
                        unsigned int jls=key&UNO_MASK;
                        key>>=UNO_SHIFT;
                        unsigned int ils=key&UNO_MASK;
                        const int *lvc1 = aux.l_vector(sq.lv1,ils);
                        const int *lvc2 = aux.l_vector(sq.lv2,jls);
                        const int *lvc3 = aux.l_vector(sq.lv3,kls);
                        const int *lvc4 = aux.l_vector(sq.lv4,lls);

                        int l1 = lvc1[0];
                        int m1 = lvc1[1];
                        int n1 = lvc1[2];

                        int l2 = lvc2[0];
                        int m2 = lvc2[1];
                        int n2 = lvc2[2];

                        int l12 = l1 + l2;
                        int m12 = m1 + m2;
                        int n12 = n1 + n2;

                        int l3 = lvc3[0];
                        int m3 = lvc3[1];
                        int n3 = lvc3[2];

                        int l4 = lvc4[0];
                        int m4 = lvc4[1];
                        int n4 = lvc4[2];

                        int l34 = l3 + l4;
                        int m34 = m3 + m4;
                        int n34 = n3 + n4;

                        double sum = 0;

                        for (int ix12=0; ix12<=l12; ++ix12) {
                            for (int iy12=0; iy12<=m12; ++iy12) {
                                for (int iz12=0; iz12<=n12; ++iz12) {
                                    double v12 =  dx12.getValue(l1,l2,ix12) *
                                                  dx12.getValue(m1,m2,iy12) *
                                                  dx12.getValue(n1,n2,iz12);
                                    for (int ix34=0; ix34<=l34; ++ix34) {
                                        for (int iy34=0; iy34<=m34; ++iy34) {
                                            const double v34 = v12 * dx34.getValue(l3,l4,ix34) *
                                                               dx34.getValue(m3,m4,iy34);
                                            const double *rzp = rfun.getRow((ix12+ix34),(iy12+iy34)) + iz12;
                                            const double *dzp = dx34.getRow(n3,n4);
                                            int sx = ((ix34+iy34)%2) ? -1:1;
                                            for (int iz34=0; iz34<=n34; ++iz34) {
                                                sum += sx * v34 *
                                                       dzp[iz34] *
                                                       rzp[iz34];
                                                sx = -sx;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        (sints+kc)->val += sum * sq.norms[kc];
                    }
                }
            }
        }
    }
}

void calc_two_electron_ints(const ShellQuartet& sq,
                            const AuxFunctions& aux,
                            MDInts& mds,
                            TwoInts* sints) {
    double p[3],q[3],pq[3];
    if ( sq.a == sq.b && sq.a == sq.c && sq.a == sq.d) return calc_two_electron_ints_one_cen(sq,aux,mds,sints);
    if ( sq.a == sq.b && sq.c == sq.d) return calc_two_electron_ints_two_cen(sq,aux,mds,sints);
    const double SRterm= 34.9868366552497250;
    const double threshold=1.e-12;
    const int lvt12=sq.lv1+sq.lv2;
    const int lvt34=sq.lv3+sq.lv4;
    const int lvt=lvt12+lvt34;
    MD_Dfunction& dx12 = mds.dx12;
    MD_Dfunction& dy12 = mds.dy12;
    MD_Dfunction& dz12 = mds.dz12;
    MD_Dfunction& dx34 = mds.dx34;
    MD_Dfunction& dy34 = mds.dy34;
    MD_Dfunction& dz34 = mds.dz34;
    MD_Rfunction& rfun = mds.rfun;
    for (int i=0; i<sq.npr1; ++i) {
        double axp=sq.al1[i];
        double c1=sq.co1[i];
        double f12=1.0;
        int jend=sq.npr2;
        if (sq.al1==sq.al2) {
            f12=2.0;
            jend=i+1;
        }
        for (int j=0; j<jend; ++j) {
            if (i==j) f12=1.0;
            double c12=c1*f12*sq.co2[j];
            double bxp=sq.al2[j];
            double pxp=axp+bxp;
            double abi=1.0/pxp;
            double s12= c12 * std::exp(-axp*bxp*sq.ab2*abi);
            p[0]=(axp*sq.a[0]+bxp*sq.b[0])*abi;
            p[1]=(axp*sq.a[1]+bxp*sq.b[1])*abi;
            p[2]=(axp*sq.a[2]+bxp*sq.b[2])*abi;
            abi *= 0.5;
            dx12.eval(abi,(p[0]-sq.a[0]),(p[0]-sq.b[0]),sq.lv1,sq.lv2);
            dy12.eval(abi,(p[1]-sq.a[1]),(p[1]-sq.b[1]),sq.lv1,sq.lv2);
            dz12.eval(abi,(p[2]-sq.a[2]),(p[2]-sq.b[2]),sq.lv1,sq.lv2);
            for (int k=0; k<sq.npr3; ++k) {
                double cxp=sq.al3[k];
                double c3=sq.co3[k];
                double f34=1.0;
                int lend=sq.npr4;
                if (sq.al3==sq.al4) {
                    f34=2.0;
                    lend=k+1;
                }
                for (int l=0; l<lend; ++l) {
                    if (k==l) f34=1.0;
                    double c34=c3*f34*sq.co4[l];
                    double dxp=sq.al4[l];
                    double qxp=cxp+dxp;
                    double cdi=1.0/qxp;
                    double s34= c34 * std::exp(-cxp*dxp*sq.cd2*cdi);
                    double txp=pxp+qxp;
                    double sr= SRterm*s12*s34*2.*abi*cdi/sqrt(txp);

                    q[0]=(cxp*sq.c[0]+dxp*sq.d[0])*cdi;
                    q[1]=(cxp*sq.c[1]+dxp*sq.d[1])*cdi;
                    q[2]=(cxp*sq.c[2]+dxp*sq.d[2])*cdi;
                    
                    cdi *= 0.5;
                    dx34.eval(cdi,(q[0]-sq.c[0]),(q[0]-sq.d[0]),sq.lv3,sq.lv4);
                    dy34.eval(cdi,(q[1]-sq.c[1]),(q[1]-sq.d[1]),sq.lv3,sq.lv4);
                    dz34.eval(cdi,(q[2]-sq.c[2]),(q[2]-sq.d[2]),sq.lv3,sq.lv4);

                    pq[0] = p[0]-q[0];
                    pq[1] = p[1]-q[1];
                    pq[2] = p[2]-q[2];
                    double pq2 = pq[0]*pq[0] + pq[1]*pq[1] + pq[2]*pq[2];
                    double w = pxp * qxp / txp;
                    double t = w * pq2;

                    rfun.eval(sr,t,w,pq,lvt);

                    for (int kc=0; kc<sq.len; ++kc) {
                        unsigned int key=sq.lstates[kc];
                        unsigned int lls=key&UNO_MASK;
                        key>>=UNO_SHIFT;
                        unsigned int kls=key&UNO_MASK;
                        key>>=UNO_SHIFT;
                        unsigned int jls=key&UNO_MASK;
                        key>>=UNO_SHIFT;
                        unsigned int ils=key&UNO_MASK;
                        const int *lvc1 = aux.l_vector(sq.lv1,ils);
                        const int *lvc2 = aux.l_vector(sq.lv2,jls);
                        const int *lvc3 = aux.l_vector(sq.lv3,kls);
                        const int *lvc4 = aux.l_vector(sq.lv4,lls);

                        int l1 = lvc1[0];
                        int m1 = lvc1[1];
                        int n1 = lvc1[2];

                        int l2 = lvc2[0];
                        int m2 = lvc2[1];
                        int n2 = lvc2[2];

                        int l12 = l1 + l2;
                        int m12 = m1 + m2;
                        int n12 = n1 + n2;

                        int l3 = lvc3[0];
                        int m3 = lvc3[1];
                        int n3 = lvc3[2];

                        int l4 = lvc4[0];
                        int m4 = lvc4[1];
                        int n4 = lvc4[2];

                        int l34 = l3 + l4;
                        int m34 = m3 + m4;
                        int n34 = n3 + n4;

                        double sum = 0;

                        for (int ix12=0; ix12<=l12; ++ix12) {
                            for (int iy12=0; iy12<=m12; ++iy12) {
                                for (int iz12=0; iz12<=n12; ++iz12) {
                                    double v12 =  dx12.getValue(l1,l2,ix12) *
                                                  dy12.getValue(m1,m2,iy12) *
                                                  dz12.getValue(n1,n2,iz12);
                                    for (int ix34=0; ix34<=l34; ++ix34) {
                                        for (int iy34=0; iy34<=m34; ++iy34) {
                                            const double v34 = v12 * dx34.getValue(l3,l4,ix34) *
                                                               dy34.getValue(m3,m4,iy34);
                                            const double *rzp = rfun.getRow((ix12+ix34),(iy12+iy34)) + iz12;
                                            const double *dzp = dz34.getRow(n3,n4);
                                            int sx = ((ix34+iy34)%2) ? -1:1;
                                            for (int iz34=0; iz34<=n34; ++iz34) {
                                                sum += sx * v34 *
                                                       dzp[iz34] *
                                                       rzp[iz34];
                                                sx = -sx;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        (sints+kc)->val += sum * sq.norms[kc];
                    }
                }
            }
        }
    }
}

#else

void calc_two_electron_ints(const ShellQuartet& sq,
                            const AuxFunctions& aux,
                            Rys& rys,
                            TwoInts* sints) {
    double p[3],q[3];
    double ab[3],cd[3];
    double pa[3],qc[3];
    const double SRterm= 34.9868366552497250;
    const double threshold=1.e-12;
    ab[0]=sq.a[0]-sq.b[0];
    ab[1]=sq.a[1]-sq.b[1];
    ab[2]=sq.a[2]-sq.b[2];
    cd[0]=sq.c[0]-sq.d[0];
    cd[1]=sq.c[1]-sq.d[1];
    cd[2]=sq.c[2]-sq.d[2];
    const int lvt12=sq.lv1+sq.lv2;
    const int lvt34=sq.lv3+sq.lv4;
    const int lvt=lvt12+lvt34;
    const int nroots=lvt/2+1;
    for (int i=0; i<sq.npr1; ++i) {
        double axp=sq.al1[i];
        double c1=sq.co1[i];
        double f12=1.0;
        int jend=sq.npr2;
        if (sq.al1==sq.al2) {
            f12=2.0;
            jend=i+1;
        }
        for (int j=0; j<jend; ++j) {
            if (i==j) f12=1.0;
            double c12=c1*f12*sq.co2[j];
            double bxp=sq.al2[j];
            double pxp=axp+bxp;
            double abi=1.0/pxp;
            double s12= std::exp(-axp*bxp*sq.ab2*abi);
            p[0]=(axp*sq.a[0]+bxp*sq.b[0])*abi;
            p[1]=(axp*sq.a[1]+bxp*sq.b[1])*abi;
            p[2]=(axp*sq.a[2]+bxp*sq.b[2])*abi;
            pa[0] = p[0] - sq.a[0];
            pa[1] = p[1] - sq.a[1];
            pa[2] = p[2] - sq.a[2];
            for (int k=0; k<sq.npr3; ++k) {
                double cxp=sq.al3[k];
                double c3=sq.co3[k];
                double f34=1.0;
                int lend=sq.npr4;
                if (sq.al3==sq.al4) {
                    f34=2.0;
                    lend=k+1;
                }
                for (int l=0; l<lend; ++l) {
                    if (k==l) f34=1.0;
                    double c34=c3*f34*sq.co4[l];
                    double dxp=sq.al4[l];
                    double qxp=cxp+dxp;
                    double cdi=1.0/qxp;
                    double s34= std::exp(-cxp*dxp*sq.cd2*cdi);
                    double txp=pxp+qxp;
                    double sr=SRterm*s12*s34*abi*cdi/sqrt(txp);
                    if (sr<threshold) continue;
                    sr *= c12 * c34;
                    q[0]=(cxp*sq.c[0]+dxp*sq.d[0])*cdi;
                    q[1]=(cxp*sq.c[1]+dxp*sq.d[1])*cdi;
                    q[2]=(cxp*sq.c[2]+dxp*sq.d[2])*cdi;
                    qc[0] = q[0] - sq.c[0];
                    qc[1] = q[1] - sq.c[1];
                    qc[2] = q[2] - sq.c[2];
                    rys.Recur(p,q,pa,qc,pxp,qxp,txp,lvt12,lvt34,nroots);
                    for (int kc=0; kc<sq.len; ++kc) {
                        unsigned int key=sq.lstates[kc];
                        int lls=key&UNO_MASK;
                        key>>=UNO_SHIFT;
                        int kls=key&UNO_MASK;
                        key>>=UNO_SHIFT;
                        int jls=key&UNO_MASK;
                        key>>=UNO_SHIFT;
                        int ils=key&UNO_MASK;
                        const int *lv1 = aux.l_vector(sq.lv1,ils);
                        const int *lv2 = aux.l_vector(sq.lv2,jls);
                        const int *lv3 = aux.l_vector(sq.lv3,kls);
                        const int *lv4 = aux.l_vector(sq.lv4,lls);
                        double sum=
                            rys.Shift(ab,cd,lv1,lv2,lv3,lv4,nroots);
                        (sints+kc)->val+=sum*sr*sq.norms[kc];
                    }
                }
            }
        }
    }
}

#endif

void
TwoElectronInts::calculate(const Basis& basis) {
    const double threshold=1.e-14;
    int pknt=0;
    const Shell* shell(basis.shell_ptr());
    const Center* center(basis.center_ptr());
    const AuxFunctions& aux(*basis.auxfun_ptr());
    const int nshell=basis.number_of_shells();
    const int maxl=basis.maxLvalue();
    const int maxlst=aux.maxLstates();
    int ml2=maxlst*maxlst;
    int ml4=ml2*ml2;
    int ir0=(0);
    int nls1=(0);
    for (int i=0; i<start; ++i,ir0+=nls1) {
        int lv1=(shell+i)->Lvalue();
        nls1=aux.number_of_lstates(lv1);
    }
#ifdef UNOMOL_MD_INTS
    MDInts mds(maxl);
#else
    Rys rys(maxl);
#endif
    ShellQuartet sq(maxl);
    TwoInts* sints=new TwoInts[ml4];
    int it;
    const double *dp;
    cache.open_for_writing();
    putils::Stopwatch timer;
    timer.start();
    int ncalc = 0;
    for (int ish=start; ish<nshell; ++ish) {
        int ir0 = basis.offset(ish);
        sq.npr1=(shell+ish)->number_of_prims();
        sq.lv1=(shell+ish)->Lvalue();
        int cen1=(shell+ish)->center();
        sq.al1=(shell+ish)->alf_ptr();
        sq.co1=(shell+ish)->cof_ptr();
        sq.a=(center+cen1)->r_vec();
        int nls1=aux.number_of_lstates(sq.lv1);
        int lv1 = sq.lv1;
        for (int jsh=0; jsh<=ish; ++jsh) {
            int jr0 = basis.offset(jsh);
            sq.npr2=(shell+jsh)->number_of_prims();
            sq.lv2=(shell+jsh)->Lvalue();
            int cen2=(shell+jsh)->center();
            sq.al2=(shell+jsh)->alf_ptr();
            sq.co2=(shell+jsh)->cof_ptr();
            sq.b=(center+cen2)->r_vec();
            sq.ab2=dist_sqr(sq.a,sq.b);
            int nls2=aux.number_of_lstates(sq.lv2);
            int lv2 = sq.lv2;
            bool switch12=sq.lv1<sq.lv2;
            if (switch12) {
                it=sq.npr1;
                sq.npr1=sq.npr2;
                sq.npr2=it;
                it=sq.lv1;
                sq.lv1=sq.lv2;
                sq.lv2=it;
                dp=sq.al1;
                sq.al1=sq.al2;
                sq.al2=dp;
                dp=sq.co1;
                sq.co1=sq.co2;
                sq.co2=dp;
                dp=sq.a;
                sq.a=sq.b;
                sq.b=dp;
            }
            for (int ksh=0; ksh<=ish; ++ksh) {
                int kr0 = basis.offset(ksh);
                sq.npr3=(shell+ksh)->number_of_prims();
                sq.lv3=(shell+ksh)->Lvalue();
                int cen3=(shell+ksh)->center();
                sq.al3=(shell+ksh)->alf_ptr();
                sq.co3=(shell+ksh)->cof_ptr();
                sq.c=(center+cen3)->r_vec();
                int nls3=aux.number_of_lstates(sq.lv3);
                int lv3 = sq.lv3;
                for (int lsh=0; lsh<=ksh; ++lsh) {
                    int lr0 = basis.offset(lsh);
                    sq.npr4=(shell+lsh)->number_of_prims();
                    sq.lv4=(shell+lsh)->Lvalue();
                    int cen4=(shell+lsh)->center();
                    sq.al4=(shell+lsh)->alf_ptr();
                    sq.co4=(shell+lsh)->cof_ptr();
                    sq.d=(center+cen4)->r_vec();
                    sq.cd2=dist_sqr(sq.c,sq.d);
                    int nls4=aux.number_of_lstates(sq.lv4);
                    int lv4 = sq.lv4;
                    bool switch34=sq.lv3<sq.lv4;
                    if (switch34) {
                        it=sq.npr3;
                        sq.npr3=sq.npr4;
                        sq.npr4=it;
                        it=sq.lv3;
                        sq.lv3=sq.lv4;
                        sq.lv4=it;
                        dp=sq.al3;
                        sq.al3=sq.al4;
                        sq.al4=dp;
                        dp=sq.co3;
                        sq.co3=sq.co4;
                        sq.co4=dp;
                        dp=sq.c;
                        sq.c=sq.d;
                        sq.d=dp;
                    }
                    int knt=0;
                    for (int ils=0; ils<nls1;++ils) {
                        int ir = ir0 + ils;
                        for (int jls=0; jls<nls2;++jls) {
                            int jr = jr0 + jls;
                            if ( jr > ir ) break;
                            for (int kls=0; kls<nls3;++kls) {
                                int kr = kr0 + kls;
                                if ( kr > ir) break;
                                for (int lls=0; lls<nls4;++lls) {
                                    int lr = lr0 + lls;
                                    if ( lr > kr || ( ir == kr && lr > jr) ) break;
                                    (sints+knt)->val=0.0;
                                    (sints+knt)->i=(unsigned int)ir;
                                    (sints+knt)->j=(unsigned int)jr;
                                    (sints+knt)->k=(unsigned int)kr;
                                    (sints+knt)->l=(unsigned int)lr;
                                    unsigned int l12 = (ils<<UNO_SHIFT) + jls;
                                    if (switch12) l12=(jls<<UNO_SHIFT)+ils;
                                    unsigned int l34=(kls<<UNO_SHIFT)+lls;
                                    if (switch34) l34=(lls<<UNO_SHIFT)+kls;
                                    sq.lstates[knt]=(l12<<UNO_SHIFT2)+l34;
                                    sq.norms[knt] =
                                         aux.normalization_factor(lv1,ils)*
                                         aux.normalization_factor(lv2,jls)*
                                         aux.normalization_factor(lv3,kls)*
                                         aux.normalization_factor(lv4,lls);                     
                                    ++knt;
                                }
                            }
                        }
                    }
                    if (!knt) continue;
                    ncalc += knt;
                    sq.len=knt;
#ifdef UNOMOL_MD_INTS                    
                    calc_two_electron_ints(sq,aux,mds,sints);
#else
                    calc_two_electron_ints(sq,aux,rys,sints);
#endif
                    for (int kc=0; kc<knt; ++kc) {
                        if (fabs((sints+kc)->val)>threshold) {
                            cache.write(sints+kc,1);
                        }
                    }
                    if (switch34) {
                        sq.npr3=sq.npr4;
                        sq.lv3=sq.lv4;
                        sq.al3=sq.al4;
                        sq.co3=sq.co4;
                        sq.c=sq.d;
                    }
                }
            }
            if (switch12) {
                sq.npr1=sq.npr2;
                sq.lv1=sq.lv2;
                sq.al1=sq.al2;
                sq.co1=sq.co2;
                sq.a=sq.b;
            }
        }
    }
    timer.stop();
    cache.close();
    std::cerr << "Time for Two Electrons Integrals = " << timer.elapsed_time() << " seconds\n";
    size_t nb = cache.total_size()/sizeof(TwoInts);
    std::cerr << " # of write integrals = " << nb << "\n";
    std::cerr << " # of calc  integrals = " << ncalc << "\n";
    delete [] sints;
}

void
TwoElectronInts::formGmatrix(const double* Pmat,double *Gmat) {
    constexpr const int BINSIZE=8196;
    TwoInts sints[BINSIZE];

    cache.open_for_reading();
    size_t nints = cache.total_size()/sizeof(TwoInts);
    size_t nbin = nints/BINSIZE;
    size_t nextra = nints%BINSIZE;
    for (int ibin=0; ibin<nbin; ++ibin) {
        cache.read(sints,BINSIZE);
        for (int ix=0; ix<BINSIZE; ++ix) {
            double val=(sints+ix)->val;
            int i=(sints+ix)->i;
            int j=(sints+ix)->j;
            int k=(sints+ix)->k;
            int l=(sints+ix)->l;
            int ii=i*(i+1)/2;
            int ij=ii+j;
            int ik=ii+k;
            int il=ii+l;
            int jk,jl;
            int kk = ( k * ( k + 1 ) ) / 2;
            int kl = kk + l;
            if (j >= k) {
                int jj = (j * ( j + 1 ) ) / 2;
                jk = jj + k;
                jl = jj + l;
            } else {
                jk = kk + j;
                if (j>l) {
                    jl=j*(j+1)/2+l;
                } else {
                    jl=l*(l+1)/2+j;
                }
            }
            double da=val*2.0*Pmat[ij];
            double db=val*2.0*Pmat[kl];
            double sjl=val*Pmat[ik];
            double sjk=val*Pmat[il];
            double sik=val*Pmat[jl];
            double sil=val*Pmat[jk];
            if (k!=l) {
                db=db+db;
                Gmat[ik]-=sik;
                if (i!=j && j>=k) Gmat[jk]-=sjk;
            }
            Gmat[il]-=sil;
            Gmat[ij]+=db;
            if (i!=j && j>=l) Gmat[jl]-=sjl;
            if (ij!=kl) {
                if (i!=j) da=da+da;
                if (j<=k) {
                    Gmat[jk]-=sjk;
                    if (i==k && i!=j) Gmat[ik]-=sik;
                    if (k!=l && j<=l) Gmat[jl]-=sjl;
                }
                Gmat[kl]+=da;
            }
        }
    }
    if (nextra) {
        cache.read(sints,nextra);
        for (int ix=0; ix<nextra; ++ix) {
            double val=(sints+ix)->val;
            int i=(sints+ix)->i;
            int j=(sints+ix)->j;
            int k=(sints+ix)->k;
            int l=(sints+ix)->l;
            int ii=i*(i+1)/2;
            int ij=ii+j;
            int ik=ii+k;
            int il=ii+l;
            int jk,jl;
            int kk = ( k * ( k + 1 ) ) / 2;
            int kl = kk + l;
            if (j >= k) {
                int jj = (j * ( j + 1 ) ) / 2;
                jk = jj + k;
                jl = jj + l;
            } else {
                jk = kk + j;
                if (j>l) {
                    jl=j*(j+1)/2+l;
                } else {
                    jl=l*(l+1)/2+j;
                }
            }
            double da=val*2.0*Pmat[ij];
            double db=val*2.0*Pmat[kl];
            double sjl=val*Pmat[ik];
            double sjk=val*Pmat[il];
            double sik=val*Pmat[jl];
            double sil=val*Pmat[jk];
            if (k!=l) {
                db=db+db;
                Gmat[ik]-=sik;
                if (i!=j && j>=k) Gmat[jk]-=sjk;
            }
            Gmat[il]-=sil;
            Gmat[ij]+=db;
            if (i!=j && j>=l) Gmat[jl]-=sjl;
            if (ij!=kl) {
                if (i!=j) da=da+da;
                if (j<=k) {
                    Gmat[jk]-=sjk;
                    if (i==k && i!=j ) Gmat[ik]-=sik;
                    if (k!=l && j<=l ) Gmat[jl]-=sjl;
                }
                Gmat[kl]+=da;
            }
        } 
    }
    cache.close();
}

void
TwoElectronInts::formGmatrix(const double* PmatA,const double *PmatB,
                             double *GmatA,double *GmatB) {
    constexpr const int BINSIZE=8196;
    TwoInts sints[BINSIZE];

    cache.open_for_reading();
    size_t nints = cache.total_size()/sizeof(TwoInts);
    size_t nbin = nints/BINSIZE;
    size_t nextra = nints%BINSIZE;
    for (size_t ibin=0; ibin<nbin; ++ibin) {
        cache.read(sints,BINSIZE);
        for (int ix=0; ix<BINSIZE; ++ix) {
            double val=(sints+ix)->val;
            int i=(sints+ix)->i;
            int j=(sints+ix)->j;
            int k=(sints+ix)->k;
            int l=(sints+ix)->l;
            int ii=i*(i+1)/2;
            int ij=ii+j;
            int ik=ii+k;
            int il=ii+l;
            int jk,jl;
            int kl=k*(k+1)/2+l;
            if (j>k) {
                jk=j*(j+1)/2+k;
            } else {
                jk=k*(k+1)/2+j;
            }
            if (j>l) {
                jl=j*(j+1)/2+l;
            } else {
                jl=l*(l+1)/2+j;
            }
            double da=val*(PmatA[ij]+PmatB[ij]);
            double db=val*(PmatA[kl]+PmatB[kl]);
            double sjlA=val*PmatA[ik];
            double sjkA=val*PmatA[il];
            double sikA=val*PmatA[jl];
            double silA=val*PmatA[jk];
            double sjlB=val*PmatB[ik];
            double sjkB=val*PmatB[il];
            double sikB=val*PmatB[jl];
            double silB=val*PmatB[jk];
            if (k!=l) {
                db=db+db;
                GmatA[ik]-=sikA;
                GmatB[ik]-=sikB;
                if (i!=j && j>=k) {
                    GmatA[jk]-=sjkA;
                    GmatB[jk]-=sjkB;
                }
            }
            GmatA[il]-=silA;
            GmatA[ij]+=db;
            GmatB[il]-=silB;
            GmatB[ij]+=db;
            if (i!=j && j>=l) {
                GmatA[jl]-=sjlA;
                GmatB[jl]-=sjlB;
            }
            if (ij!=kl) {
                if (i!=j) da=da+da;
                if (j<=k) {
                    GmatA[jk]-=sjkA;
                    GmatB[jk]-=sjkB;
                    if (i!=j && i<=k) {
                        GmatA[ik]-=sikA;
                        GmatB[ik]-=sikB;
                    }
                    if (k!=l && j<=l) {
                        GmatA[jl]-=sjlA;
                        GmatB[jl]-=sjlB;
                    }
                }
                GmatA[kl]+=da;
                GmatB[kl]+=da;
            }
        }
    }
    if (nextra!=0) {
        cache.read(sints,nextra);
        for (int ix=0; ix<nextra; ++ix) {
            double val=(sints+ix)->val;
            int i=(sints+ix)->i;
            int j=(sints+ix)->j;
            int k=(sints+ix)->k;
            int l=(sints+ix)->l;
            int ii=i*(i+1)/2;
            int ij=ii+j;
            int ik=ii+k;
            int il=ii+l;
            int jk,jl;
            int kl=k*(k+1)/2+l;
            if (j>k) {
                jk=j*(j+1)/2+k;
            } else {
                jk=k*(k+1)/2+j;
            }
            if (j>l) {
                jl=j*(j+1)/2+l;
            } else {
                jl=l*(l+1)/2+j;
            }
            double da=val*(PmatA[ij]+PmatB[ij]);
            double db=val*(PmatA[kl]+PmatB[kl]);
            double sjlA=val*PmatA[ik];
            double sjkA=val*PmatA[il];
            double sikA=val*PmatA[jl];
            double silA=val*PmatA[jk];
            double sjlB=val*PmatB[ik];
            double sjkB=val*PmatB[il];
            double sikB=val*PmatB[jl];
            double silB=val*PmatB[jk];
            if (k!=l) {
                db=db+db;
                GmatA[ik]-=sikA;
                GmatB[ik]-=sikB;
                if (i!=j && j>=k) {
                    GmatA[jk]-=sjkA;
                    GmatB[jk]-=sjkB;
                }
            }
            GmatA[il]-=silA;
            GmatA[ij]+=db;
            GmatB[il]-=silB;
            GmatB[ij]+=db;
            if (i!=j && j>=l) {
                GmatA[jl]-=sjlA;
                GmatB[jl]-=sjlB;
            }
            if (ij!=kl) {
                if (i!=j) da=da+da;
                if (j<=k) {
                    GmatA[jk]-=sjkA;
                    GmatB[jk]-=sjkB;
                    if (i!=j && i<=k) {
                        GmatA[ik]-=sikA;
                        GmatB[ik]-=sikB;
                    }
                    if (k!=l && j<=l) {
                        GmatA[jl]-=sjlA;
                        GmatB[jl]-=sjlB;
                    }
                }
                GmatA[kl]+=da;
                GmatB[kl]+=da;
            }
        }
    }
    cache.close();
}

}
