
#include "TwoElectronInts.hpp"

namespace unomol {

void TwoElectronInts::calc_two_electron_ints(const ShellQuartet& sq,
        const AuxFunctions& aux,
        Rys& rys,
        Sints* sints) {
    double p[3],q[3];
    double ab[3],cd[3];
    const double SRterm= 34.9868366552497250;
    const double threshold=1.e-15;
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
            double s12=std::exp(-axp*bxp*sq.ab2*abi);
            p[0]=(axp*sq.a[0]+bxp*sq.b[0])*abi;
            p[1]=(axp*sq.a[1]+bxp*sq.b[1])*abi;
            p[2]=(axp*sq.a[2]+bxp*sq.b[2])*abi;
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
                    double s34=std::exp(-cxp*dxp*sq.cd2*cdi);
                    double txp=pxp+qxp;
                    double sr=SRterm*s12*s34*abi*cdi/sqrt(txp);
                    if (sr<threshold) continue;
                    q[0]=(cxp*sq.c[0]+dxp*sq.d[0])*cdi;
                    q[1]=(cxp*sq.c[1]+dxp*sq.d[1])*cdi;
                    q[2]=(cxp*sq.c[2]+dxp*sq.d[2])*cdi;
                    rys.Recur(p,q,sq.a,sq.c,pxp,qxp,txp,lvt12,lvt34,nroots);
                    for (register int kc=0; kc<sq.len; ++kc) {
                        unsigned int key=sq.lstates[kc];
                        int lls=key&0xF;
                        key>>=4;
                        int kls=key&0xF;
                        key>>=4;
                        int jls=key&0xF;
                        key>>=4;
                        int ils=key&0xF;
                        double nfact=c12*c34*sr*
                                     aux.normalization_factor(sq.lv1,ils)*
                                     aux.normalization_factor(sq.lv2,jls)*
                                     aux.normalization_factor(sq.lv3,kls)*
                                     aux.normalization_factor(sq.lv4,lls);
                        const int *lv1 = aux.l_vector(sq.lv1,ils);
                        const int *lv2 = aux.l_vector(sq.lv2,jls);
                        const int *lv3 = aux.l_vector(sq.lv3,kls);
                        const int *lv4 = aux.l_vector(sq.lv4,lls);
                        double sum=
                            rys.Shift(ab,cd,lv1,lv2,lv3,lv4,nroots);
                        (sints+kc)->val+=sum*nfact;
                    }
                }
            }
        }
    }
}


void
TwoElectronInts::calculate(const Basis& basis) {
    const double threshold=1.e-14;
    int pknt=0;
    numints.clear();
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
    nfiles=0;
    Rys rys(maxl);
    ShellQuartet sq(maxl);
    Sints* sints=new Sints[ml4];
    int it;
    const double *dp;
    FILE *out=create_ints_file(0);
    int maxfsize=MAXFILESIZE/sizeof(Sints);
    int fsize=0;
    clock_t s=clock();
    for (int ish=start; ish<nshell; ++ish,ir0+=nls1) {
        sq.npr1=(shell+ish)->number_of_prims();
        sq.lv1=(shell+ish)->Lvalue();
        int cen1=(shell+ish)->center();
        sq.al1=(shell+ish)->alf_ptr();
        sq.co1=(shell+ish)->cof_ptr();
        sq.a=(center+cen1)->r_vec();
        nls1=aux.number_of_lstates(sq.lv1);
        int jr0=0;
        int nls2=0;
        for (int jsh=0; jsh<=ish; ++jsh,jr0+=nls2) {
            sq.npr2=(shell+jsh)->number_of_prims();
            sq.lv2=(shell+jsh)->Lvalue();
            int cen2=(shell+jsh)->center();
            sq.al2=(shell+jsh)->alf_ptr();
            sq.co2=(shell+jsh)->cof_ptr();
            sq.b=(center+cen2)->r_vec();
            sq.ab2=dist_sqr(sq.a,sq.b);
            nls2=aux.number_of_lstates(sq.lv2);
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
            int kr0=0;
            int nls3=0;
            for (int ksh=0; ksh<=ish; ++ksh,kr0+=nls3) {
                sq.npr3=(shell+ksh)->number_of_prims();
                sq.lv3=(shell+ksh)->Lvalue();
                int cen3=(shell+ksh)->center();
                sq.al3=(shell+ksh)->alf_ptr();
                sq.co3=(shell+ksh)->cof_ptr();
                sq.c=(center+cen3)->r_vec();
                nls3=aux.number_of_lstates(sq.lv3);
                int lr0=0;
                int nls4=0;
                for (int lsh=0; lsh<=ksh; ++lsh,lr0+=nls4) {
                    sq.npr4=(shell+lsh)->number_of_prims();
                    sq.lv4=(shell+lsh)->Lvalue();
                    int cen4=(shell+lsh)->center();
                    sq.al4=(shell+lsh)->alf_ptr();
                    sq.co4=(shell+lsh)->cof_ptr();
                    sq.d=(center+cen4)->r_vec();
                    sq.cd2=dist_sqr(sq.c,sq.d);
                    nls4=aux.number_of_lstates(sq.lv4);
#ifdef UNOMOL_MPI_API_
                    ++pknt;
                    pknt%=psize;
                    if (pknt!=rank) continue;
#endif
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
                    int ir=ir0;
                    for (register int ils=0; ils<nls1; ++ils,++ir) {
                        int ijr=ir*(ir+1)/2+jr0;
                        int jr=jr0;
                        int jend=nls2;
                        if (ish==jsh) jend=ils+1;
                        for (register int jls=0; jls<jend; ++ijr,++jr,++jls) {
                            int kr=kr0;
                            int kend=nls3;
                            if (ish==ksh) kend=ils+1;
                            for (register int kls=0; kls<kend; ++kls,++kr) {
                                int klr=kr*(kr+1)/2+lr0;
                                int lr=lr0;
                                int lend=nls4;
                                if (ksh==lsh) lend=kls+1;
                                for (register int lls=0; lls<lend; ++lr,++klr,++lls) {
                                    if (klr>ijr) break;
                                    (sints+knt)->val=0.0;
                                    (sints+knt)->i=(unsigned int)ir;
                                    (sints+knt)->j=(unsigned int)jr;
                                    (sints+knt)->k=(unsigned int)kr;
                                    (sints+knt)->l=(unsigned int)lr;
                                    unsigned int l12{0U};
                                    unsigned int l34{0U};
                                    l12 = (ils<<4) + jls;
                                    if (switch12) l12=(jls<<4)+ils;
                                    l34=(kls<<4)+lls;
                                    if (switch34) l34=(lls<<4)+kls;
                                    sq.lstates[knt] = 0U;
                                    sq.lstates[knt]=(l12<<8)+l34;
                                    ++knt;
                                }
                            }
                        }
                    }
                    if (!knt) continue;
                    sq.len=knt;
                    calc_two_electron_ints(sq,aux,rys,sints);
                    for (register int kc=0; kc<knt; ++kc) {
                        if (fabs((sints+kc)->val)>threshold) {
                            fwrite(sints+kc,sizeof(Sints),1,out);
                            ++fsize;
                        }
                    }
                    if (fsize>maxfsize) {
                        fclose(out);
                        numints.push_back(fsize);
                        ++nfiles;
                        out=create_ints_file(nfiles);
                        fsize=0;
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
    clock_t f=clock();
    double etime=((double)f-s)/CLOCKS_PER_SEC;
    fprintf(stderr,"Time for Two Electrons Integrals = %10.2le seconds\n",etime);
    if (fsize) {
        fclose(out);
        numints.push_back(fsize);
        ++nfiles;
    }
    delete [] sints;
}

void
TwoElectronInts::formGmatrix(const double* Pmat,double *Gmat) const {
    const int BINSIZE=1000;
    Sints sints[BINSIZE];

    for (int ifile=0; ifile<nfiles; ++ifile) {
        int nints=numints[ifile];
        int nbin=nints/BINSIZE;
        int nex=nints-BINSIZE*nbin;
        FILE *in=open_ints_file(ifile);
        if (nex) {
            fread(sints,sizeof(Sints),nex,in);
            for (register int ix=0; ix<nex; ++ix) {
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
                        if (i!=j && i<=k) Gmat[ik]-=sik;
                        if (k!=l && j<=l) Gmat[jl]-=sjl;
                    }
                    Gmat[kl]+=da;
                }
            }
        }
        for (int ibin=0; ibin<nbin; ++ibin) {
            fread(sints,sizeof(Sints),BINSIZE,in);
            for (register int ix=0; ix<BINSIZE; ++ix) {
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
                        if (i!=j && i<=k) Gmat[ik]-=sik;
                        if (k!=l && j<=l) Gmat[jl]-=sjl;
                    }
                    Gmat[kl]+=da;
                }
            }
        }
        fclose(in);
    }
}

void
TwoElectronInts::formGmatrix(const double* PmatA,const double *PmatB,
                             double *GmatA,double *GmatB) const {
    const int BINSIZE=1000;
    Sints sints[BINSIZE];

    for (int ifile=0; ifile<nfiles; ++ifile) {
        int nints=numints[ifile];
        int nbin=nints/BINSIZE;
        int nex=nints-BINSIZE*nbin;
        FILE *in=open_ints_file(ifile);
        if (nex) {
            fread(sints,sizeof(Sints),nex,in);
            for (register int ix=0; ix<nex; ++ix) {
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
        for (int ibin=0; ibin<nbin; ++ibin) {
            fread(sints,sizeof(Sints),BINSIZE,in);
            for (register int ix=0; ix<BINSIZE; ++ix) {
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
        fclose(in);
    }
}

}
