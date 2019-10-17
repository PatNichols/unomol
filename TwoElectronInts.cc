#ifndef _TWOELEC_CC_
#define _TWOELEC_CC_
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "pconfig.h"
#include "Util.cc"
#include "Basis.cc"
#include "Rys.cc"
#ifdef _MPI_API_
#include <mpi.h>
#endif
#define MAXFILESIZE 134217728

struct Sints
{
    double val;
    unsigned short i,j,k,l;
};

struct ShellQuartet
{
    double ab2,cd2;
    const double *a,*b,*c,*d;
    int npr1,lv1;
    const double *al1,*co1;
    int npr2,lv2;
    const double *al2,*co2;
    int npr3,lv3;
    const double *al3,*co3;
    int npr4,lv4;
    const double *al4,*co4;
    unsigned short *lstates;
    unsigned int len;
};

#ifdef _64BIT_API_
inline unsigned long int find_max_files(int no_new,int no_old,int nprocs)
{
    unsigned long int non=no_new;
    non=non*(non+1UL)/2UL;
    non=non*(non+1UL)/2UL;
    unsigned long int noo=no_old;
    noo=noo*(noo+1UL)/2UL;
    noo=noo*(noo+1UL)/2UL;
    non=non-noo;
    non*=sizeof(Sints);
    non/=MAXFILESIZE;
    non/=nprocs;
    ++non;
    return (unsigned long int)non;
}
#else
inline unsigned long int find_max_files(int no_new,int no_old,int nprocs)
{
    unsigned long long int non=no_new;
    non=non*(non+1ULL)/2ULL;
    non=non*(non+1ULL)/2ULL;
    unsigned long long int noo=no_old;
    noo=noo*(noo+1ULL)/2ULL;
    noo=noo*(noo+1ULL)/2ULL;
    non=non-noo;
    non*=sizeof(Sints);
    non/=MAXFILESIZE;
    non/=nprocs;
    ++non;
    return static_cast<unsigned long>(non);
}
#endif

class TwoElectronInts
{
public:
    TwoElectronInts(const Basis& base,int start_shell,const string& base_str);
    void recalculate(const Basis& base);
    void formGmatrix(const double* Pmat,double *Gmat) const;
    void formGmatrix(const double* PmatA,const double* PmatB,
                     double* GmatA,double* GmatB) const;
private:
    char basename[72];
    int start;
    int rank,psize;
    int *numints;
    int nfiles;
    void calc_two_electron_ints(const ShellQuartet& sq,
                                const AuxFunctions& aux,
                                Rys& rys,
                                Sints* sints,
                                double ***Gx,
                                double ***Gy,
                                double ***Gz);
    TwoElectronInts();
};

inline void TwoElectronInts::calc_two_electron_ints(const ShellQuartet& sq,
        const AuxFunctions& aux,
        Rys& rys,
        Sints* sints,
        double ***Gx,
        double ***Gy,
        double ***Gz)
{
    double p[3],q[3],pq[3],pq2;
    const double SRterm= 34.9868366552497250;
    const double threshold=1.e-15;
    const double abx=sq.a[0]-sq.b[0];
    const double aby=sq.a[1]-sq.b[1];
    const double abz=sq.a[2]-sq.b[2];
    const double cdx=sq.c[0]-sq.d[0];
    const double cdy=sq.c[1]-sq.d[1];
    const double cdz=sq.c[2]-sq.d[2];
    const int lvt12=sq.lv1+sq.lv2;
    const int lvt34=sq.lv3+sq.lv4;
    const int lvt=lvt12+lvt34;
    const int nroots=lvt/2+1;
    for (int i=0;i<sq.npr1;++i)
    {
        double axp=sq.al1[i];
        double c1=sq.co1[i];
        double f12=1.0;
        int jend=sq.npr2;
        if (sq.al1==sq.al2)
        {
            f12=2.0;
            jend=i+1;
        }
        for (int j=0;j<jend;++j)
        {
            if (i==j) f12=1.0;
            double c12=c1*f12*sq.co2[j];
            double bxp=sq.al2[j];
            double pxp=axp+bxp;
            double abi=1.0/pxp;
            double s12=std::exp(-axp*bxp*sq.ab2*abi);
            p[0]=(axp*sq.a[0]+bxp*sq.b[0])*abi;
            p[1]=(axp*sq.a[1]+bxp*sq.b[1])*abi;
            p[2]=(axp*sq.a[2]+bxp*sq.b[2])*abi;
            for (int k=0;k<sq.npr3;++k)
            {
                double cxp=sq.al3[k];
                double c3=sq.co3[k];
                double f34=1.0;
                int lend=sq.npr4;
                if (sq.al3==sq.al4)
                {
                    f34=2.0;
                    lend=k+1;
                }
                for (int l=0;l<lend;++l)
                {
                    if (k==l) f34=1.0;
                    double c34=c3*f34*sq.co4[l];
                    double dxp=sq.al4[l];
                    double qxp=cxp+dxp;
                    double cdi=1.0/qxp;
                    double s34=std::exp(-cxp*dxp*sq.cd2*cdi);
                    double txp=pxp+qxp;
                    double sr=SRterm*s12*s34*abi*cdi/sqrt(txp);
                    if (sr<threshold) continue;
                    double w=pxp*qxp/txp;
                    q[0]=(cxp*sq.c[0]+dxp*sq.d[0])*cdi;
                    q[1]=(cxp*sq.c[1]+dxp*sq.d[1])*cdi;
                    q[2]=(cxp*sq.c[2]+dxp*sq.d[2])*cdi;
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
                        rys.C = (p[0]-sq.a[0]) - qxp * (pq[0]) * fff;
                        rys.Cp = (q[0] - sq.c[0]) + pxp * (pq[0]) * fff;
                        rys.Recur (Gx[iroot], lvt12, lvt34);
                        rys.C = (p[1] - sq.a[1]) - qxp * (pq[1]) * fff;
                        rys.Cp = (q[1] - sq.c[1]) + pxp * (pq[1]) * fff;
                        rys.Recur (Gy[iroot], lvt12, lvt34);
                        rys.C = (p[2] - sq.a[2]) - qxp * (pq[2]) * fff;
                        rys.Cp = (q[2] - sq.c[2]) + pxp * (pq[2]) * fff;
                        rys.Recur (Gz[iroot], lvt12, lvt34);
                    }
                    for (register int kc=0;kc<sq.len;++kc)
                    {
                        unsigned short key=sq.lstates[kc];
                        int lls=key&0xF;
                        key>>=4;
                        int kls=key&0xF;
                        key>>=4;
                        int jls=key&0xF;
                        key>>=4;
                        int ils=key;
                        double nfact=c12*c34*sr*
                                     aux.normalization_factor(sq.lv1,ils)*
                                     aux.normalization_factor(sq.lv2,jls)*
                                     aux.normalization_factor(sq.lv3,kls)*
                                     aux.normalization_factor(sq.lv4,lls);
                        int l1=aux.Lxyz(sq.lv1,ils,0);
                        int m1=aux.Lxyz(sq.lv1,ils,1);
                        int n1=aux.Lxyz(sq.lv1,ils,2);
                        int l2=aux.Lxyz(sq.lv2,jls,0);
                        int m2=aux.Lxyz(sq.lv2,jls,1);
                        int n2=aux.Lxyz(sq.lv2,jls,2);
                        int l3=aux.Lxyz(sq.lv3,kls,0);
                        int m3=aux.Lxyz(sq.lv3,kls,1);
                        int n3=aux.Lxyz(sq.lv3,kls,2);
                        int l4=aux.Lxyz(sq.lv4,lls,0);
                        int m4=aux.Lxyz(sq.lv4,lls,1);
                        int n4=aux.Lxyz(sq.lv4,lls,2);
                        int l12=l1+l2;
                        int m12=m1+m2;
                        int n12=n1+n2;
                        int l34=l3+l4;
                        int m34=m3+m4;
                        int n34=n3+n4;
                        register double sum=0.0;
                        for (register int iroot=0;iroot<nroots;++iroot)
                        {
                            sum+=rys.Shift(abx,cdx,Gx[iroot],l12,l2,l34,l4)*
                                 rys.Shift(aby,cdy,Gy[iroot],m12,m2,m34,m4)*
                                 rys.Shift(abz,cdz,Gz[iroot],n12,n2,n34,n4)*
                                 rys.Weight(iroot);
                        }
                        (sints+kc)->val+=sum*nfact;
                    }
                }
            }
        }
    }
}

TwoElectronInts::TwoElectronInts(const Basis& basis,
                                 int start_shell,const string& base_str)
{
    const double threshold=1.e-14;
    size_t nstrlen = base_str.size() + 1;
    strncpy(basename,base_str.c_str(),nstrlen);
    start=start_shell;
    rank=0;
    psize=1;
#ifdef _MPI_API_
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&psize);
    int pknt=0;
#endif    
    char fname[80];
    const Shell* shell(basis.shell_ptr());
    const Center* center(basis.center_ptr());
    const AuxFunctions& aux(*basis.auxfun_ptr());
    const int nshell=basis.number_of_shells();
    const int maxl=basis.maxLvalue();
    const int maxlst=aux.maxLstates();
    int ml2=maxlst*maxlst;
    int ml4=ml2*ml2;
    int ir0=0;
    int nls1=0;
    for (int i=0;i<start;++i,ir0+=nls1)
    {
        int lv1=(shell+i)->Lvalue();
        nls1=aux.number_of_lstates(lv1);
    }
    unsigned long asize=find_max_files(basis.number_of_orbitals(),ir0,psize);
    numints=new int[asize];
    nfiles=0;
    Rys rys(maxl);
    int maxr=2*maxl+1;
    maxr=(maxr>2)?maxr:2;
    double ***Gx=new_tensor3<double>(maxr,maxr,maxr);
    double ***Gy=new_tensor3<double>(maxr,maxr,maxr);
    double ***Gz=new_tensor3<double>(maxr,maxr,maxr);
    Sints* sints=new Sints[ml4];
    ShellQuartet sq;
    sq.lstates=new unsigned short[ml4];
    make_filename(nfiles,rank,basename,fname);
    int it;
    const double *dp;
    FILE *out=create_file(fname);
    int maxfsize=MAXFILESIZE/sizeof(Sints);
    int fsize=0;
    clock_t s=clock();
    for (int ish=start;ish<nshell;++ish,ir0+=nls1)
    {
        sq.npr1=(shell+ish)->number_of_prims();
        sq.lv1=(shell+ish)->Lvalue();
        int cen1=(shell+ish)->center();
        sq.al1=(shell+ish)->alf_ptr();
        sq.co1=(shell+ish)->cof_ptr();
        sq.a=(center+cen1)->r_vec();
        nls1=aux.number_of_lstates(sq.lv1);
        int jr0=0;
        int nls2=0;
        for (int jsh=0;jsh<=ish;++jsh,jr0+=nls2)
        {
            sq.npr2=(shell+jsh)->number_of_prims();
            sq.lv2=(shell+jsh)->Lvalue();
            int cen2=(shell+jsh)->center();
            sq.al2=(shell+jsh)->alf_ptr();
            sq.co2=(shell+jsh)->cof_ptr();
            sq.b=(center+cen2)->r_vec();
            sq.ab2=dist_sqr(sq.a,sq.b);
            nls2=aux.number_of_lstates(sq.lv2);
            bool switch12=sq.lv1<sq.lv2;
            if (switch12)
            {
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
            for (int ksh=0;ksh<=ish;++ksh,kr0+=nls3)
            {
                sq.npr3=(shell+ksh)->number_of_prims();
                sq.lv3=(shell+ksh)->Lvalue();
                int cen3=(shell+ksh)->center();
                sq.al3=(shell+ksh)->alf_ptr();
                sq.co3=(shell+ksh)->cof_ptr();
                sq.c=(center+cen3)->r_vec();
                nls3=aux.number_of_lstates(sq.lv3);
                int lr0=0;
                int nls4=0;
                for (int lsh=0;lsh<=ksh;++lsh,lr0+=nls4)
                {
                    sq.npr4=(shell+lsh)->number_of_prims();
                    sq.lv4=(shell+lsh)->Lvalue();
                    int cen4=(shell+lsh)->center();
                    sq.al4=(shell+lsh)->alf_ptr();
                    sq.co4=(shell+lsh)->cof_ptr();
                    sq.d=(center+cen4)->r_vec();
                    sq.cd2=dist_sqr(sq.c,sq.d);
                    nls4=aux.number_of_lstates(sq.lv4);
#ifdef _MPI_API_
                    ++pknt;
                    pknt%=psize;
                    if (pknt!=rank) continue;
                    pknt=rank;
#endif
                    bool switch34=sq.lv3<sq.lv4;
                    if (switch34)
                    {
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
                    for (register int ils=0;ils<nls1;++ils,++ir)
                    {
                        int ijr=ir*(ir+1)/2+jr0;
                        int jr=jr0;
                        int jend=nls2;
                        if (ish==jsh) jend=ils+1;
                        for (register int jls=0;jls<jend;++ijr,++jr,++jls)
                        {
                            int kr=kr0;
                            int kend=nls3;
                            if (ish==ksh) kend=ils+1;
                            for (register int kls=0;kls<kend;++kls,++kr)
                            {
                                int klr=kr*(kr+1)/2+lr0;
                                int lr=lr0;
                                int lend=nls4;
                                if (ksh==lsh) lend=kls+1;
                                for (register int lls=0;lls<lend;++lr,++klr,++lls)
                                {
                                    if (klr>ijr) break;
                                    (sints+knt)->val=0.0;
                                    (sints+knt)->i=(unsigned short)ir;
                                    (sints+knt)->j=(unsigned short)jr;
                                    (sints+knt)->k=(unsigned short)kr;
                                    (sints+knt)->l=(unsigned short)lr;
                                    unsigned short l12=(ils<<4)+jls;
                                    if (switch12) l12=(jls<<4)+ils;
                                    unsigned short l34=(kls<<4)+lls;
                                    if (switch34) l34=(lls<<4)+kls;
                                    sq.lstates[knt]=(l12<<8)+l34;
                                    ++knt;
                                }
                            }
                        }
                    }
                    if (!knt) continue;
                    sq.len=knt;
                    calc_two_electron_ints(sq,aux,rys,sints,Gx,Gy,Gz);
                    for (register int kc=0;kc<knt;++kc)
                    {
                        if (fabs((sints+kc)->val)>threshold)
                        {
                            fwrite(sints+kc,sizeof(Sints),1,out);
                            ++fsize;
                        }
                    }
                    if (fsize>maxfsize)
                    {
                        fclose(out);
                        numints[nfiles]=fsize;
                        ++nfiles;
                        make_filename(nfiles,rank,basename,fname);
                        out=create_file(fname);
                        fsize=0;
                    }
                    if (switch34)
                    {
                        sq.npr3=sq.npr4;
                        sq.lv3=sq.lv4;
                        sq.al3=sq.al4;
                        sq.co3=sq.co4;
                        sq.c=sq.d;
                    }
                }
            }
            if (switch12)
            {
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
    if (fsize)
    {
        fclose(out);
        numints[nfiles]=fsize;
        ++nfiles;
    }
    delete [] sq.lstates;
    delete [] sints;
    delete_tensor3<double>(Gz,maxr,maxr);
    delete_tensor3<double>(Gy,maxr,maxr);
    delete_tensor3<double>(Gx,maxr,maxr);
}

inline void
TwoElectronInts::recalculate(const Basis& basis)
{
    const double threshold=1.e-14;
#ifdef _MPI_API_
    int pknt=0;
#endif    
    char fname[80];
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
    for (int i=0;i<start;++i,ir0+=nls1)
    {
        int lv1=(shell+i)->Lvalue();
        nls1=aux.number_of_lstates(lv1);
    }
    nfiles=0;
    Rys rys(maxl);
    int maxr=2*maxl+1;
    double ***Gx=new_tensor3<double>(maxr,maxr,maxr);
    double ***Gy=new_tensor3<double>(maxr,maxr,maxr);
    double ***Gz=new_tensor3<double>(maxr,maxr,maxr);
    Sints* sints=new Sints[ml4];
    ShellQuartet sq;
    int it;
    const double *dp;
    sq.lstates=new unsigned short[ml4];
    make_filename(nfiles,rank,basename,fname);
    FILE *out=create_file(fname);
    int maxfsize=MAXFILESIZE/sizeof(Sints);
    int fsize=0;
    clock_t s=clock();
    for (int ish=start;ish<nshell;++ish,ir0+=nls1)
    {
        sq.npr1=(shell+ish)->number_of_prims();
        sq.lv1=(shell+ish)->Lvalue();
        int cen1=(shell+ish)->center();
        sq.al1=(shell+ish)->alf_ptr();
        sq.co1=(shell+ish)->cof_ptr();
        sq.a=(center+cen1)->r_vec();
        nls1=aux.number_of_lstates(sq.lv1);
        int jr0=0;
        int nls2=0;
        for (int jsh=0;jsh<=ish;++jsh,jr0+=nls2)
        {
            sq.npr2=(shell+jsh)->number_of_prims();
            sq.lv2=(shell+jsh)->Lvalue();
            int cen2=(shell+jsh)->center();
            sq.al2=(shell+jsh)->alf_ptr();
            sq.co2=(shell+jsh)->cof_ptr();
            sq.b=(center+cen2)->r_vec();
            sq.ab2=dist_sqr(sq.a,sq.b);
            nls2=aux.number_of_lstates(sq.lv2);
            bool switch12=sq.lv1<sq.lv2;
            if (switch12)
            {
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
            for (int ksh=0;ksh<=ish;++ksh,kr0+=nls3)
            {
                sq.npr3=(shell+ksh)->number_of_prims();
                sq.lv3=(shell+ksh)->Lvalue();
                int cen3=(shell+ksh)->center();
                sq.al3=(shell+ksh)->alf_ptr();
                sq.co3=(shell+ksh)->cof_ptr();
                sq.c=(center+cen3)->r_vec();
                nls3=aux.number_of_lstates(sq.lv3);
                int lr0=0;
                int nls4=0;
                for (int lsh=0;lsh<=ksh;++lsh,lr0+=nls4)
                {
                    sq.npr4=(shell+lsh)->number_of_prims();
                    sq.lv4=(shell+lsh)->Lvalue();
                    int cen4=(shell+lsh)->center();
                    sq.al4=(shell+lsh)->alf_ptr();
                    sq.co4=(shell+lsh)->cof_ptr();
                    sq.d=(center+cen4)->r_vec();
                    sq.cd2=dist_sqr(sq.c,sq.d);
                    nls4=aux.number_of_lstates(sq.lv4);
#ifdef _MPI_API_
                    ++pknt;
                    pknt%=psize;
                    if (pknt!=rank) continue;
#endif
                    bool switch34=sq.lv3<sq.lv4;
                    if (switch34)
                    {
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
                    for (register int ils=0;ils<nls1;++ils,++ir)
                    {
                        int ijr=ir*(ir+1)/2+jr0;
                        int jr=jr0;
                        int jend=nls2;
                        if (ish==jsh) jend=ils+1;
                        for (register int jls=0;jls<jend;++ijr,++jr,++jls)
                        {
                            int kr=kr0;
                            int kend=nls3;
                            if (ish==ksh) kend=ils+1;
                            for (register int kls=0;kls<kend;++kls,++kr)
                            {
                                int klr=kr*(kr+1)/2+lr0;
                                int lr=lr0;
                                int lend=nls4;
                                if (ksh==lsh) lend=kls+1;
                                for (register int lls=0;lls<lend;++lr,++klr,++lls)
                                {
                                    if (klr>ijr) break;
                                    (sints+knt)->val=0.0;
                                    (sints+knt)->i=(unsigned short)ir;
                                    (sints+knt)->j=(unsigned short)jr;
                                    (sints+knt)->k=(unsigned short)kr;
                                    (sints+knt)->l=(unsigned short)lr;
                                    unsigned short l12=(ils<<4)+jls;
                                    if (switch12) l12=(jls<<4)+ils;
                                    unsigned short l34=(kls<<4)+lls;
                                    if (switch34) l34=(lls<<4)+kls;
                                    sq.lstates[knt]=(l12<<8)+l34;
                                    ++knt;
                                }
                            }
                        }
                    }
                    if (!knt) continue;
                    sq.len=knt;
                    calc_two_electron_ints(sq,aux,rys,sints,Gx,Gy,Gz);
                    for (register int kc=0;kc<knt;++kc)
                    {
                        if (fabs((sints+kc)->val)>threshold)
                        {
                            fwrite(sints+kc,sizeof(Sints),1,out);
                            ++fsize;
                        }
                    }
                    if (fsize>maxfsize)
                    {
                        fclose(out);
                        numints[nfiles]=fsize;
                        ++nfiles;
                        make_filename(nfiles,rank,basename,fname);
                        out=create_file(fname);
                        fsize=0;
                    }
                    if (switch34)
                    {
                        sq.npr3=sq.npr4;
                        sq.lv3=sq.lv4;
                        sq.al3=sq.al4;
                        sq.co3=sq.co4;
                        sq.c=sq.d;
                    }
                }
            }
            if (switch12)
            {
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
    if (fsize)
    {
        fclose(out);
        numints[nfiles]=fsize;
        ++nfiles;
    }
    delete [] sq.lstates;
    delete [] sints;
    delete_tensor3<double>(Gz,maxr,maxr);
    delete_tensor3<double>(Gy,maxr,maxr);
    delete_tensor3<double>(Gx,maxr,maxr);
}

inline void
TwoElectronInts::formGmatrix(const double* Pmat,double *Gmat) const
{
    const int BINSIZE=1000;
    Sints sints[BINSIZE];
    char fname[80];

    for (int ifile=0;ifile<nfiles;++ifile)
    {
        int nints=numints[ifile];
        int nbin=nints/BINSIZE;
        int nex=nints-BINSIZE*nbin;
        make_filename(ifile,rank,basename,fname);
        FILE *in=open_file(fname);
        if (nex)
        {
            fread(sints,sizeof(Sints),nex,in);
            for (register int ix=0;ix<nex;++ix)
            {
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
                if (j>k)
                {
                    jk=j*(j+1)/2+k;
                }else{
                    jk=k*(k+1)/2+j;
                }
                if (j>l)
                {
                    jl=j*(j+1)/2+l;
                }else{
                    jl=l*(l+1)/2+j;
                }
                double da=val*2.0*Pmat[ij];
                double db=val*2.0*Pmat[kl];
                double sjl=val*Pmat[ik];
                double sjk=val*Pmat[il];
                double sik=val*Pmat[jl];
                double sil=val*Pmat[jk];
                if (k!=l)
                {
                    db=db+db;
                    Gmat[ik]-=sik;
                    if (i!=j && j>=k) Gmat[jk]-=sjk;
                }
                Gmat[il]-=sil;
                Gmat[ij]+=db;
                if (i!=j && j>=l) Gmat[jl]-=sjl;
                if (ij!=kl) {
                    if (i!=j) da=da+da;
                    if (j<=k)
                    {
                        Gmat[jk]-=sjk;
                        if (i!=j && i<=k) Gmat[ik]-=sik;
                        if (k!=l && j<=l) Gmat[jl]-=sjl;
                    }
                    Gmat[kl]+=da;
                }
            }
        }
        for (int ibin=0;ibin<nbin;++ibin)
        {
            fread(sints,sizeof(Sints),BINSIZE,in);
            for (register int ix=0;ix<BINSIZE;++ix)
            {
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
                if (j>k)
                {
                    jk=j*(j+1)/2+k;
                }else{
                    jk=k*(k+1)/2+j;
                }
                if (j>l)
                {
                    jl=j*(j+1)/2+l;
                }else{
                    jl=l*(l+1)/2+j;
                }
                double da=val*2.0*Pmat[ij];
                double db=val*2.0*Pmat[kl];
                double sjl=val*Pmat[ik];
                double sjk=val*Pmat[il];
                double sik=val*Pmat[jl];
                double sil=val*Pmat[jk];
                if (k!=l)
                {
                    db=db+db;
                    Gmat[ik]-=sik;
                    if (i!=j && j>=k) Gmat[jk]-=sjk;
                }
                Gmat[il]-=sil;
                Gmat[ij]+=db;
                if (i!=j && j>=l) Gmat[jl]-=sjl;
                if (ij!=kl) {
                    if (i!=j) da=da+da;
                    if (j<=k)
                    {
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

inline void
TwoElectronInts::formGmatrix(const double* PmatA,const double *PmatB,
                             double *GmatA,double *GmatB) const
{
    const int BINSIZE=1000;
    Sints sints[BINSIZE];
    char fname[80];

    for (int ifile=0;ifile<nfiles;++ifile)
    {
        int nints=numints[ifile];
        int nbin=nints/BINSIZE;
        int nex=nints-BINSIZE*nbin;
        make_filename(ifile,rank,basename,fname);
        FILE *in=open_file(fname);
        if (nex)
        {
            fread(sints,sizeof(Sints),nex,in);
            for (register int ix=0;ix<nex;++ix)
            {
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
                if (j>k)
                {
                    jk=j*(j+1)/2+k;
                }else{
                    jk=k*(k+1)/2+j;
                }
                if (j>l)
                {
                    jl=j*(j+1)/2+l;
                }else{
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
                if (k!=l)
                {
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
                if (i!=j && j>=l)
                {
                    GmatA[jl]-=sjlA;
                    GmatB[jl]-=sjlB;
                }
                if (ij!=kl) {
                    if (i!=j) da=da+da;
                    if (j<=k)
                    {
                        GmatA[jk]-=sjkA;
                        GmatB[jk]-=sjkB;
                        if (i!=j && i<=k)
                        {
                            GmatA[ik]-=sikA;
                            GmatB[ik]-=sikB;
                        }
                        if (k!=l && j<=l)
                        {
                            GmatA[jl]-=sjlA;
                            GmatB[jl]-=sjlB;
                        }
                    }
                    GmatA[kl]+=da;
                    GmatB[kl]+=da;
                }
            }
        }
        for (int ibin=0;ibin<nbin;++ibin)
        {
            fread(sints,sizeof(Sints),BINSIZE,in);
            for (register int ix=0;ix<BINSIZE;++ix)
            {
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
                if (j>k)
                {
                    jk=j*(j+1)/2+k;
                }else{
                    jk=k*(k+1)/2+j;
                }
                if (j>l)
                {
                    jl=j*(j+1)/2+l;
                }else{
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
                if (k!=l)
                {
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
                if (i!=j && j>=l)
                {
                    GmatA[jl]-=sjlA;
                    GmatB[jl]-=sjlB;
                }
                if (ij!=kl) {
                    if (i!=j) da=da+da;
                    if (j<=k)
                    {
                        GmatA[jk]-=sjkA;
                        GmatB[jk]-=sjkB;
                        if (i!=j && i<=k)
                        {
                            GmatA[ik]-=sikA;
                            GmatB[ik]-=sikB;
                        }
                        if (k!=l && j<=l)
                        {
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

#endif
