#ifndef _BASIS_CC_
#define _BASIS_CC_
#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>
using namespace std;
#include "Util.cc"
#include "Shell.cc"
#include "Center.cc"
#include "AuxFunctions.cc"
#include "SymmPack.cc"

class Basis
{
private:
    double eps,ffield;
    double *trans;
    int nshell,ncen,norb,maxl,nelec,maxits;
    int tnshell,tncen,tnorb,skipcen,tnelec,dpmgrid;
    int scf_flag[3];
    int int_flag[3];
    int prt_flag[4];
    Center* centers;
    Shell* shells;
    AuxFunctions* aux;
    void find_com(double*);
public:
    Basis(const std::string& infile=string("patin.dat"));
    ~Basis()
    {
        delete [] centers;
        delete [] shells;
    }
    int maxLvalue() const { return maxl;};
    int number_of_electrons() const { return nelec;};
    int number_of_shells() const { return nshell;};
    int number_of_centers() const { return ncen;};
    int number_of_orbitals() const { return norb;};
    int total_number_of_shells() const { return tnshell;};
    int total_number_of_centers() const { return tncen;};
    int total_number_of_orbitals() const { return tnorb;};
    int maximum_iterations() const { return maxits;};
    int scf_flags(int i) const { return scf_flag[i];};
    int prt_flags(int i) const { return prt_flag[i];};
    int int_flags(int i) const { return int_flag[i];};
    const Shell* shell_ptr() const { return shells;};
    const Center* center_ptr() const { return centers;};
    const AuxFunctions* auxfun_ptr() const { return aux;};
    double scf_eps() const { return eps;};
    int skip_center() const { return skipcen;};
    void dpm_augment();
    void SetCenterPosition(double x,double y,double z,int n)
    {
        (centers+n)->setPosition(x,y,z);
    }
    int number_of_dpm_points() { return dpmgrid;};
    double FiniteFieldValue() { return ffield;};
    const double *transformation_matrix() const { return trans;};
};

inline Basis::Basis(const std::string& infile)
{
    int xsh,xno,xmaxl;
    ifstream ain;
    ifstream in(infile.c_str());
    if (!in) {
        std::string errmsg="could not open file "+infile;
        fatal_error(errmsg.c_str());
    }
    in>>nshell;
    in>>norb;
    in>>ncen;
    in>>maxl;
    in>>nelec;
    in>>maxits;
    in>>eps;
    in>>int_flag[0];
    in>>int_flag[1];
    in>>scf_flag[0];
    in>>scf_flag[1];
    in>>scf_flag[2];
    in>>prt_flag[0];
    in>>prt_flag[1];
    in>>prt_flag[2];
    tnshell=nshell;
    tncen=ncen;
    tnorb=norb;
    tnelec=nelec;
    if (int_flag[0]==1)
    {
        ain.open("posin.dat");
        if (!ain) {
            fatal_error("could not open posin.dat");
        }
        ain>>xsh;
        ain>>xno;
        ain>>xmaxl;
        tnshell+=xsh;
        tnorb+=xno;
        ++tncen;
        if (xmaxl>maxl) maxl=xmaxl;
    }
    if (maxl>4) fatal_error("Angular Momentum is too large for present program\n");
    centers=new Center[tncen];
    shells=new Shell[tnshell];
    for (int i=0;i<ncen;++i) in>>centers[i];
    for (int i=0;i<nshell;++i) in>>shells[i];
    ffield=5.e-3;
    in.close();
    if (int_flag[0]==1)
    {
        (centers+ncen)->setCharge(1.0000000000000000000);
        (centers+ncen)->setPosition(0.0,0.0,0.0);
        for (int i=nshell;i<tnshell;++i) {
            ain>>shells[i];
            (shells+i)->setCenter(ncen);
        }
        ain>>dpmgrid;
        ain.close();
    }
    for (int i=0;i<tnshell;i++) (shells+i)->normalize();
    skipcen=ncen;
    aux=new AuxFunctions(maxl);
    trans=new double[9];
    find_com(trans);
    double xeps=DBL_EPSILON*norb*norb;
    if (eps<xeps) {
	fprintf(stderr,"SCF convergence of %15.6le is way too low!\n",eps); 
	eps=xeps;
	fprintf(stderr,"SCF Convergence set at %15.6le\n",eps);
    }
}

void Basis::dpm_augment()
{
    nshell=tnshell;
    norb=tnorb;
    ++ncen;
}

inline void Basis::find_com(double *trans)
{
    int i,icen;
    const double *ri;
    double q[6];
    double tmp[3];
    double qvals[3];
    double sumx,sumy,sumz,sumc;
    double qi,cx,cy,cz,rx,ry,rz,x,y,z;

    sumx=sumy=sumz=sumc=0.0;
    for (i=0;i<ncen;++i) {
        qi=(centers+i)->charge();
        ri=(centers+i)->r_vec();
        sumx+=ri[0]*qi;
        sumy+=ri[1]*qi;
        sumz+=ri[2]*qi;
        sumc+=qi;
    }
    cx=sumx/sumc;
    cy=sumy/sumc;
    cz=sumz/sumc;
    q[5]=q[4]=q[3]=q[2]=q[1]=q[0]=0.0;
    for (icen=0;icen<ncen;++icen)
    {
        qi=(centers+icen)->charge();
        ri=(centers+icen)->r_vec();
        x=ri[0]-cx;
        y=ri[1]-cy;
        z=ri[2]-cz;
        q[0]+= (-2.0*x*x+y*y+z*z)*qi;
        q[1]+= (-3.0*x*y)*qi;
        q[2]+= (-2.0*y*y+x*x+z*z)*qi;
        q[3]+= (-3.0*x*z)*qi;
        q[4]+= (-3.0*y*z)*qi;
        q[5]+= (-2.0*z*z+x*x+y*y)*qi;
    }
    for (i=0;i<6;++i) q[i]/=sumc;
    SymmPack::rsp(3,q,trans,qvals,tmp);
}


#endif

