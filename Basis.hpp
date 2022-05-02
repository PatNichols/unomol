#ifndef UNOMOL_BASIS_hpp_
#define UNOMOL_BASIS_hpp_
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cfloat>
#include <iomanip>
using namespace std;
#include "Util.hpp"
#include "AuxFunctions.hpp"
#include "SymmPack.hpp"

namespace unomol {

class Shell {
  private:
    int npr,lsh,cen;
    double* al;
    double* co;
  public:
    Shell():npr(0),lsh(0),cen(0),al(0),co(0) {};
    Shell(const Shell& rs):npr(rs.npr),lsh(rs.lsh),cen(rs.cen),
        al(new double[npr]),co(new double[npr]) {
        for (int i=0; i<npr; ++i) {
            al[i]=rs.al[i];
            co[i]=rs.co[i];
        }
    }
    constexpr double alf(int i) const noexcept {
        return al[i];
    };
    constexpr double cof(int i) const noexcept {
        return co[i];
    };
    constexpr const double* alf_ptr() const noexcept {
        return al;
    };
    constexpr const double* cof_ptr() const noexcept {
        return co;
    };
    constexpr int number_of_prims() const noexcept {
        return npr;
    };
    constexpr int Lvalue() const noexcept {
        return lsh;
    };
    constexpr int center() const noexcept {
        return cen;
    };
    constexpr void setCenter(int icen) {
        cen=icen;
    };

    void normalize() noexcept {
        const double twofact= 2.8284271247461903;
        const double piterm=5.568327996831707;
        double lpow=1.5+lsh;
        double sum=0.0;
        for (int i=0; i<npr; i++) {
            double a1=al[i];
            double c1=co[i];
            for (int j=0; j<npr; j++) {
                double a2=al[j];
                double c2=co[j];
                sum+=c1*c2*pow(sqrt(a1*a2)/(a1+a2),lpow);
            }
        }
        sum*=twofact;
        sum=1.0/sqrt(sum);
        for (int i=0; i<npr; i++) {
            co[i]=co[i]*sum*sqrt(pow(2*al[i],lpow)/piterm);
        }
    };

    friend std::ostream& operator << ( std::ostream& os,
                                       const Shell& sh) {
        os<<std::setw(3)<<sh.npr<<sh.lsh<<sh.cen<<endl;
        for (int i=0; i<sh.npr; ++i) {
            os.setf(std::ios::showpoint);
            os<< std::setw(20) << std::setprecision(10)
              << sh.al[i] << " " << sh.co[i]<<endl;
        }
        return os;
    }

    friend std::istream& operator >> ( std::istream& is,
                                       Shell& sh) {
        is>>sh.npr;
        is>>sh.lsh;
        is>>sh.cen;
        if (sh.co) delete [] sh.co;
        if (sh.al) delete [] sh.al;
        sh.al=new double[sh.npr];
        sh.co=new double[sh.npr];
        for (int i=0; i<sh.npr; ++i) {
            is>>sh.al[i];
            is>>sh.co[i];
        }
        return is;
    }
};

class Center {
  private:
    double chg;
    double r[3];
  public:
    Center() {
        chg=0.0;
        r[2]=r[1]=r[0]=0.0;
    }

    constexpr double charge() const noexcept {
        return chg;
    };

    constexpr const double* r_vec() const noexcept {
        return r;
    };

    constexpr double position(int i) const noexcept {
        return r[i];
    };

    constexpr void setPosition(double rx,double ry,double rz) noexcept {
        r[0]=rx;
        r[1]=ry;
        r[2]=rz;
    }

    constexpr void setCharge(double q) noexcept {
        chg=q;
    }

    friend ostream& operator << (std::ostream& os,const Center& c) {
        os.setf(std::ios::showpoint);
        os<<std::setprecision(10)<<std::setw(15)
          <<c.chg<<" "
          <<c.r[0]<<" "
          <<c.r[1]<<" "
          <<c.r[2]<<" "
          <<endl;
        return os;
    }

    friend istream& operator >> (std::istream& is,Center& c) {
        is>>c.chg;
        is>>c.r[0];
        is>>c.r[1];
        is>>c.r[2];
        return is;
    }
};

class Basis {
  private:
    double eps,ffield;
    double *trans;
    int nshell,ncen,norb,maxl,nelec,maxits;
    int tnshell,tncen,tnorb,skipcen,tnelec;
    int scf_flag[3];
    int int_flag[3];
    int prt_flag[4];
    Center* centers;
    Shell* shells;
    AuxFunctions* aux;
  public:
    Basis(const std::string& infile=string("patin.dat")) {
        int xsh,xno,xmaxl;
        std::ifstream ain;
        std::ifstream in(infile.c_str());
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
        if (int_flag[0]==1) {
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
        for (int i=0; i<ncen; ++i) in>>centers[i];
        for (int i=0; i<nshell; ++i) in>>shells[i];
        ffield=5.e-3;
        in.close();
        if (int_flag[0]==1) {
            (centers+ncen)->setCharge(1.0000000000000000000);
            (centers+ncen)->setPosition(0.0,0.0,0.0);
            for (int i=nshell; i<tnshell; ++i) {
                ain>>shells[i];
                (shells+i)->setCenter(ncen);
            }
            ain.close();
        }
        for (int i=0; i<tnshell; i++) (shells+i)->normalize();
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

    ~Basis() {
        delete [] centers;
        delete [] shells;
    }

    constexpr int maxLvalue() const noexcept {
        return maxl;
    }
    constexpr int number_of_electrons() const noexcept {
        return nelec;
    }
    constexpr int number_of_shells() const noexcept {
        return nshell;
    }
    constexpr int number_of_centers() const noexcept {
        return ncen;
    }
    constexpr int number_of_orbitals() const noexcept {
        return norb;
    }
    constexpr int total_number_of_shells() const noexcept {
        return tnshell;
    }
    constexpr int total_number_of_centers() const noexcept {
        return tncen;
    }
    constexpr int total_number_of_orbitals() const noexcept {
        return tnorb;
    }
    constexpr int maximum_iterations() const noexcept {
        return maxits;
    }
    constexpr int scf_flags(int i) const noexcept {
        return scf_flag[i];
    }
    constexpr int prt_flags(int i) const noexcept {
        return prt_flag[i];
    }
    constexpr int int_flags(int i) const noexcept {
        return int_flag[i];
    }
    constexpr const Shell* shell_ptr() const noexcept {
        return shells;
    }
    constexpr const Center* center_ptr() const noexcept {
        return centers;
    }
    constexpr const AuxFunctions* auxfun_ptr() const noexcept {
        return aux;
    }
    constexpr double scf_eps() const noexcept {
        return eps;
    }
    constexpr int skip_center() const noexcept {
        return skipcen;
    }
    constexpr void dpm_augment() noexcept {
        nshell=tnshell;
        norb=tnorb;
        ++ncen;
    }
    constexpr void SetCenterPosition(double x,double y,double z,int n) noexcept {
        (centers+n)->setPosition(x,y,z);
    }
    constexpr double FiniteFieldValue() const noexcept {
        return ffield;
    }
    constexpr const double *transformation_matrix() const noexcept {
        return trans;
    }
    void find_com(double *trans) noexcept {
        int i,icen;
        const double *ri;
        double q[6];
        double tmp[3];
        double qvals[3];
        double sumx,sumy,sumz,sumc;
        double qi,cx,cy,cz,rx,ry,rz,x,y,z;

        sumx=sumy=sumz=sumc=0.0;
        for (i=0; i<ncen; ++i) {
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
        for (icen=0; icen<ncen; ++icen) {
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
        for (i=0; i<6; ++i) q[i]/=sumc;
        SymmPack::rsp(3,q,trans,qvals,tmp);
    }
};

}
#endif

