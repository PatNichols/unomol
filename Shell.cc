#ifndef _SHELL_CC_
#define _SHELL_CC_
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
using namespace std;

class Shell
{
private:
    int npr,lsh,cen;
    double* al;
    double* co;
public:
    Shell():npr(0),lsh(0),cen(0),al(0),co(0) {};
    Shell(const Shell& rs):npr(rs.npr),lsh(rs.lsh),cen(rs.cen),
            al(new double[npr]),co(new double[npr])
    {
        for (register int i=0;i<npr;++i)
        {
            al[i]=rs.al[i];
            co[i]=rs.co[i];
        }
    }
    double alf(int i) const {  return al[i];};
    double cof(int i) const {  return co[i];};
    const double* alf_ptr() const { return al;};
    const double* cof_ptr() const { return co;};
    int number_of_prims() const { return npr;};
    int Lvalue() const { return lsh;};
    int center() const { return cen;};
    istream& init_from_stream(istream& in);
    ostream& write_to_stream(ostream& os) const;
    void setCenter(int icen) { cen=icen; };
    inline void normalize()
    {
        const double twofact= 2.8284271247461903;
        const double piterm=5.568327996831707;
        double lpow=1.5+lsh;
        double sum=0.0;
        for (int i=0;i<npr;i++)
        {
            double a1=al[i];
            double c1=co[i];
            for (int j=0;j<npr;j++)
            {
                double a2=al[j];
                double c2=co[j];
                sum+=c1*c2*pow(sqrt(a1*a2)/(a1+a2),lpow);
            }
        }
        sum*=twofact;
        sum=1.0/sqrt(sum);
        for (int i=0;i<npr;i++)
        {
            co[i]=co[i]*sum*sqrt(pow(2*al[i],lpow)/piterm);
        }
    };
};

ostream& Shell::write_to_stream(ostream& os) const
{
    os<<setw(3)<<npr<<lsh<<cen<<endl;
    for (register int i=0;i<npr;++i)
    {
        os.setf(ios::showpoint);
        os<<setw(20)<<setprecision(10)
        <<al[i]<<" "<<co[i]<<endl;
    }
    return os;
}


istream& Shell::init_from_stream(istream& is)
{
    is>>npr;
    is>>lsh;
    is>>cen;
    if (co) delete [] co;
    if (al) delete [] al;
    this->al=new double[npr];
    this->co=new double[npr];
    for (int i=0;i<npr;++i) {
        is>>al[i];
        is>>co[i];
    }
    return is;
}


ostream& operator<<(ostream& os,const Shell& sh)
{
    return sh.write_to_stream(os);
}


istream& operator>>(istream& is,Shell& sh)
{
    return sh.init_from_stream(is);
}
#endif
