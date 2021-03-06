#ifndef _CENTER_CC_
#define _CENTER_CC_
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

class Center
{
private:
    double chg;
    double r[3];

public:
    Center()
    {
        chg=0.0;
        r[2]=r[1]=r[0]=0.0;
    }

    double charge() const { return chg; };

    const double* r_vec() const { return r;};

    double position(int i) const { return r[i];};

    void setPosition(double rx,double ry,double rz)
    {
        r[0]=rx;
        r[1]=ry;
        r[2]=rz;
    }

    void setCharge(double q)
    {
        chg=q;
    }

    istream& init_from_stream(istream& in);
    ostream& write_to_stream(ostream& os) const;
};

ostream& Center::write_to_stream(ostream& os) const
{
    os.setf(ios::showpoint);
    os<<setprecision(10)<<setw(15)
    <<chg<<" "
    <<r[0]<<" "
    <<r[1]<<" "
    <<r[2]<<" "
    <<endl;
    return os;
}


istream& Center::init_from_stream(istream& is)
{
    is>>chg;
    is>>r[0];
    is>>r[1];
    is>>r[2];
    return is;
}


ostream& operator<<(ostream& os,const Center& c)
{
    return c.write_to_stream(os);
}


istream& operator>>(istream& is,Center& c)
{
    return c.init_from_stream(is);
}
#endif
