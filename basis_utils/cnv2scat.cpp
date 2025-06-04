#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <cmath>

int main()
{
    std::ifstream in("pos.grid.dat");
    int npts;
    in >> npts;
    double * r = new double[npts];
    double * vstat =  new double[npts];
    double * vpol =  new double[npts];
    double * alf = new double[npts];

    for (int i=0;i<npts;++i) {
        double x,y,z;
        in >> x;
        in >> y;
        in >> z;
        r[i] = sqrt(x*x+y*y+z*z);
    }
    in.close();
    in.open("spol.out");
    for (int i=0;i<npts;++i) {
        double eg,ef,ei,de,vs;
        int ix;
        in >> ix;
        in >> eg;
        in >> ei;
        in >> ef;
        in >> vs;
        in >> de;
        vstat[i] = vs;
    }
    in.close();
    in.open("vpol.out");
    for (int i=0;i<npts;++i) {
        double px,py,pz;
        double vp,alfa;
        in >> px;
        in >> py;
        in >> pz;
        in >> vp;
        in >> alfa;
        vpol[i] = vp;
        alf[i] = alfa;
    }
    in.close();

    std::ofstream out("vsph_scat.dat");
    out << " " << npts << " ";
    out << std::setw(16) << std::setprecision(8) << std::fixed << r[npts-1] << " ";
    out << std::setw(16) << std::setprecision(8) << std::fixed << alf[npts-1] << "\n";
    for (int i=0;i<npts;++i)
    {
        out << std::setw(16) << std::setprecision(8) << std::fixed;
        out << r[i] << " ";
        out << std::setw(20) << std::setprecision(12) << std::scientific;         
        out << (vstat[i] + vpol[i]) << "\n";    
    }
    out.close();

    std::ofstream out2("vsp.out");
    out2 << " " << npts << " ";
    out2 << std::setw(16) << std::setprecision(8) << std::fixed << r[npts-1] << " ";
    out2 << std::setw(16) << std::setprecision(8) << std::fixed << alf[npts-1] << "\n";
    for (int i=0;i<npts;++i)
    {
        out2 << std::setw(16) << std::setprecision(8) << std::fixed;
        out2 << r[i] << " ";
        out2 << std::setw(20) << std::setprecision(12) << std::scientific;         
        out2 << (vstat[i]) << " ";    
        out2 << std::setw(20) << std::setprecision(12) << std::scientific;         
        out2 << (vpol[i]) << "\n";    
    }
    out2.close();


    delete [] alf;
    delete [] vpol;
    delete [] vstat;
    delete [] r;
}