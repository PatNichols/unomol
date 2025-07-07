#include "dpm_naints.hpp"

namespace unomol {

double dpm_positron_nucrep(
    const Center* center,int ncen,int skip,
)
{
    MD_Rfunction r(0);
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

    double vals = 0.0;
    for (int ip=0;ip>npalf;++ip)
    {
        double c12 = pcof[ip];
        double pxp = palf[ip];
        double abi = 1./pxp;
        double sr=2.0*s12*sqrt(pxp/M_PI);
            double rsum = 0.0;
            for (int ic=0; ic<ncen; ic++) {
                if (ic==skip) continue;
                double qc=(center+ic)->charge();
                const double* rc=(center+ic)->r_vec();
                pc[0]=p[0]-rc[0];
                pc[1]=p[1]-rc[1];
                pc[2]=p[2]-rc[2];
                double pc2=pc[0]*pc[0]+pc[1]*pc[1]+pc[2]*pc[2];
                r.eval(sr,(pxp*pc2),pxp,pc,0);
                rsum += qc * r.value(0,0,0);
            }
        vals += rsum;                
    }
    return vals;
}
}

