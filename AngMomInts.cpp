#include "AngMomInts.hpp"
#include <iostream>
namespace unomol {

struct AngMomVec {
    double x;
    double y;
    double z;
};

void calc_angmom_ints(
    const ShellPairData & sp,
    AngMomVec * aints,
    const AuxFunctions& aux,
    MD_Dfunction& dx,
    MD_Dfunction& dy,
    MD_Dfunction& dz,
    const double *cen)
{
    int ix,iy,iz;
    const double piterm=5.5683279968317079;
    double p[3];
    double pc[3];
    int lvt=sp.lv1+sp.lv2;
    for (int ip=0; ip<sp.npr1; ip++) {
        double axp=sp.al1[ip];
        double c1=sp.co1[ip];
        for (int jp=0; jp<sp.npr2; jp++) {
            double c2 = sp.co2[jp];
            double c12=c1*c2;
            double bxp=sp.al2[jp];
            double pxp=axp+bxp;
            double abi=1.0/pxp;
            double s12=exp(-axp*bxp*abi*sp.ab2);
            p[0]=(axp*sp.a[0]+bxp*sp.b[0])*abi;
            p[1]=(axp*sp.a[1]+bxp*sp.b[1])*abi;
            p[2]=(axp*sp.a[2]+bxp*sp.b[2])*abi;
            pc[0] = p[0] - cen[0];
            pc[1] = p[1] - cen[1];
            pc[2] = p[2] - cen[2];
            abi=abi*0.5;
            dx.eval(abi,p[0]-sp.a[0],p[0]-sp.b[0],sp.lv1+1,sp.lv2+1);
            dy.eval(abi,p[1]-sp.a[1],p[1]-sp.b[1],sp.lv1+1,sp.lv2+1);
            dz.eval(abi,p[2]-sp.a[2],p[2]-sp.b[2],sp.lv1+1,sp.lv2+1);
            for (int kc=0; kc<sp.len; kc++) {
                unsigned short key=sp.lstates[kc];
                int ls2=key&0xF;
                int ls1=key>>4;
                int lx1=aux.Lxyz(sp.lv1,ls1,0);
                int ly1=aux.Lxyz(sp.lv1,ls1,1);
                int lz1=aux.Lxyz(sp.lv1,ls1,2);
                int lx2=aux.Lxyz(sp.lv2,ls2,0);
                int ly2=aux.Lxyz(sp.lv2,ls2,1);
                int lz2=aux.Lxyz(sp.lv2,ls2,2);
                double nfct=aux.normalization_factor(sp.lv1,ls1)*
                            aux.normalization_factor(sp.lv2,ls2)*c12*s12;                
                double s0x = dx(lx1,lx2,0);
                double s0y = dy(ly1,ly2,0);
                double s0z = dz(lz1,lz2,0);
                double s1x = dx(lx1,lx2,1) + pc[0] * s0x;
                double s1y = dy(ly1,ly2,1) + pc[1] * s0y;
                double s1z = dz(lz1,lz2,1) + pc[2] * s0z;
                double twoc = -2. * c2;
                double d1x = twoc * dx(lx1,lx2+1,0);
                double d1y = twoc * dy(ly1,ly2+1,0);
                double d1z = twoc * dz(lz1,lz2+1,0);
                if ( lx2 ) d1x += lx2 * dx(lx1,lx2-1,0);
                if ( ly2 ) d1y += ly2 * dy(ly1,ly2-1,0);
                if ( lz2 ) d1z += lz2 * dz(lz1,lz2-1,0);
                double ax = -s0x * ( s1y * d1z - s1z * d1y);
                double ay = -s0y * ( s1z * d1x - s1x * d1z);
                double az = -s0z * ( s1x * d1y - s1y * d1x);
                aints[kc].x += ax * nfct;
                aints[kc].y += ay * nfct;
                aints[kc].z += az * nfct;
            }
        }
    }
};

void AngMomInts(const Basis& bas) {
    int maxl=bas.maxLvalue();
    MD_Dfunction dx(maxl+2);
    MD_Dfunction dy(maxl+2);
    MD_Dfunction dz(maxl+2);
    double cen[3];
    int nshell=bas.number_of_shells();
    const Shell* shell=bas.shell_ptr();
    const Center* center=bas.center_ptr();
    const AuxFunctions& aux(*bas.auxfun_ptr());
    int ml2=aux.number_of_lstates( maxl );
    ml2=ml2*ml2;
    AngMomVec *Lvals=new AngMomVec[ml2];
    ShellPairData sp;
    sp.lstates=new unsigned short[ml2];
    int * ostates=new int[ml2];
    cen[0] = cen[1] = cen[2] = 0.0;
    std::ofstream out("angmom.out");
    for (int ishell=0; ishell<nshell; ++ishell) {
        int ir0 = bas.offset(ishell);
        sp.npr1=(shell+ishell)->number_of_prims();
        sp.lv1=(shell+ishell)->Lvalue();
        int cn1=(shell+ishell)->center();
        sp.al1=(shell+ishell)->alf_ptr();
        sp.co1=(shell+ishell)->cof_ptr();
        sp.a=(center+cn1)->r_vec();
        int nls1=aux.number_of_lstates(sp.lv1);
        for (int jshell=0; jshell<=ishell; ++jshell) {
            int jr0 = bas.offset(jshell);
            sp.npr2=(shell+jshell)->number_of_prims();
            sp.lv2=(shell+jshell)->Lvalue();
            int cn2=(shell+jshell)->center();
            sp.al2=(shell+jshell)->alf_ptr();
            sp.co2=(shell+jshell)->cof_ptr();
            sp.b=(center+cn2)->r_vec();
            sp.ab2=dist_sqr(sp.a,sp.b);
            int nls2=aux.number_of_lstates(sp.lv2);
            int knt=0;
            for (int ils=0; ils<nls1; ++ils) {
                int ir = ir0 + ils;
                int iir = ir * ( ir + 1 ) / 2;
                for (int jls=0; jls<nls2; ++jls) {
                    int jr = jr0 + jls;
                    if (jr>ir) break;
                    sp.lstates[knt]=(unsigned short)((ils<<4)+jls);
                    ostates[knt]= iir + jr;
                    Lvals[knt].x=0.0;
                    Lvals[knt].y=0.0;
                    Lvals[knt].z=0.0;
                    ++knt;
                }
            }
            if (!knt) continue;
            sp.len=knt;
            calc_angmom_ints(sp,Lvals,aux,dx,dy,dz,cen);
            for (int kc=0; kc<knt; kc++) {
                int ijr=ostates[kc];
                out.write((const char*)&ijr,sizeof(int));
                out.write((const char*)(Lvals+kc),sizeof(AngMomVec));

            }
        }
    }
    out.close();
    delete [] ostates;
    delete [] sp.lstates;
    delete [] Lvals;
};

#define __check_input

void calculateAngMomentum(const double * pmat,size_t no2)
{
    int ij;
    double lx = 0;
    double ly = 0;
    double lz = 0;
    std::ifstream in("angmom.out");
//    __check_input(in,"angmom.out");
    for (size_t i=0;i<no2;++i)
    {
        AngMomVec lvec;
        in.read((char*)&ij,sizeof(int));
        in.read((char*)&lvec,sizeof(AngMomVec));
        lx += pmat[i] * lvec.x;
        ly += pmat[i] * lvec.y;
        lz += pmat[i] * lvec.z;
    }
    in.close();
    std::cout << "Angular Momentum\n";
    std::cout << "Lx = " << lx << "\n";
    std::cout << "Ly = " << ly << "\n";
    std::cout << "Lz = " << lz << "\n";    
}

void calculateAngMomentum(const double * pmatA, const double * pmatB, size_t no2)
{
    int ij;
    AngMomVec lvec;
    double lx = 0;
    double ly = 0;
    double lz = 0;
    std::ifstream in("angmom.out");
//    __check_input(in,"angmom.out");
    for (size_t i=0;i<no2;++i)
    {
        in.read((char*)&ij,sizeof(int));
        in.read((char*)&lvec,sizeof(AngMomVec));
        double ptot = pmatA[i] + pmatB[i];
        lx += ptot * lvec.x;
        ly += ptot * lvec.y;
        lz += ptot * lvec.z;
    }
    in.close();
    std::cout << "Angular Momentum\n";
    std::cout << "Lx = " << lx << "\n";
    std::cout << "Ly = " << ly << "\n";
    std::cout << "Lz = " << lz << "\n";    
}

} // end namespace
