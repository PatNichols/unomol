#ifndef UNOMOL_STRUCTS_HPP
#define UNOMOL_STRUCTS_HPP

namespace unomol {

namespace utils {
inline int number_of_lstates(int l) { return ((l+1)*(l+2))/2;}
}

struct MomInts {
    double dx,dy,dz,qxx,qxy,qxz,qyy,qyz,qzz;
    unsigned int ijr;
};

struct Moments {
    double dx,dy,dz,qxx,qxy,qxz,qyy,qyz,qzz;
};


struct ShellPair
{
    double * norms;
    unsigned int * lstates;
    unsigned int * orbs;
    //
    const double * a;
    const double * b;
    const double * al1;
    const double * co1;
    const double * al2;
    const double * co2;
    double ab2;
    int npr1,lv1;
    int npr2,lv2;
    unsigned int len;
    
    ShellPair() = delete;
    
    ShellPair(int maxl)
    {
        int maxlst = utils::number_of_lstates(maxl);
        int maxints = maxlst * maxlst;
        orbs = new unsigned int[maxints];
        norms = new double[maxints];
        lstates = new unsigned int[maxints];
        memset(lstates,0x0,sizeof(unsigned int)*maxints);
        len = 0;
    }
    
    ~ShellPair()
    {
        delete [] lstates;
        delete [] norms;
        delete [] orbs;
    }
    
    
    void assign1(const Shell& sh, const Center * center)
    {
        npr1 = sh.number_of_prims();
        lv1 = sh.Lvalue();
        int cn1 = sh.center();
        al1 = sh.alf_ptr();
        co1 = sh.cof_ptr();
        a = (center+cn1)->r_vec();
    }
    void assign2(const Shell& sh, const Center * center)
    {
        npr2 = sh.number_of_prims();
        lv2 = sh.Lvalue();
        int cn2 = sh.center();
        al2 = sh.alf_ptr();
        co2 = sh.cof_ptr();
        b = (center+cn2)->r_vec();
        ab2 = dist_sqr(a,b);
    }

    unsigned int precalc(const AuxFunctions& aux,
        unsigned int off1, unsigned int off2)
    {
        int nls1 = aux.number_of_lstates(lv1);
        int nls2 = aux.number_of_lstates(lv2);
        len = 0;
        for (unsigned int ls1=0;ls1<nls1;++ls1)
        {
            unsigned int ir1 = off1 + ls1;
            unsigned int index0 = (ir1 * (ir1+1))/2;
            for (unsigned int ls2=0;ls2<nls2;++ls2)
            {
                unsigned int ir2 = off2 + ls2;
                if ( ir2 > ir1 ) break;
                norms[len] = aux.normalization_factor(lv1,ls1) * aux.normalization_factor(lv2,ls2);
                lstates[len] = (ls1 << 4U) + ls2;
                orbs[len] = index0 + ir2;
                ++len;
            }
        }
        return len;    
    }
    unsigned int precalc(const AuxFunctions& aux,
        unsigned int off1, unsigned int off2, double *factors)
    {
        int nls1 = aux.number_of_lstates(lv1);
        int nls2 = aux.number_of_lstates(lv2);
        len = 0;
        for (unsigned int ls1=0;ls1<nls1;++ls1)
        {
            unsigned int ir1 = off1 + ls1;
            unsigned int index0 = (ir1 * (ir1+1))/2;
            for (unsigned int ls2=0;ls2<nls2;++ls2)
            {
                unsigned int ir2 = off2 + ls2;
                if ( ir2 > ir1 ) break;
                norms[len] = aux.normalization_factor(lv1,ls1) * aux.normalization_factor(lv2,ls2);
                lstates[len] = (ls1 << 4U) + ls2;
                orbs[len] = index0 + ir2;
                factors[len] = (ir1!=ir2) ? 4.:2.;
                ++len;
                
            }
        }
        return len;    
    }
};

/**
struct ShellQuartet
{
    double * norms;
    unsigned int * lstates;
    unsigned int ** orbs;
    //
    const double * a;
    const double * b;
    const double * c;
    const double * d;
    const double * al1;
    const double * co1;
    const double * al2;
    const double * co2;
    const double * al3;
    const double * co3;
    const double * al4;
    const double * co4;
    double ab2;
    double cd2;
    int npr1,lv1;
    int npr2,lv2;
    int npr3,lv3;
    int npr4,lv4;
    int switch12,switch34;
    unsigned int len;    
    unsigned int maxints;

    ShellQuartet() = delete;
    
    ShellQuartet(int maxl)
    {
        int maxlst = utils::number_of_lstates(maxl);
        maxints = maxlst * maxlst;
        maxints = maxints * maxints;
        orbs = new_matrix<unsigned int>(maxints,4);
        norms = new double[maxints];
        lstates = new unsigned int[maxints];
        switch12 = 0;
        switch34 = 0;
        len = 0;
    }
    
    ~ShellQuartet()
    {
        delete [] lstates;
        delete [] norms;
        delete_matrix(orbs,maxints);
    }
    
    void assign1(const Shell& sh, const Center * center)
    {
        npr1 = sh.number_of_prims();
        lv1 = sh.Lvalue();
        int cn1 = sh.center();
        al1 = sh.alf_ptr();
        co1 = sh.cof_ptr();
        a = (center+cn1)->r_vec();
    }
    void assign3(const Shell& sh, const Center * center)
    {
        npr3 = sh.number_of_prims();
        lv3 = sh.Lvalue();
        int cn3 = sh.center();
        al3 = sh.alf_ptr();
        co3 = sh.cof_ptr();
        c = (center+cn3)->r_vec();
    }
    void assign2(const Shell& sh, const Center * center)
    {
        lv2 = sh.Lvalue();
        if (lv2 > lv1) {
            switch12 = 1;
            int t = lv2;
            lv2 = lv1;
            lv1 = t;
            npr2 = npr1;
            al2 = al1;
            co2 = co1;
            b = a;
            npr1 = sh.number_of_prims();
            int cn1 = sh.center();
            al1 = sh.alf_ptr();
            co1 = sh.cof_ptr();
            a = (center+cn1)->r_vec();
        }else{
            npr2 = sh.number_of_prims();
            lv2 = sh.Lvalue();
            int cn2 = sh.center();
            al2 = sh.alf_ptr();
            co2 = sh.cof_ptr();
            b = (center+cn2)->r_vec();
        }
        ab2 = dist_sqr(a,b);
    }
    void assign4(const Shell& sh, const Center * center)
    {
        lv4 = sh.Lvalue();
        if (lv4 > lv3) {
            switch34 = 1;
            int t = lv4;
            lv4 = lv3;
            lv3 = t;
            npr4 = npr3;
            al4 = al3;
            co4 = co3;
            d = c;
            npr3 = sh.number_of_prims();
            int cn3 = sh.center();
            al3 = sh.alf_ptr();
            co3 = sh.cof_ptr();
            c = (center+cn3)->r_vec();
        }else{
            npr4 = sh.number_of_prims();
            int cn4 = sh.center();
            lv4 = sh.Lvalue();
            al4 = sh.alf_ptr();
            co4 = sh.cof_ptr();
            d = (center+cn4)->r_vec();
        }
        cd2 = dist_sqr(c,d);
    }
    void unswitch34()
    {
        if ( switch34) {
            c = d;
            npr3 = npr4;
            lv3 = lv4;
            al3 = al4;
            co3 = co4;
        }
        switch34 = 0;
    }
    void unswitch12()
    {
        if ( switch12) {
            a = b;
            npr1 = npr2;
            lv1 = lv2;
            al1 = al2;
            co1 = co2;
        }
        switch12 = 0;
    }
    unsigned int precalc(const AuxFunctions& aux,
        unsigned int off1, unsigned int off2, unsigned int off3, unsigned int off4)
    {
        int nls1 = aux.number_of_lstates(lv1);
        int nls2 = aux.number_of_lstates(lv2);
        int nls3 = aux.number_of_lstates(lv3);
        int nls4 = aux.number_of_lstates(lv4);
        len = 0;
        for (unsigned int ls1=0;ls1<nls1;++ls1)
        {
            unsigned int ir1 = off1 + ls1;
            for (unsigned int ls2=0;ls2<nls2;++ls2)
            {
                unsigned int ir2 = off2 + ls2;
                if ( ir2 > ir1 ) break;
                double n12 = aux.normalization_factor(lv1,ls1) * aux.normalization_factor(lv2,ls2);
                unsigned int ls12 = (ls1 << 4U) + ls2;
                for (unsigned int ls3=0;ls3<nls3;++ls3)
                {
                    unsigned int ir3 = off3 + ls3;
                    if (ir3 > ir1) break;
                    for (unsigned int ls4=0;ls4<nls4;++ls4)
                    {
                        unsigned int ir4 = off4 + ls4;
                        if ( ir4 > ir3 ) break;
                        if ( ir1 = ir3 && ir4 > ir3) break;
                        norms[len] = n12 * aux.normalization_factor(lv1,ls1) * aux.normalization_factor(lv2,ls2);
                        unsigned int ls34 = (ls3 << 4U) + ls4;
                        lstates[len] = (ls12 << 8U) + ls34;
                        orbs[len][0] = ir1;
                        orbs[len][1] = ir2;
                        orbs[len][2] = ir3;
                        orbs[len][3] = ir4;
                        ++len;
                    }
                }
            }          
        }	
        return len;    
    }
};
*/

}
#endif

