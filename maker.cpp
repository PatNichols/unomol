#include <cstdlib>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

struct rterm
{
    long v;
    int m;
    int px;
    int py;
    int pz;
    int lx;
    int ly;
    int lz;
};

typedef struct rterm rterm;

std::ostream& operator << ( std::ostream& os, const rterm& r)
{
    if ( r.v > 1L) os << r.v << "*";
    if ( r.px ) os << " x" << r.px << "*";
    if ( r.py ) os << " y" << r.py << "*";
    if ( r.pz ) os << " z" << r.pz << "*";
    os << " r["<<r.lx<<"]["<<r.ly<<"]["<<r.lz<<"]["<<r.m<<"]";
    return os;
}

bool recur_term(const rterm& rc,std::vector<rterm>& ro)
{
    rterm r1 = rc;
    rterm r2 = rc;
    r1.m +=1;
    r2.m +=1;
    std::cout << "recur :"<<rc<< " => ";
    if ( rc.lz ) {
        r1.lz -= 1;
        r1.pz += 1;
        ro.push_back(r1);
        std::cout << r1;
        if ( rc.lz > 1) {
            r2.lz -= 2;
            r2.v *= ( rc.lz-1);
            ro.push_back(r2);
            std::cout << " + " << r2;
        }
        std::cout << "\n";
        return true;
    }
    if ( rc.ly ) {
        r1.ly -= 1;
        r1.py += 1;
        ro.push_back(r1);
        std::cout << r1;
        if ( rc.ly > 1) {
            r2.ly -= 2;
            r2.v *= ( rc.ly-1);
            ro.push_back(r2);
            std::cout << " + " << r2;
        }
        std::cout << "\n";
        return true;
    }
    if ( rc.lx ) {
        r1.lx -= 1;
        r1.px += 1;
        ro.push_back(r1);
        std::cout << r1;
        if ( rc.lx > 1) {
            r2.lx -= 2;
            r2.v *= ( rc.lx-1);
            ro.push_back(r2);
            std::cout << " + " << r2;
        }
        std::cout << "\n";
        return true;
    }
    ro.push_back(rc);
    std::cout << " prim \n";
    return false;
}

void print_rterms(int lx,int ly,int lz,const std::vector<rterm>& rvecs)
{
    std::cout <<"------- FINAL -----------\n";
    std::cout <<"    " << lx << " " << ly << " " << lz << "\n";
    for ( const rterm& r : rvecs)
    {
        std::cout << r << "\n";
    }
}

std::vector<rterm> form_rterms(int lx,int ly,int lz)
{
    rterm rc;
    std::vector<rterm> ro;
    std::vector<rterm> rn;

    rc.v = 1L;
    rc.lx = lx;
    rc.ly = ly;
    rc.lz = lz;
    rc.m = 0;
    rc.px = 0;
    rc.py = 0;
    rc.pz = 0;
    ro.push_back(rc);
    if (!(lx+ly+lz)) return ro;
    std::cout <<"-----------------------\n";
    std::cout <<"    " << lx << " " << ly << " " << lz << "\n";
    while (1)
    {
        std::cout <<" __________________ \n";
        rn.clear();
        bool done = true;
        for (const rterm& r : ro)
        {
            std::cout << "oooooooooooooooooooooooooooooooooo  \n";
            std::cout << "    new " << r << "\n";
            bool not_done = recur_term(r,rn);
            if (done && not_done) done = false;
        }
        std::cout << ro.size() << " " << rn.size() << "\n";
        if (done) {
            print_rterms(lx,ly,lz,ro);
            return ro;
        }
        ro = rn;
    }
    return ro;
}

void consolidate_rterms(std::vector<rterm>& rvecs)
{
    std::vector<rterm> rnew;
    std::vector<int> taken(rvecs.size(),0);
    for (size_t k=0; k<rvecs.size(); ++k)
    {
        if ( taken[k] ) continue;
        rterm rn = rvecs[k];
        for (size_t j=k+1; j<rvecs.size(); ++j)
        {
            const rterm& rj = rvecs[j];
            if ( rn.m == rj.m) {
                if ( rn.px == rj.px && rn.py == rj.py && rn.pz == rj.pz)
                {
                    rn.v += rj.v;
                    taken[j] = 1;
                }
            }
        }
        rnew.push_back(rn);
    }
    rvecs.swap(rnew);
}

void print_expr0(const std::vector<rterm>& rvecs,std::ostream& out)
{
    for (size_t i=0; i<rvecs.size(); ++i)
    {
        rterm r = rvecs[i];
        int mi = r.m;
        if ( i!=0) out << " + ";
        if ( r.v != 1L) out << r.v << "*";
        if ( r.px ) out << "x" << r.px << "*";
        if ( r.py ) out << "y" << r.py << "*";
        if ( r.pz ) out << "z" << r.pz << "*";
        out << "rm[" << mi << "]";
    }
    out << ";\n";
}

void print_expr(const std::vector<rterm>& rvecs,std::ostream& out)
{
    std::vector<int> taken(rvecs.size(),0);
    std::vector<std::string> sterms;
    std::ostringstream sout;
    for (size_t i=0; i<rvecs.size(); ++i)
    {
        if ( taken[i] ) continue;
        sterms.clear();
        bool need_paran = false;
        int cnt = 0;
        rterm r = rvecs[i];
        int mi = r.m;
        if ( i!=0) out << " + ";
        if ( r.v!= 1L) sterms.push_back(std::to_string(r.v));
        if ( r.px!= 0) {
            if ( r.v!=1L) sterms.push_back(std::string("*"));
            sterms.push_back(std::string("x"));
            sterms.push_back(std::to_string(r.px));
        }
        if ( r.py!= 0) {
            if ( r.v!=1L || r.px ) sterms.push_back(std::string("*"));
            sterms.push_back(std::string("y"));
            sterms.push_back(std::to_string(r.py));
        }
        if ( r.pz!= 0) {
            if ( r.v!=1L || r.px || r.py ) sterms.push_back(std::string("*"));
            sterms.push_back(std::string("z"));
            sterms.push_back(std::to_string(r.pz));
        }
//        if ( r.v == 1L && r.px==0 && r.py==0 && r.pz==0)
//        {
//            sterms.push_back("1");
//        }  
        for (size_t j=i+1; j<rvecs.size(); ++j)
        {
            if ( rvecs[j].m == mi )
            {
                const rterm& rj = rvecs[j];
                taken[j] = 1;
                sterms.push_back("+");
                if ( rj.v!= 1L) sterms.push_back(std::to_string(rj.v));
                if ( rj.px!= 0) {
                    if ( rj.v!=1L) sterms.push_back(std::string("*"));
                    sterms.push_back(std::string("x"));
                    sterms.push_back(std::to_string(rj.px));
                }
                if ( rj.py!= 0) {
                    if ( rj.v!=1L || rj.px ) sterms.push_back(std::string("*"));
                    sterms.push_back(std::string("y"));
                    sterms.push_back(std::to_string(rj.py));
                }
                if ( rj.pz!= 0) {
                    if ( rj.v!=1L || rj.px || rj.py ) sterms.push_back(std::string("*"));
                    sterms.push_back(std::string("z"));
                    sterms.push_back(std::to_string(rj.pz));
                }
                if ( rj.v == 1L && rj.px==0 && rj.py==0 && rj.pz==0)
                {
                    sterms.push_back("1");
                } 
                need_paran = true;
            }
        }
        if (need_paran)
        {
            sterms.push_back(std::string(")"));
            out << "(";
        }
        for (const std::string& s  : sterms)
        {
            out << s;
        }       
        if (r.v != 1L || (r.px + r.py + r.pz)!=0) out << "*";
        out << "rm[" << mi << "]";
    }
    out << ";\n";
}


int main(int argc, char **argv)
{
    int lmax = 4;
    std::ofstream out("rfunx.cpp");

    if ( argc != 1) {
        lmax = std::stoi(argv[1]);
    }
    out << "void rfunction::unroll_eval(double ***r, const double * rm, int ltot)\n{\n";
    out << "   r[0][0][0] = rm[0];\n";
    out << "   if ( ltot == " << 0 << ") return;\n";
    for (int ltot = 1; ltot <= lmax; ++ltot) {
        out <<"\n";
        if ( ltot!=1) {
            out << "   double x" <<(ltot) << "= x" << (ltot-1) << "*x1;\n";
            out << "   double y" <<(ltot) << "= y" << (ltot-1) << "*y1;\n";
            out << "   double z" <<(ltot) << "= z" << (ltot-1) << "*z1;\n";
        }else{
            out << "   double x1 = pq[0];\n";
            out << "   double y1 = pq[1];\n";
            out << "   double z1 = pq[2];\n";
        }
        for (int lx=0; lx<=lmax; ++lx)
        {
            int lyend = lmax - lx;
            for (int ly=0; ly<=lyend; ++ly)
            {
                int lzend = lyend - ly;
                for (int lz=0; lz<=lzend; ++lz)
                {
                    int lsum = lx + ly + lz;
                    if ( lsum!=ltot) continue;
                    std::vector<rterm> rvecs = std::move(form_rterms(lx,ly,lz));
                    std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n";
                    std::cout << "rvecs size = " << rvecs.size() << "\n";
                    consolidate_rterms(rvecs);
                    out << "   r[";
                    out << lx << "][" << ly << "][" << lz << "] = ";
                    print_expr(rvecs,out);
                }
            }
        }
        out <<"   if ( ltot == " << ltot << ") return;\n";
    } 
    out << "}\n";
    out.close();
    return EXIT_SUCCESS;
}