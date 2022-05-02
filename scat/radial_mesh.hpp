#ifndef UNOMOL_RADIAL_MESH_HPP
#define UNOMOL_RADIAL_MESH_HPP
#include <cstdlib>
#include <vector>
#include <iostream>
#include <cmath>

namespace unomol {

class radial_mesh {
    std::vector<double> rst;
    std::vector<double> delr;
    std::vector<std::size_t> nst;
    std::size_t npts;
    std::size_t nreg;
public:

    constexpr double rstart(std::size_t i) const noexcept { return rst[i];}
    constexpr double deltar(std::size_t i) const noexcept { return delr[i];}
    constexpr double begin_rad() const noexcept { return rst[0];}
    constexpr double end_rad() const noexcept { return rst[nreg];}
    constexpr std::size_t nsteps(std::size_t i) const noexcept { return nst[i];}
    constexpr std::size_t num_points() const noexcept { return npts;}
    constexpr std::size_t num_regions() const noexcept { return nreg;}    
    
    std::size_t fill_points(std::vector<double>& rpts) const noexcept
    {
        rpts.clear();
        for (std::size_t k=0;k<nreg;++k) {
            double rs = rst[k];
            double dr = delr[k];
            std::size_t ns = nst[k];
            for (std::size_t j=0;j<ns;++j) {
                double r = rs + j * dr;
                rpts.push_back(r);
            }
        }	
        rpts.push_back(rst[nreg]);
        return rpts.size();
    }
    
    friend std::istream& operator >> ( std::istream& is, radial_mesh& mesh)
    {
        mesh.rst.clear();
        mesh.delr.clear();
        mesh.nst.clear();
        mesh.npts = 0;
        is >> mesh.nreg;
        double rs,re,dr;
        std::size_t ns;
        is >> rs;
        for (std::size_t k=0;k<mesh.nreg;++k) {
            is >> re;
            is >> dr;
            if ( re <= rs) {
                std::cerr << "your mesh is out of order in mesh.in!\n";
                std::cerr << " end = " << re << "\n";
                std::cerr << " begin = " << rs << "\n";
                exit(-1);
            }
            ns = static_cast<std::size_t>(rint((re-rs)/dr));
            mesh.rst.push_back(rs);
            mesh.delr.push_back(dr);
            mesh.nst.push_back(ns);
            mesh.npts+=ns;
            rs = re;
        }
        is >> rs;   
        mesh.rst.push_back(rs);
        mesh.npts+=1;
        return is; 
    }
    
    friend std::ostream& operator << ( std::ostream& os, const radial_mesh& mesh)
    {
        os << mesh.nreg;
        os << " " << mesh.rst[0] << "\n";
        for (std::size_t k=0;k<mesh.nreg;++k) {
            os << " " << mesh.rst[k+1] << " " << mesh.delr[k] << "\n";
        }    
        return os;
    }
};

}
#endif
