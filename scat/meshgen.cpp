#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "radial_mesh.hpp"

int main()
{
    std::ifstream in("mesh.in");
    unomol::radial_mesh mesh;
    in >> mesh;
    in.close();
    std::cout << "# of points = " << mesh.num_points() << "\n";
    std::cout << "# of regions = " << mesh.num_regions() << "\n";
    std::ofstream out("pos.grid.in");
    std::vector<double> rpts;
    mesh.fill_points(rpts);
    out << std::setw(12) << mesh.num_points() << "\n";
    for (double r : rpts) {
        out << std::setw(15) << std::setprecision(8) << std::fixed << r << " ";
        out << std::setw(15) << std::setprecision(8) << std::fixed << 0 << " ";
        out << std::setw(15) << std::setprecision(8) << std::fixed << 0 << "\n"; 
    }
    out.close();
}




