#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <string>
#include <cmath>
using namespace std;

int main()
{
int i;
double fnrg,snrg;
char *file;
string filename;
string molcase;
string first="short.dat.";
string dot=".";
string sec[3]={"3g","431","631"};
string mol[7]={"h2","n2","co","hf","h2o","nh3","ch4"};
string spstr[7]={" "," "," "," ","","",""};
double nrgs[7][3]={ {-1.1167,-1.127,-1.131},
                    {-107.496,-108.754,-108.943},
                    {-111.225,-112.552,-112.737},
                    {-98.571,-99.887,-100.011},
                    {-74.963,-75.907,-76.023},
                    {-55.454,-56.102,-56.195},
                    {-39.727,-40.140,-40.202}};

    ofstream out("nrgs.out");
    ifstream in;
    out<<"MOL"<<" Basis "<<"           energy  "<<"    Szabo&Ostlund"<<endl;
    for (i=0;i<7;i++)
    {
     molcase=mol[i];
     filename=((first+sec[0])+dot)+molcase;
     file=(char*)filename.c_str();
     in.open(file);
     in>>fnrg;
     in>>snrg;
     in.close();
     out.setf(ios::showpoint);
     out.setf(ios::uppercase);
     out<<molcase<<spstr[i]<<"  ";
     out<<sec[0]<<" "<<setw(20)<<setprecision(10)<<snrg;
     out<<setw(15)<<setprecision(7)<<nrgs[i][0]<< " ";
     out<<setw(20)<<setprecision(7)<<std::scientific << fabs((snrg - nrgs[i][0])) << "\n";
     filename=((first+sec[1])+dot)+molcase;
     file=(char*)filename.c_str();
     in.open(file);
     in>>fnrg;
     in>>snrg;
     in.close();
     out.setf(ios::showpoint);
     out.setf(ios::uppercase);
     out<<molcase<<spstr[i]<<"  ";
     out<<sec[1]<<setw(20)<<setprecision(10)<<snrg;
     out<<setw(15)<<setprecision(7)<<nrgs[i][1]<< " ";
     out<<setw(20)<<setprecision(7)<<std::scientific << fabs((snrg - nrgs[i][1])) << "\n";
     filename=((first+sec[2])+dot)+molcase;
     file=(char*)filename.c_str();
     in.open(file);
     in>>fnrg;
     in>>snrg;
     out.setf(ios::showpoint);
     out.setf(ios::uppercase);
     out<<molcase<<spstr[i]<<"  ";
     out<<sec[2]<<setw(20)<<setprecision(10)<<snrg;
     out<<setw(15)<<setprecision(7)<<nrgs[i][2]<< " ";
     out<<setw(20)<<setprecision(7)<<std::scientific << fabs((snrg - nrgs[i][2])) << "\n";
     in.close();
    }   
    out.close();
    return 0;
} 


