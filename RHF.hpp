#ifndef UNOMOL_RHF_HPP
#define UNOMOL_RHF_HPP
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include "Util.hpp"
#include "Basis.hpp"
#include "TwoElectronInts.hpp"
#include "AuxFunctions.hpp"
#include "GDPMInts.hpp"
#include "Moments.hpp"
#include "OneElectronInts.hpp"
#include "SymmPack.hpp"
#include "FField.hpp"
#include "Stopwatch.hpp"
using namespace std;

namespace unomol {

class RestrictedHartreeFock {
  public:
    RestrictedHartreeFock() = delete;
    RestrictedHartreeFock(Basis* b,
                          TwoElectronInts* t):basis(*b),tints(*t) {
        nshell=basis.number_of_shells();
        no=basis.number_of_orbitals();
        no2=no*(no+1)/2;
        ncen=basis.number_of_centers();
        tnshell=basis.total_number_of_shells();
        tno=basis.total_number_of_orbitals();
        tno2=tno*(tno+1)/2;
        tncen=basis.total_number_of_centers();
        nocc=basis.number_of_electrons();
        if (nocc%2) {
            fatal_error("Odd number of electrons in RHF");
        }
        nocc/=2;
        maxits=basis.maximum_iterations();
        eps=basis.scf_eps();
        scf_accel=basis.scf_flags(1);
        cflag=basis.scf_flags(0);
        ediff=10.0;
        pdiff=10.0;
        eold=0.0;
        nucrep=0.0;
        energy=0.0;
        energyGs=0.0;
        int tnsqr=tno*tno;
        PmatGs=new double[tno2];
        Pold2=new double[tno2];
        Pold=new double[tno2];
        Pmat=new double[tno2];
        Gmat=new double[tno2];
        Hmat=new double[tno2];
        Fock=new double[tno2];
        Xmat=new double[tnsqr];
        Wmat=new double[tnsqr];
        Cmat=new double[tnsqr];
        Tmat=new double[tno2];
        Smat=new double[tno2];
        Evals=new double[tno];
        sEvals=new double[tno];
        Wrka=new double[tnsqr];
    }

    ~RestrictedHartreeFock() {
        delete [] Wrka;
        delete [] sEvals;
        delete [] Evals;
        delete [] Smat;
        delete [] Tmat;
        delete [] Cmat;
        delete [] Wmat;
        delete [] Xmat;
        delete [] Fock;
        delete [] Hmat;
        delete [] Gmat;
        delete [] Pmat;
        delete [] Pold;
        delete [] Pold2;
        delete [] PmatGs;
    }

    void update() noexcept {
        int i;
        for (i=0; i<no2; ++i) Gmat[i]=0.0;
        tints.formGmatrix(Pmat,Gmat);
        energy=SymmPack::TraceSymmPackProduct(Pmat,Hmat,no)*2.0+
               SymmPack::TraceSymmPackProduct(Pmat,Gmat,no);
        ediff=energy-eold;
        eold=energy;
        for (i=0; i<no2; ++i) Fock[i]=Hmat[i]+Gmat[i];
        SymmPack::sp_trans(no,Fock,Xmat,Wrka);
        SymmPack::rsp(no,Fock,Wmat,Evals,Wrka);
        formCmatrix(Cmat);
        if (scf_accel==1) memcpy(Pold2,Pold,sizeof(double)*no2);
        memcpy(Pold,Pmat,sizeof(double)*no2);
        formPmatrix(Pmat,Cmat,nocc);
        pdiff=SymmPack::SymmPackDiffNorm(Pmat,Pold,no);
        ++iteration;
    }

    void dpm_update(TwoElectronInts& xints) noexcept {
        int i;
        for (i=0; i<tno2; ++i) Gmat[i]=0.0;
        xints.formGmatrix(Pmat,Gmat);
        tints.formGmatrix(Pmat,Gmat);
        energy=SymmPack::TraceSymmPackProduct(Pmat,Hmat,no)*2.0+
               SymmPack::TraceSymmPackProduct(Pmat,Gmat,no);
        ediff=energy-eold;
        eold=energy;
        for (i=0; i<no2; ++i) Fock[i]=Hmat[i]+Gmat[i];
        SymmPack::sp_trans(no,Fock,Xmat,Wrka);
        SymmPack::rsp(no,Fock,Wmat,Evals,Wrka);
        formCmatrix(Cmat);
        if (scf_accel==1) memcpy(Pold2,Pold,sizeof(double)*no2);
        memcpy(Pold,Pmat,sizeof(double)*no2);
        formPmatrix(Pmat,Cmat,nocc);
        pdiff=SymmPack::SymmPackDiffNorm(Pmat,Pold,no);
        ++iteration;
    }

    void findEnergy() noexcept {
        nucrep=nuclear_repulsion_energy(ncen,basis.center_ptr());
        OneElectronInts(basis,Smat,Tmat,Hmat);
        MomentInts(basis);
        OintsOutput();
        formXmatrix();
        if (basis.scf_flags(2)) {
            FILE* in=open_file("PMATRIX.DAT");
            Fread(Pmat,sizeof(double),no2,in);
            fclose(in);
        } else {
            PmatrixGuess();
        }

        fprintf(stderr,"convergence based upon ");
        switch (cflag) {
        case 0:
            fprintf(stderr,"energy\n");
            break;
        case 1:
            fprintf(stderr,"density\n");
            break;
        default:
            fprintf(stderr,"both energy and density\n");
            break;
        }
        putils::Stopwatch timer;
        timer.start();
        iteration=0;
        extrap=0;
        eold=0.0;
        update();
        fprintf(stderr,"Iteration     =  %6d\n",iteration);
        fprintf(stderr,"Energy        =  %25.16le\n",(eold+nucrep));
        fprintf(stderr,"Delta Energy  =  %25.15le\n",ediff);
        fprintf(stderr,"Delta Density =  %25.15le\n\n",pdiff);
        double init_energy=energy+nucrep;
        iteration=1;
        update();
        fprintf(stderr,"Iteration     =  %6d\n",iteration);
        fprintf(stderr,"Energy        =  %25.16le\n",(eold+nucrep));
        fprintf(stderr,"Delta Energy  =  %25.15le\n",ediff);
        fprintf(stderr,"Delta Density =  %25.15le\n\n",pdiff);
        iteration=2;
        while (iteration<maxits) {
            scf_converger();
            update();
            if (is_converged()) break;
            fprintf(stderr,"Iteration     =  %6d\n",iteration);
            fprintf(stderr,"Energy        =  %25.16le\n",(eold+nucrep));
            fprintf(stderr,"Delta Energy  =  %25.15le\n",ediff);
            fprintf(stderr,"Delta Density =  %25.15le\n\n",pdiff);
        }
        timer.stop();
        std::cerr << "SCF time = " << timer.elapsed_time() << " s\n";
        memcpy(PmatGs,Pmat,sizeof(double)*no2);
        for (int j=no2; j<tno2; ++j) PmatGs[j]=0.0;
        energyGs=energy+nucrep;
        FILE *fp=create_file("PMATRIX.DAT");
        Fwrite(Pmat,sizeof(double),no2,fp);
        fclose(fp);
        final_output(init_energy);
    }

    void formCmatrix(double* c) {
        copy_trans(no,Wmat);
        for (int i=0; i<no; ++i) {
            const double *Xi=Xmat+i*no;
            for (int j=0; j<no; ++j) {
                const double *xp=Xi;
                const double *wp=Wmat+j*no;
                double sum=0.0;
                for (int k=0; k<no; ++k) sum+=xp[k]*wp[k];
                (*(c+i*no+j))=sum;
            }
        }
    }

    void formPmatrix(double* p,double *c,int noc) {
        int ij=0;
        for (int i=0; i<no; ++i) {
            const double *ci=c+i*no;
            for (int j=0; j<=i; ++j,++ij) {
                const double *cj=c+j*no;
                double sum=0.0;
                for (int k=0; k<noc; ++k) sum+=ci[k]*cj[k];
                p[ij]=sum;
            }
        }
    }

    void PmatrixGuess() {
        for (int i=0; i<no2; ++i) Fock[i]=Hmat[i];
        SymmPack::sp_trans(no,Fock,Xmat,Wrka);
        SymmPack::rsp(no,Fock,Wmat,Evals,Wrka);
        formCmatrix(Cmat);
        formPmatrix(Pmat,Cmat,nocc);
    }


    void formXmatrix() {
        memcpy(Fock,Smat,sizeof(double)*no2);
        SymmPack::rsp(no,Fock,Wmat,Evals,Wrka);
        for (int i=0; i<no; ++i) {
            double f=Evals[i];
            if (fabs(f)<1.e-7) {
                if (f<0.0) f=-f;
            }
            if (f==0.0) fatal_error("Zero Eigenvalue in Overlap Matrix");
            if (f<0.0) fatal_error("Negative Eigenvalue in Overlap Matrix");
            sEvals[i]=1.0/sqrt(f);
        }
        for (int i=0; i<no; ++i) {
            double *Xi=Xmat+i*no;
            const double *Wi=Wmat+i*no;
            for (int j=0; j<no; ++j) {
                Xi[j]=Wi[j]*sEvals[j];
            }
        }
    }

    void FiniteFieldAnalysis() {
        double polnrg[3],alfpol[3];
        double Emag=basis.FiniteFieldValue();
        double Efield[3];
        for (int ix=0; ix<3; ++ix) {
            Efield[2]=Efield[1]=Efield[0]=0.0;
            Efield[ix]=Emag;
            OneElectronInts(basis,Smat,Tmat,Hmat);
            FiniteFieldMatrix(basis,Hmat,Efield);
            formXmatrix();
            memcpy(Pmat,PmatGs,sizeof(double)*no2);
            iteration=0;
            extrap=0;
            eold=0.0;
            update();
            double init_energy=energy+nucrep;
            while (iteration<maxits) {
                scf_converger();
                update();
                if (is_converged()) break;
            }
            polnrg[ix]=energy+nucrep-energyGs;
            alfpol[ix]= -polnrg[ix]*2.0/Emag/Emag;
        }
        FILE *out=create_file("finitefield.out");
        fprintf(out,"        Unomol Finite Field Analysis\n");
        fprintf(out,"        Electric field magnitude = %20.12le\n",Emag);
        fprintf(out,"     alfa                pol energy\n");
        fprintf(out," x   %20.12le %20.12le\n",alfpol[0],polnrg[0]);
        fprintf(out," y   %20.12le %20.12le\n",alfpol[1],polnrg[1]);
        fprintf(out," z   %20.12le %20.12le\n",alfpol[2],polnrg[2]);
        double alfiso= (alfpol[0]+alfpol[1]+alfpol[2])/3.0;
        double alfan=  (2.*alfpol[2]-alfpol[1]-alfpol[0])/3.0;
        fprintf(out," alf isotropic = %20.12le\n",alfiso);
        fprintf(out," alf aniso     = %20.12le\n",alfan);
        fclose(out);
    }

    double nuclear_repulsion_energy(int ncen,
                                    const Center* center) {
        double sum=0.0;
        double q1,q2,r12;
        const double *r1,*r2;

        for (int i=0; i<ncen; ++i) {
            q1=(center+i)->charge();
            r1=(center+i)->r_vec();
            for (int j=i+1; j<ncen; ++j) {
                q2=(center+j)->charge();
                r2=(center+j)->r_vec();
                r12=sqrt(dist_sqr(r1,r2));
                sum+=q1*q2/r12;
            }
        }
        return sum;
    }

    void findPolarizationPotential() {
        double px,py,pz;
        const string xstr("XINTS.DAT");
        int sshell=nshell;
        nshell=tnshell;
        no=tno;
        no2=tno2;
        ncen=tncen;
        int pcen=basis.skip_center();
        int npts;
        basis.dpm_augment();
        eps=1.e-12;
        ////////////////////////////////////////////////////////
        // start calculation
        std::ifstream in("pos.grid.dat");
        in >> npts;
        in >> px;
        in >> py;
        in >> pz;
        basis.SetCenterPosition(px,py,pz,pcen);
        nucrep=nuclear_repulsion_energy(ncen,
                                        basis.center_ptr());
        TwoElectronInts xints(basis,sshell,xstr);
        OneElectronInts(basis,Smat,Tmat,Hmat);
        GDPMInts(basis,Hmat);
        OintsOutput();
        formXmatrix();
        memcpy(Pmat,PmatGs,sizeof(double)*no2);
        iteration=0;
        extrap=0;
        eold=0.0;
        dpm_update(xints);
        double init_energy=energy+nucrep;
        while (iteration<maxits) {
            scf_converger();
            dpm_update(xints);
            if (is_converged()) break;
        }

        double final_energy=energy+nucrep;
        double vpol=final_energy-init_energy;
        double vstat=init_energy-energyGs+nucrep;
        double r2=px*px+py*py+pz*pz;
        double alfa= -2.0*vpol*r2*r2;
        FILE *vout = create_file("vpol.out");
        FILE *sout = create_file("spol.out");
        fprintf(vout,"%15.10lf %15.10lf %15.10lf %25.15le %25.15le\n",px,py,pz,
                vpol,alfa);
        fflush(vout);
        fprintf(sout,"%3d %20.10le %20.10le %20.10le %20.10le %25.15le\n",
                0,energyGs,init_energy,final_energy,vstat,ediff);
        fflush(sout);
        for (int i=1; i<npts; ++i) {
            in >> px;
            in >> py;
            in >> pz;
            basis.SetCenterPosition(px,py,pz,pcen);
            nucrep=nuclear_repulsion_energy(basis.number_of_centers(),
                                            basis.center_ptr());
            xints.recalculate(basis);
            OneElectronInts(basis,Smat,Tmat,Hmat);
            GDPMInts(basis,Hmat);
            formXmatrix();
            memcpy(Pmat,PmatGs,sizeof(double)*no2);
            iteration=0;
            extrap=0;
            eold=0.0;
            dpm_update(xints);
            double init_energy=energy+nucrep;
            while (iteration<maxits) {
                scf_converger();
                dpm_update(xints);
                if (is_converged()) break;
            }
            final_energy=energy+nucrep;
            vpol=final_energy-init_energy;
            vstat=init_energy-energyGs;
            r2=px*px+py*py+pz*pz;
            alfa= -2.0*vpol*(r2*r2);
            fprintf(vout,"%15.10lf %15.10lf %15.10lf %25.15le %25.15le\n",px,py,pz,
                    vpol,alfa);
            fflush(vout);
            fprintf(sout,"%3d %20.10le %20.10le %20.10le %20.10le %25.15le\n",
                    i,energyGs,init_energy,final_energy,vstat,ediff);
            fflush(sout);
        }
        fclose(vout);
        fclose(sout);
        in.close();
    }

    constexpr bool is_converged() noexcept {
        switch (cflag) {
        case 0:
            if (ediff<eps) return true;
            return false;
        case 1:
            if (pdiff<eps) return true;
            return false;
        default:
            if (pdiff<eps && ediff<eps) return true;
            return false;
        }
    }

    void final_output(double init_energy) {
        FILE *out=create_file("short.gs.out");
        fprintf(out,"%25.16le\n",init_energy);
        fprintf(out,"%25.16le\n",energy+nucrep);
        fprintf(out,"%25.16le\n",ediff);
        fclose(out);
        out=create_file("scfout.gs.out");
        double trace_t=2.0*SymmPack::TraceSymmPackProduct(Pmat,Tmat,no);
        double virial=fabs((energy+nucrep-trace_t)/(trace_t)/2.0);
        if (maxits<=iteration) fprintf(out,"WARNING CONVERGENCE _NOT_ REACHED!\n");
        fprintf(out,"Final Iteration              = %12u\n",iteration);
        fprintf(out,"Hartree Fock Energy          = %25.15le Hartree\n",energy+nucrep);
        fprintf(out,"Electronic Energy            = %25.15le Hartree\n",energy);
        fprintf(out,"Nuclear Rep. Energy          = %25.15le Hartree\n",nucrep);
        fprintf(out,"Kinetic Energy               = %25.15le Hartree\n",trace_t);
        fprintf(out,"Difference in Final Energies = %25.16le Hartree\n",ediff);
        fprintf(out,"Pmatrix diff norm            = %25.16le\n",pdiff);
        fprintf(out,"virial                       = %25.16le\n",virial);
        fprintf(out,"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
        fprintf(out,"               Orbital Energies\n");
        fprintf(out,"Orbital Energy          Occupancy\n");
        for (int i=0; i<nocc; i++) {
            fprintf(out,"%7u %25.16le %12u\n",i+1,Evals[i],2);
        }
        for (int i=nocc; i<no; i++) {
            fprintf(out,"%7u %25.16le %12u\n",i+1,Evals[i],0);
        }
        fprintf(out,"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
        PopulationAnalysis(out);
        if (basis.prt_flags(0)) {
            fprintf(out,"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
            fprintf(out,"              Fock Eigenvectors\n");
            fprintf(out,"\n");
            for (int i=0; i<no; ++i) {
                int occup=(i<nocc)?2:0;
                fprintf(out,"Column %5d Occupancy= %2d EigenEnergy= %25.15le\n",i,
                        occup,Evals[i]);
                int nex=no-4*no/4;
                int end=4*(no/4);
                const double *Cm=Cmat+i;
                for (int j=0; j<no; j+=4) {
                    fprintf(out,"%18.10le %18.10le %18.10le %18.10le\n",
                            (*Cm),(*(Cm+no)),(*(Cm+2*no)),(*(Cm+3*no)));
                    Cm+=4*no;
                }
                for (int j=end; j<no; ++j) {
                    fprintf(out,"%18.10le ",*Cm);
                    Cm+=no;
                }
                fprintf(out,"\n");
            }
        }
        fclose(out);
        AnalyzeMoments(Pmat,basis.center_ptr(),ncen,no2);
        AnalyzeMOMoments(Cmat,no2,no);
    }


    void OintsOutput() {
        FILE* out=create_file("intsout.out");
        fprintf(out,"  UnoMol alpha version \n");
        fprintf(out," If you are using this for serious research you are insane!\n");
        fprintf(out,"\n");
        fprintf(stderr,"Using ");
        switch (scf_accel) {
        case 0:
            fprintf(stderr,"Simple Averaging scf\n");
            break;
        case 1:
            fprintf(stderr,"Anderson Averaging scf\n");
            break;
        }
        fprintf(out,"Number of Shells   = %12u\n",nshell);
        fprintf(out,"Number of Orbitals = %12u\n",no);
        fprintf(out,"Number of Centers  = %12u\n",ncen);
        fprintf(out,"Max L value        = %12u\n",basis.maxLvalue());
        fprintf(out,"Max SCF Iterations = %12u\n",maxits);
        fprintf(out,"Number of Electrons= %12u\n",nocc*2);
        fprintf(out," Centers \n");
        fprintf(out,"Charge     Position(x,y,z  in bohr)\n");
        const Center *center=basis.center_ptr();
        const Shell* shell=basis.shell_ptr();
        for (int i=0; i<ncen; ++i) {
            const double* rx;
            double qx=(center+i)->charge();
            rx=(center+i)->r_vec();
            fprintf(out,"%15.10lf %15.10lf %15.10lf %15.10lf\n",
                    qx,rx[0],rx[1],rx[2]);
        }
        fprintf(out,"\n");
        fprintf(out," Shells\n");
        for (int i=0; i<nshell; ++i) {
            int npr=(shell+i)->number_of_prims();
            int lsh=(shell+i)->Lvalue();
            int cen=(shell+i)->center();
            const double* al=(shell+i)->alf_ptr();
            const double* co=(shell+i)->cof_ptr();
            fprintf(out,"%3d %3d %3d\n",npr,lsh,cen);
            for (int j=0; j<npr; ++j)
                fprintf(out," %20.10lf %20.10lf\n",al[j],co[j]);
        }
        fprintf(out,"\n");
        fprintf(out,"Overlap Matrix \n");
        int ij=0;
        for (int i=0; i<no; ++i) {
            for (int j=0; j<=i; ++j) {
                fprintf(out,"%7u %7u %24.16le\n",i,j,Smat[ij]);
                ++ij;
            }
        }
        fprintf(out,"\n");
        fprintf(out,"Kinetic Energy Matrix \n");
        ij=0;
        for (int i=0; i<no; ++i) {
            for (int j=0; j<=i; ++j) {
                fprintf(out,"%7u %7u %24.16le\n",i,j,Tmat[ij]);
                ++ij;
            }
        }
        fprintf(out,"\n");
        fprintf(out,"Core Hamiltonian Matrix \n");
        ij=0;
        for (int i=0; i<no; ++i) {
            for (int j=0; j<=i; ++j) {
                fprintf(out,"%7u %7u %24.16le\n",i,j,Hmat[ij]);
                ++ij;
            }
        }
        fclose(out);
    }


    void scf_converger() {
        int i;
        double p00,p11,p01,beta;
        if (scf_accel==1) {
            if (extrap) {
                p00=p11=p01=0.0;
                for (i=0; i<no2; ++i) {
                    double s0=Pmat[i]-Pold[i];
                    double s1=Pold[i]-Pold2[i];
                    p00+=s0*s0;
                    p11+=s1*s1;
                    p01+=s0*s1;
                }
                beta=(p00-p01)/(p00-2.0*p01+p11);
                if (beta<0.001) beta=0.001;
                if (beta>1.500) beta=1.5;
                double betam1=1.0-beta;
                for (i=0; i<no2; ++i) {
                    Pold2[i]=Pold[i];
                    Pold[i]=Pmat[i];
                    Pmat[i]=Pold[i]*beta+Pold2[i]*betam1;
                }
                return;
            }
            extrap=1;
            for (i=0; i<no2; ++i) {
                Pold2[i]=Pold[i];
                Pold[i]=Pmat[i];
                Pmat[i]=(Pold2[i]+Pold[i])*0.5;
            }
            return;
        } else {
            for (int i=0; i<no2; ++i)
                Pmat[i]=(Pmat[i]+Pold[i])*0.5;
            return;
        }
    }

    void PopulationAnalysis(FILE *out) {
        const Center *center=basis.center_ptr();
        int ncen=basis.number_of_centers();
        double *Overlap=new double[no*no];
        double *NetChg=new double[ncen];
        for (int i=0; i<no; ++i) {
            int ii=i*(i+1)/2;
            for (int j=0; j<=i; ++j) {
                int ik=ii;
                int jk=j*(j+1)/2;
                double sum=0.0;
                for (int k=0; k<=j; ++ik,++jk,++k) {
                    sum+=Pmat[ik]*Smat[jk];
                }
                --jk;
                for (int k=j+1; k<=i; ++ik,++k) {
                    jk=jk+k;
                    sum+=Pmat[ik]*Smat[jk];
                }
                --ik;
                for (int k=i+1; k<no; ++k) {
                    jk=jk+k;
                    ik=ik+k;
                    sum+=Pmat[ik]*Smat[jk];
                }
                Overlap[i*no+j]=(sum+sum);
            }
            for (int j=i+1; j<no; ++j) {
                int ik=ii;
                int jk=j*(j+1)/2;
                double sum=0.0;
                for (int k=0; k<=i; ++ik,++jk,++k) {
                    sum+=Pmat[ik]*Smat[jk];
                }
                --ik;
                for (int k=i+1; k<=j; ++jk,++k) {
                    ik=ik+k;
                    sum+=Pmat[ik]*Smat[jk];
                }
                --jk;
                for (int k=j+1; k<no; ++k) {
                    ik=ik+k;
                    jk=jk+k;
                    sum+=Pmat[ik]*Smat[jk];
                }
                Overlap[i*no+j]=(sum+sum);
            }
        }
        for (int i=0; i<ncen; ++i) {
            NetChg[i]=(center+i)->charge();
        }
        const Shell* shell=basis.shell_ptr();
        const AuxFunctions& aux(*basis.auxfun_ptr());
        int ir=0;
        for (int i=0; i<nshell; ++i) {
            int lv=(shell+i)->Lvalue();
            int cn=(shell+i)->center();
            int nls=aux.number_of_lstates(lv);
            for (int ils=0; ils<nls; ++ils) {
                NetChg[cn]-=Overlap[ir*no+ir];
                ++ir;
            }
        }
        fprintf(out,"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
        fprintf(out,"           Mulliken Populations \n");
        fprintf(out,"Orbital  Net Population\n");
        for (int i=0; i<no; ++i) {
            fprintf(out,"%7u %25.16le\n",i+1,Overlap[i*no+i]);
        }
        fprintf(out,"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
        fprintf(out,"Mulliken Atomic Charges \n");
        fprintf(out,"Orbital   Nuclear Charge   Net Charge \n");
        for (int i=0; i<ncen; ++i) {
            double qc=(center+i)->charge();
            fprintf(out,"%7u %15.10lf %15.10lf \n",i+1,qc,NetChg[i]);
        }
        if (basis.prt_flags(2)) {
            FILE* pout=create_file("overlaps.out");
            fprintf(pout,"   MULLIKEN OVERLAP POPULATIONS\n");
            fprintf(pout,"  i      j                   value \n");
            for (int i=0; i<no; ++i) {
                for (int j=0; j<=i; ++j) {
                    fprintf(pout,"%5d %5d  %24.16le\n",i,j,Overlap[i*no+j]);
                }
            }
            fclose(pout);
        }
        delete [] NetChg;
        delete [] Overlap;
    }
  private:
    Basis& basis;
    TwoElectronInts& tints;
    double eps,ediff,pdiff,eold,nucrep,energy;
    double energyGs;
    double *PmatGs;
    double *Pold2;
    double *Pold;
    double *Pmat;
    double *Gmat;
    double *Hmat;
    double *Fock;
    double *Xmat;
    double *Wmat;
    double *Cmat;
    double *Tmat;
    double *Smat;
    double *Evals;
    double *sEvals;
    double *Wrka;
    int maxits,iteration,extrap;
    int no,no2,ncen,nshell,nocc;
    int tno,tno2,tncen,tnshell,scf_accel;
    int cflag;
};
}
#endif
