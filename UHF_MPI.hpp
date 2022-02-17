#ifndef _UHF_MPI_CC_
#define _UHF_MPI_CC_
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <mpi.h>
#include "Util.hpp"
#include "Basis.hpp"
#include "TwoElectronInts.hpp"
#include "AuxFunctions.hpp"
#include "GDPMInts.hpp"
#include "Moments.hpp"
#include "OneElectronInts.hpp"
#include "SymmPack.hpp"
#include "FField.hpp"
using namespace std;

namespace unomol
{

class UnRestrictedHartreeFockMPI
    {
    public:
        UnRestrictedHartreeFockMPI(Basis* b,const TwoElectronInts* x);
        ~UnRestrictedHartreeFockMPI();
        void update();
        void dpm_update(const TwoElectronInts& xints);
        void findEnergy();
        void formCmatrix(double* c);
        void formPmatrix(double* p,double *c,int noc);
        void PmatrixGuess();
        void formXmatrix();
        void FiniteFieldAnalysis();
        double nuclear_repulsion_energy(int ncen,const Center* center);
        void findPolarizationPotential();
        bool is_converged();
        void final_output(double init_energy);
        void OintsOutput();
        void scf_converger();
        void PopulationAnalysis(FILE* fp);
    private:
        Basis& basis;
        const TwoElectronInts& tints;
        double eps,ediff,pdiff,eold,nucrep,energy;
        double energyGs;
        double *PmatGsA;
        double *Pold2A;
        double *PoldA;
        double *PmatA;
        double *GmatA;
        double *PmatGsB;
        double *Pold2B;
        double *PoldB;
        double *PmatB;
        double *GmatB;
        double *Hmat;
        double *Fock;
        double *Xmat;
        double *Wmat;
        double *CmatA;
        double *CmatB;
        double *Tmat;
        double *Smat;
        double *EvalsA;
        double *EvalsB;
        double *sEvals;
        double *Wrka;
        double *GbufA;
        double *GbufB;
        int maxits,iteration,extrap;
        int no,no2,ncen,nshell,noccA,noccB;
        int tno,tno2,tncen,tnshell,scf_accel;
        int cflag,rank,im_done;
    };

UnRestrictedHartreeFockMPI::UnRestrictedHartreeFockMPI(Basis* b,
        const TwoElectronInts* t):basis(*b),tints(*t)
    {
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    nshell=basis.number_of_shells();
    no=basis.number_of_orbitals();
    no2=no*(no+1)/2;
    ncen=basis.number_of_centers();
    tnshell=basis.total_number_of_shells();
    tno=basis.total_number_of_orbitals();
    tno2=tno*(tno+1)/2;
    tncen=basis.total_number_of_centers();
    noccB=basis.number_of_electrons();
    noccA=noccB-noccB/2;
    noccB/=2;
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
    PmatA=new double[tno2];
    GmatA=new double[tno2];
    PmatB=new double[tno2];
    GmatB=new double[tno2];
    GbufA=new double[tno2];
    GbufB=new double[tno2];
    if (!rank)
        {
        PmatGsA=new double[tno2];
        PmatGsB=new double[tno2];
        Pold2A=new double[tno2];
        PoldA=new double[tno2];
        Pold2B=new double[tno2];
        PoldB=new double[tno2];
        Hmat=new double[tno2];
        Fock=new double[tno2];
        Xmat=new double[tnsqr];
        Wmat=new double[tnsqr];
        CmatA=new double[tnsqr];
        CmatB=new double[tnsqr];
        Tmat=new double[tno2];
        Smat=new double[tno2];
        EvalsA=new double[tno];
        EvalsB=new double[tno];
        sEvals=new double[tno];
        Wrka=new double[tnsqr];
        }
    }

UnRestrictedHartreeFockMPI::~UnRestrictedHartreeFockMPI()
    {
    if (!rank)
        {
        delete [] Wrka;
        delete [] sEvals;
        delete [] EvalsB;
        delete [] EvalsA;
        delete [] Smat;
        delete [] Tmat;
        delete [] CmatB;
        delete [] CmatA;
        delete [] Wmat;
        delete [] Xmat;
        delete [] Fock;
        delete [] Hmat;
        delete [] PoldB;
        delete [] Pold2B;
        delete [] PmatGsB;
        delete [] PoldA;
        delete [] Pold2A;
        delete [] PmatGsA;
        }
    delete [] GbufB;
    delete [] GbufA;
    delete [] GmatB;
    delete [] PmatB;
    delete [] GmatA;
    delete [] PmatA;

    }

inline void UnRestrictedHartreeFockMPI::update()
    {
    register int i;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(PmatA,no2,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(PmatB,no2,MPI_DOUBLE,0,MPI_COMM_WORLD);
    for (i=0; i<no2; ++i)
        {
        GbufA[i]=0.0;
        GbufB[i]=0.0;
        }
    tints.formGmatrix(PmatA,PmatB,GmatA,GmatB);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(GbufA,GmatA,no2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(GbufB,GmatB,no2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    if (rank)
        {
        ++iteration;
        return;
        }
    energy=SymmPack::TraceSymmPackProduct(PmatA,Hmat,no)+
           SymmPack::TraceSymmPackProduct(PmatB,Hmat,no)+
           0.5*(SymmPack::TraceSymmPackProduct(PmatA,GmatA,no)+
                SymmPack::TraceSymmPackProduct(PmatB,GmatB,no));
    ediff=energy-eold;
    eold=energy;
    for (i=0; i<no2; ++i) Fock[i]=Hmat[i]+GmatA[i];
    SymmPack::sp_trans(no,Fock,Xmat,Wrka);
    SymmPack::rsp(no,Fock,Wmat,EvalsA,Wrka);
    formCmatrix(CmatA);
    if (scf_accel==1) memcpy(Pold2A,PoldA,sizeof(double)*no2);
    memcpy(PoldA,PmatA,sizeof(double)*no2);
    formPmatrix(PmatA,CmatA,noccA);
    pdiff=SymmPack::SymmPackDiffNorm(PmatA,PoldA,no);
    for (i=0; i<no2; ++i) Fock[i]=Hmat[i]+GmatB[i];
    SymmPack::sp_trans(no,Fock,Xmat,Wrka);
    SymmPack::rsp(no,Fock,Wmat,EvalsB,Wrka);
    formCmatrix(CmatB);
    if (scf_accel==1) memcpy(Pold2B,PoldB,sizeof(double)*no2);
    memcpy(PoldB,PmatB,sizeof(double)*no2);
    formPmatrix(PmatB,CmatB,noccB);
    pdiff+=SymmPack::SymmPackDiffNorm(PmatB,PoldB,no);
    ++iteration;
    }

inline void UnRestrictedHartreeFockMPI::dpm_update(const TwoElectronInts& xints)
    {
    register int i;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(PmatA,no2,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(PmatB,no2,MPI_DOUBLE,0,MPI_COMM_WORLD);
    for (i=0; i<no2; ++i)
        {
        GbufA[i]=0.0;
        GbufB[i]=0.0;
        }
    tints.formGmatrix(PmatA,PmatB,GmatA,GmatB);
    xints.formGmatrix(PmatA,PmatB,GmatA,GmatB);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(GbufA,GmatA,no2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(GbufB,GmatB,no2,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    if (rank)
        {
        ++iteration;
        return;
        }
    energy=SymmPack::TraceSymmPackProduct(PmatA,Hmat,no)+
           SymmPack::TraceSymmPackProduct(PmatB,Hmat,no)+
           0.5*(SymmPack::TraceSymmPackProduct(PmatA,GmatA,no)+
                SymmPack::TraceSymmPackProduct(PmatB,GmatB,no));
    ediff=energy-eold;
    eold=energy;
    for (i=0; i<no2; ++i) Fock[i]=Hmat[i]+GmatA[i];
    SymmPack::sp_trans(no,Fock,Xmat,Wrka);
    SymmPack::rsp(no,Fock,Wmat,EvalsA,Wrka);
    formCmatrix(CmatA);
    if (scf_accel==1) memcpy(Pold2A,PoldA,sizeof(double)*no2);
    memcpy(PoldA,PmatA,sizeof(double)*no2);
    formPmatrix(PmatA,CmatA,noccA);
    pdiff=SymmPack::SymmPackDiffNorm(PmatA,PoldA,no);
    for (i=0; i<no2; ++i) Fock[i]=Hmat[i]+GmatB[i];
    SymmPack::sp_trans(no,Fock,Xmat,Wrka);
    SymmPack::rsp(no,Fock,Wmat,EvalsB,Wrka);
    formCmatrix(CmatB);
    if (scf_accel==1) memcpy(Pold2B,PoldB,sizeof(double)*no2);
    memcpy(PoldB,PmatB,sizeof(double)*no2);
    formPmatrix(PmatB,CmatB,noccB);
    pdiff+=SymmPack::SymmPackDiffNorm(PmatB,PoldB,no);
    ++iteration;
    }

inline void UnRestrictedHartreeFockMPI::findEnergy()
    {
    if (!rank)
        {
        nucrep=nuclear_repulsion_energy(ncen,basis.center_ptr());
        OneElectronInts(basis,Smat,Tmat,Hmat);
        MomentInts(basis);
        OintsOutput();
        formXmatrix();
        if (basis.scf_flags(2))
            {
            FILE* in=open_file("PMATRIX.DAT");
            fread(PmatA,sizeof(double),no2,in);
            fread(PmatB,sizeof(double),no2,in);
            fclose(in);
            }
        else
            {
            PmatrixGuess();
            }
        }
    im_done=0;
    iteration=0;
    extrap=0;
    eold=0.0;
    update();
    double init_energy=energy+nucrep;
    while (iteration<maxits)
        {
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&im_done,1,MPI_INT,0,MPI_COMM_WORLD);
        if (im_done) break;
        if (!rank) scf_converger();
        update();
        if (!rank)
            {
            if (is_converged()) im_done=1;
            }
        fprintf(stderr,"Iteration    =     %5d\n",iteration);
        fprintf(stderr,"Delta Energy =  %25.15le\n",ediff);
        }
    if (rank) return;
    memcpy(PmatGsA,PmatA,sizeof(double)*no2);
    for (register int j=no2; j<tno2; ++j) PmatGsA[j]=0.0;
    memcpy(PmatGsB,PmatB,sizeof(double)*no2);
    for (register int j=no2; j<tno2; ++j) PmatGsB[j]=0.0;
    energyGs=energy+nucrep;
    FILE* fp=create_file("PMATRIX.DAT");
    fwrite(PmatA,sizeof(double),no2,fp);
    fwrite(PmatB,sizeof(double),no2,fp);
    fclose(fp);
    final_output(init_energy);
    }

inline void UnRestrictedHartreeFockMPI::formCmatrix(double* c)
    {
    copy_trans(no,Wmat);
    for (register int i=0; i<no; ++i)
        {
        const double *Xi=Xmat+i*no;
        for (register int j=0; j<no; ++j)
            {
            const double *xp=Xi;
            const double *wp=Wmat+j*no;
            register double sum=0.0;
            for (register int k=0; k<no; ++k) sum+=xp[k]*wp[k];
            (*(c+i*no+j))=sum;
            }
        }
    }

inline void UnRestrictedHartreeFockMPI::formPmatrix(double* p,double *c,int noc)
    {
    register int ij=0;
    for (register int i=0; i<no; ++i)
        {
        const double *ci=c+i*no;
        for (register int j=0; j<=i; ++j,++ij)
            {
            const double *cj=c+j*no;
            register double sum=0.0;
            for (register int k=0; k<noc; ++k) sum+=ci[k]*cj[k];
            p[ij]=sum;
            }
        }
    }

inline void UnRestrictedHartreeFockMPI::PmatrixGuess()
    {
    for (register int i=0; i<no2; ++i) Fock[i]=Hmat[i];
    SymmPack::sp_trans(no,Fock,Xmat,Wrka);
    SymmPack::rsp(no,Fock,Wmat,EvalsA,Wrka);
    formCmatrix(CmatA);
    formPmatrix(PmatA,CmatA,noccA);
    int nsqr=no*no;
    memcpy(CmatB,CmatA,sizeof(double)*nsqr);
    formPmatrix(PmatB,CmatB,noccB);
    }


inline void UnRestrictedHartreeFockMPI::formXmatrix()
    {
    memcpy(Fock,Smat,sizeof(double)*no2);
    SymmPack::rsp(no,Fock,Wmat,sEvals,Wrka);
    for (int i=0; i<no; ++i)
        {
        double f=sEvals[i];
        if (fabs(f)<1.e-7)
            {
            if (f<0.0) f=-f;
            }
        if (f==0.0) fatal_error("Zero Eigenvalue in Overlap Matrix");
        if (f<0.0) fatal_error("Negative Eigenvalue in Overlap Matrix");
        sEvals[i]=1.0/sqrt(f);
        }
    for (int i=0; i<no; ++i)
        {
        double *Xi=Xmat+i*no;
        const double *Wi=Wmat+i*no;
        for (int j=0; j<no; ++j)
            {
            Xi[j]=Wi[j]*sEvals[j];
            }
        }
    }

inline void UnRestrictedHartreeFockMPI::FiniteFieldAnalysis()
    {
    double polnrg[3],alfpol[3];
    double Emag=basis.FiniteFieldValue();
    double Efield[3];
    for (int ix=0; ix<3; ++ix)
        {
        if (!rank)
            {
            Efield[2]=Efield[1]=Efield[0]=0.0;
            Efield[ix]=Emag;
            OneElectronInts(basis,Smat,Tmat,Hmat);
            FiniteFieldMatrix(basis,Hmat,Efield);
            formXmatrix();
            memcpy(PmatA,PmatGsA,sizeof(double)*no2);
            memcpy(PmatB,PmatGsB,sizeof(double)*no2);
            }
        im_done=0;
        iteration=0;
        extrap=0;
        eold=0.0;
        update();
        double init_energy=energy+nucrep;
        while (iteration<maxits)
            {
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(&im_done,1,MPI_INT,0,MPI_COMM_WORLD);
            if (im_done) break;
            if (!rank) scf_converger();
            update();
            if (!rank)
                {
                if (is_converged()) im_done=1;
                }
            }
        if (rank) continue;
        polnrg[ix]=energy+nucrep-energyGs;
        alfpol[ix]= -polnrg[ix]*2.0/Emag/Emag;
        }
    if (rank) return;
    FILE *out=create_file("finitefield.out");
    fprintf(out,"        Unomol Finite Field Analysis\n");
    fprintf(out,"        Electric field magnitude = %20.12le\n",Emag);
    fprintf(out,"     alfa                pol energy\n");
    fprintf(out," x   %20.12le %20.12le\n",alfpol[0],polnrg[0]);
    fprintf(out," y   %20.12le %20.12le\n",alfpol[1],polnrg[1]);
    fprintf(out," z   %20.12le %20.12le\n",alfpol[2],polnrg[2]);
    fclose(out);
    }

inline double UnRestrictedHartreeFockMPI::nuclear_repulsion_energy(int ncen,
        const Center* center)
    {
    register double sum=0.0;
    double q1,q2,r12;
    const double *r1,*r2;

    for (register int i=0; i<ncen; ++i)
        {
        q1=(center+i)->charge();
        r1=(center+i)->r_vec();
        for (register int j=i+1; j<ncen; ++j)
            {
            q2=(center+j)->charge();
            r2=(center+j)->r_vec();
            r12=sqrt(dist_sqr(r1,r2));
            sum+=q1*q2/r12;
            }
        }
    return sum;
    }

inline void UnRestrictedHartreeFockMPI::findPolarizationPotential()
    {
    double px,py,pz;
    const string xstr("XINTS.DAT");
    int sshell=nshell;
    nshell=tnshell;
    no=tno;
    no2=tno2;
    ncen=tncen;
    int pcen=basis.skip_center();
    int npts=basis.number_of_dpm_points();
    basis.dpm_augment();
    eps=1.e-12;
    ////////////////////////////////////////////////////////
    // start calculation
    FILE* in=open_file("pos.grid.dat");
    fscanf(in,"%15lf%15lf%15lf",&px,&py,&pz);
    basis.SetCenterPosition(px,py,pz,pcen);
    TwoElectronInts xints(basis,sshell,xstr);
    if (!rank)
        {
        nucrep=nuclear_repulsion_energy(ncen,
                                        basis.center_ptr());
        OneElectronInts(basis,Smat,Tmat,Hmat);
        GDPMInts(basis,Hmat);
        formXmatrix();
        memcpy(PmatA,PmatGsA,sizeof(double)*no2);
        memcpy(PmatB,PmatGsB,sizeof(double)*no2);
        }
    im_done=0;
    iteration=0;
    extrap=0;
    eold=0.0;
    dpm_update(xints);
    double init_energy=energy+nucrep;
    while (iteration<maxits)
        {
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&im_done,1,MPI_INT,0,MPI_COMM_WORLD);
        if (im_done) break;
        if (!rank) scf_converger();
        dpm_update(xints);
        if (!rank)
            if (is_converged()) im_done=1;
        }
    FILE *vout;
    FILE *sout;
    if (!rank)
        {
        double final_energy=energy+nucrep;
        double vpol=final_energy-init_energy;
        double vstat=init_energy-energyGs+nucrep;
        double r2=px*px+py*py+pz*pz;
        double alfa= -2.0*vpol*r2*r2;
        vout = create_file("vpol.out");
        sout = create_file("spol.out");
        fprintf(vout,"%15.10lf%15.10lf%15.10lf %25.15le %25.15le\n",px,py,pz,
                vpol,alfa);
        fflush(vout);
        fprintf(sout,"%3d %20.10le %20.10le %20.10le %20.10le %25.15le\n",
                0,energyGs,init_energy,final_energy,vstat,ediff);
        fflush(sout);
        }
    for (int i=1; i<npts; ++i)
        {
        fscanf(in,"%15lf%15lf%15lf",&px,&py,&pz);
        basis.SetCenterPosition(px,py,pz,pcen);
        xints.recalculate(basis);
        if (!rank)
            {
            nucrep=nuclear_repulsion_energy(basis.number_of_centers(),
                                            basis.center_ptr());
            OneElectronInts(basis,Smat,Tmat,Hmat);
            GDPMInts(basis,Hmat);
            formXmatrix();
            memcpy(PmatA,PmatGsA,sizeof(double)*no2);
            memcpy(PmatB,PmatGsB,sizeof(double)*no2);
            }
        im_done=0;
        iteration=0;
        extrap=0;
        eold=0.0;
        dpm_update(xints);
        double init_energy=energy+nucrep;
        while (iteration<maxits)
            {
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(&im_done,1,MPI_INT,0,MPI_COMM_WORLD);
            if (im_done) break;
            if (!rank) scf_converger();
            dpm_update(xints);
            if (!rank)
                if (is_converged()) im_done=1;
            }
        if (rank) continue;
        double final_energy=energy+nucrep;
        double vpol=final_energy-init_energy;
        double vstat=init_energy-energyGs;
        double r2=px*px+py*py+pz*pz;
        double alfa= -2.0*vpol*(r2*r2);
        fprintf(vout,"%15.10lf%15.10lf%15.10lf %25.15le %25.15le\n",px,py,pz,
                vpol,alfa);
        fflush(vout);
        fprintf(sout,"%3d %20.10le %20.10le %20.10le %20.10le %25.15le\n",
                i,energyGs,init_energy,final_energy,vstat,ediff);
        fflush(sout);
        }
    if (rank) return;
    fclose(vout);
    fclose(sout);
    }

inline bool UnRestrictedHartreeFockMPI::is_converged()
    {
    switch (cflag)
        {
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

inline void UnRestrictedHartreeFockMPI::final_output(double init_energy)
    {
    FILE *out=create_file("short.gs.dat");
    fprintf(out,"%25.16le\n",init_energy);
    fprintf(out,"%25.16le\n",energy+nucrep);
    fprintf(out,"%25.16le\n",ediff);
    fclose(out);
    out=create_file("scfout.gs.dat");
    double trace_t=SymmPack::TraceSymmPackProduct(PmatA,Tmat,no)+
                   SymmPack::TraceSymmPackProduct(PmatB,Tmat,no);
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
    fprintf(out,"        Alpha Orbital Energies\n");
    fprintf(out,"Orbital Energy          Occupancy\n");
    for (int i=0; i<noccA; i++)
        {
        fprintf(out,"%7u %25.16le %12u\n",i+1,EvalsA[i],1);
        }
    for (int i=noccA; i<no; i++)
        {
        fprintf(out,"%7u %25.16le %12u\n",i+1,EvalsA[i],0);
        }
    fprintf(out,"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
    fprintf(out,"        Beta Orbital Energies\n");
    fprintf(out,"Orbital Energy          Occupancy\n");
    for (int i=0; i<noccB; i++)
        {
        fprintf(out,"%7u %25.16le %12u\n",i+1,EvalsB[i],1);
        }
    for (int i=noccB; i<no; i++)
        {
        fprintf(out,"%7u %25.16le %12u\n",i+1,EvalsB[i],0);
        }
    fprintf(out,"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
    PopulationAnalysis(out);
    if (basis.prt_flags(0))
        {
        fprintf(out,"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
        fprintf(out,"              Alpha Fock Eigenvectors\n");
        fprintf(out,"\n");
        for (int i=0; i<no; ++i)
            {
            int occup=(i<noccA)?1:0;
            fprintf(out,"Column %5d Occupancy= %2d EigenEnergy= %25.15le\n",i,
                    occup,EvalsA[i]);
            int nex=no-4*no/4;
            int end=4*(no/4);
            const double *Cm=CmatA+i;
            for (int j=0; j<no; j+=4)
                {
                fprintf(out,"%18.10le %18.10le %18.10le %18.10le\n",
                        (*Cm),(*(Cm+no)),(*(Cm+2*no)),(*(Cm+3*no)));
                Cm+=4*no;
                }
            for (int j=end; j<no; ++j)
                {
                fprintf(out,"%18.10le ",*Cm);
                Cm+=no;
                }
            }
        fprintf(out,"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
        fprintf(out,"              Beta Fock Eigenvectors\n");
        fprintf(out,"\n");
        for (int i=0; i<no; ++i)
            {
            int occup=(i<noccB)?1:0;
            fprintf(out,"Column %5d Occupancy= %2d EigenEnergy= %25.15le\n",i,
                    occup,EvalsB[i]);
            int nex=no-4*no/4;
            int end=4*(no/4);
            const double *Cm=CmatB+i;
            for (int j=0; j<no; j+=4)
                {
                fprintf(out,"%18.10le %18.10le %18.10le %18.10le\n",
                        (*Cm),(*(Cm+no)),(*(Cm+2*no)),(*(Cm+3*no)));
                Cm+=4*no;
                }
            for (int j=end; j<no; ++j)
                {
                fprintf(out,"%18.10le ",*Cm);
                Cm+=no;
                }
            }
        }
    fclose(out);
    AnalyzeMoments(PmatA,PmatB,basis.center_ptr(),ncen,no2);
    };


inline void UnRestrictedHartreeFockMPI::OintsOutput()
    {
    FILE* out=create_file("intsout.dat");
    fprintf(out,"  UnoMol alpha version \n");
    fprintf(out," If you are using this for serious research you are insane!\n");
    fprintf(out,"\n");
    fprintf(stderr,"Using ");
    switch (scf_accel)
        {
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
    fprintf(out,"Number of Electrons= %12u\n",noccA+noccB);
    fprintf(out," Centers \n");
    fprintf(out,"Charge     Position(x,y,z  in bohr)\n");
    const Center *center=basis.center_ptr();
    const Shell* shell=basis.shell_ptr();
    for (int i=0; i<ncen; ++i)
        {
        const double* rx;
        double qx=(center+i)->charge();
        rx=(center+i)->r_vec();
        fprintf(out,"%15.10lf %15.10lf %15.10lf %15.10lf\n",
                qx,rx[0],rx[1],rx[2]);
        }
    fprintf(out,"\n");
    fprintf(out," Shells\n");
    for (int i=0; i<nshell; ++i)
        {
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
    for (int i=0; i<no; ++i)
        {
        for (int j=0; j<=i; ++j)
            {
            fprintf(out,"%7u %7u %24.16le\n",i,j,Smat[ij]);
            ++ij;
            }
        }
    fprintf(out,"\n");
    fprintf(out,"Kinetic Energy Matrix \n");
    ij=0;
    for (int i=0; i<no; ++i)
        {
        for (int j=0; j<=i; ++j)
            {
            fprintf(out,"%7u %7u %24.16le\n",i,j,Tmat[ij]);
            ++ij;
            }
        }
    fprintf(out,"\n");
    fprintf(out,"Core Hamiltonian Matrix \n");
    ij=0;
    for (int i=0; i<no; ++i)
        {
        for (int j=0; j<=i; ++j)
            {
            fprintf(out,"%7u %7u %24.16le\n",i,j,Hmat[ij]);
            ++ij;
            }
        }
    fclose(out);
    }


inline void
UnRestrictedHartreeFockMPI::scf_converger()
    {
    register int i;
    double p00,p11,p01,beta;
    if (scf_accel==1)
        {
        if (extrap)
            {
            p00=p11=p01=0.0;
            for (i=0; i<no2; ++i)
                {
                double s0=PmatA[i]-PoldA[i];
                double s1=PoldA[i]-Pold2A[i];
                p00+=s0*s0;
                p11+=s1*s1;
                p01+=s0*s1;
                }
            beta=(p00-p01)/(p00-2.0*p01+p11);
            if (beta<0.001) beta=0.001;
            if (beta>1.500) beta=1.5;
            double betam1=1.0-beta;
            for (i=0; i<no2; ++i)
                {
                Pold2A[i]=PoldA[i];
                PoldA[i]=PmatA[i];
                PmatA[i]=PoldA[i]*beta+Pold2A[i]*betam1;
                }
            p00=p11=p01=0.0;
            for (i=0; i<no2; ++i)
                {
                double s0=PmatB[i]-PoldB[i];
                double s1=PoldB[i]-Pold2B[i];
                p00+=s0*s0;
                p11+=s1*s1;
                p01+=s0*s1;
                }
            beta=(p00-p01)/(p00-2.0*p01+p11);
            if (beta<0.001) beta=0.001;
            if (beta>1.500) beta=1.5;
            betam1=1.0-beta;
            for (i=0; i<no2; ++i)
                {
                Pold2B[i]=PoldB[i];
                PoldB[i]=PmatB[i];
                PmatB[i]=PoldB[i]*beta+Pold2B[i]*betam1;
                }
            return;
            }
        extrap=1;
        for (i=0; i<no2; ++i)
            {
            Pold2A[i]=PoldA[i];
            PoldA[i]=PmatA[i];
            PmatA[i]=(Pold2A[i]+PoldA[i])*0.5;
            }
        for (i=0; i<no2; ++i)
            {
            Pold2B[i]=PoldB[i];
            PoldB[i]=PmatB[i];
            PmatB[i]=(Pold2B[i]+PoldB[i])*0.5;
            }
        return;
        }
    else
        {
        for (register int i=0; i<no2; ++i)
            {
            PmatA[i]=(PmatA[i]+PoldA[i])*0.5;
            PmatB[i]=(PmatB[i]+PoldB[i])*0.5;
            }
        return;
        }
    }

inline void UnRestrictedHartreeFockMPI::PopulationAnalysis(FILE *out)
    {
    const Center *center=basis.center_ptr();
    int ncen=basis.number_of_centers();
    double *OverlapA=new double[no*no];
    double *OverlapB=new double[no*no];
    double *NetChg=new double[ncen];
    for (int i=0; i<no; ++i)
        {
        int ii=i*(i+1)/2;
        for (int j=0; j<=i; ++j)
            {
            int ik=ii;
            int jk=j*(j+1)/2;
            double sum=0.0;
            for (int k=0; k<=j; ++ik,++jk,++k)
                {
                sum+=PmatA[ik]*Smat[jk];
                }
            --jk;
            for (int k=j+1; k<=i; ++ik,++k)
                {
                jk=jk+k;
                sum+=PmatA[ik]*Smat[jk];
                }
            --ik;
            for (int k=i+1; k<no; ++k)
                {
                jk=jk+k;
                ik=ik+k;
                sum+=PmatA[ik]*Smat[jk];
                }
            OverlapA[i*no+j]=(sum);
            }
        for (int j=i+1; j<no; ++j)
            {
            int ik=ii;
            int jk=j*(j+1)/2;
            double sum=0.0;
            for (int k=0; k<=i; ++ik,++jk,++k)
                {
                sum+=PmatA[ik]*Smat[jk];
                }
            --ik;
            for (int k=i+1; k<=j; ++jk,++k)
                {
                ik=ik+k;
                sum+=PmatA[ik]*Smat[jk];
                }
            --jk;
            for (int k=j+1; k<no; ++k)
                {
                ik=ik+k;
                jk=jk+k;
                sum+=PmatA[ik]*Smat[jk];
                }
            OverlapA[i*no+j]=(sum);
            }
        }
    for (int i=0; i<no; ++i)
        {
        int ii=i*(i+1)/2;
        for (int j=0; j<=i; ++j)
            {
            int ik=ii;
            int jk=j*(j+1)/2;
            double sum=0.0;
            for (int k=0; k<=j; ++ik,++jk,++k)
                {
                sum+=PmatB[ik]*Smat[jk];
                }
            --jk;
            for (int k=j+1; k<=i; ++ik,++k)
                {
                jk=jk+k;
                sum+=PmatB[ik]*Smat[jk];
                }
            --ik;
            for (int k=i+1; k<no; ++k)
                {
                jk=jk+k;
                ik=ik+k;
                sum+=PmatB[ik]*Smat[jk];
                }
            OverlapB[i*no+j]=(sum);
            }
        for (int j=i+1; j<no; ++j)
            {
            int ik=ii;
            int jk=j*(j+1)/2;
            double sum=0.0;
            for (int k=0; k<=i; ++ik,++jk,++k)
                {
                sum+=PmatB[ik]*Smat[jk];
                }
            --ik;
            for (int k=i+1; k<=j; ++jk,++k)
                {
                ik=ik+k;
                sum+=PmatB[ik]*Smat[jk];
                }
            --jk;
            for (int k=j+1; k<no; ++k)
                {
                ik=ik+k;
                jk=jk+k;
                sum+=PmatB[ik]*Smat[jk];
                }
            OverlapB[i*no+j]=(sum);
            }
        }
    for (int i=0; i<ncen; ++i)
        {
        NetChg[i]=(center+i)->charge();
        }
    const Shell* shell=basis.shell_ptr();
    const AuxFunctions& aux(*basis.auxfun_ptr());
    int ir=0;
    for (int i=0; i<nshell; ++i)
        {
        int lv=(shell+i)->Lvalue();
        int cn=(shell+i)->center();
        int nls=aux.number_of_lstates(lv);
        for (int ils=0; ils<nls; ++ils)
            {
            NetChg[cn]-=(OverlapA[ir*no+ir]+OverlapB[ir*no+ir]);
            ++ir;
            }
        }
    fprintf(out,"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
    fprintf(out,"           Mulliken Populations \n");
    fprintf(out,"Orbital  Alpha Population          Beta Population\n");
    for (int i=0; i<no; ++i)
        {
        fprintf(out,"%7u %25.16le %25.16le\n",i+1,OverlapA[i*no+i],OverlapB[i*no+i]);
        }
    fprintf(out,"xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
    fprintf(out,"Mulliken Atomic Charges \n");
    fprintf(out,"Orbital   Nuclear Charge   Net Charge \n");
    for (int i=0; i<ncen; ++i)
        {
        double qc=(center+i)->charge();
        fprintf(out,"%7u %15.10lf %15.10lf \n",i+1,qc,NetChg[i]);
        }
    if (basis.prt_flags(2))
        {
        FILE* pout=create_file("overlaps.dat");
        fprintf(pout,"   MULLIKEN OVERLAP POPULATIONS\n");
        fprintf(pout,"  i      j        Alpha value                 Beta Value\n");
        for (int i=0; i<no; ++i)
            {
            for (int j=0; j<=i; ++j)
                {
                fprintf(pout,"%5d %5d  %24.16le  %24.16le\n",i,j,
                        OverlapA[i*no+j],OverlapB[i*no+j]);
                }
            }
        fclose(pout);
        }
    delete [] NetChg;
    delete [] OverlapB;
    delete [] OverlapA;
    }

}
#endif

