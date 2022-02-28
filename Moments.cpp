#include "Moments.hpp"

namespace unomol {

void calc_moments (MomInts * mvals,
                   ShellPairData & sp,
                   const AuxFunctions & aux,
                   MD_Dfunction & dx, MD_Dfunction & dy, MD_Dfunction & dz) {
    int ipr, jpr, kc;
    int l1, m1, n1, l2, m2, n2, ls1, ls2;
    double c1, c12, axp, bxp, pxp, abi, p[3], s12, nfact, ovl, sov, pvx, pvy,
           pvz;
    const double *a, *b;
    const double piterm=5.5683279968317079;
    a = sp.a;
    b = sp.b;
    for (ipr = 0; ipr < sp.npr1; ++ipr) {
        axp = sp.al1[ipr];
        c1 = sp.co1[ipr];
        for (jpr = 0; jpr < sp.npr2; ++jpr) {
            c12 = c1 * sp.co2[jpr];
            bxp = sp.al2[jpr];
            pxp = axp + bxp;
            abi = 1.0 / pxp;
            s12 = piterm * exp (-axp * bxp * abi * sp.ab2) / (pxp * sqrt (pxp));
            p[0] = (axp * a[0] + bxp * b[0]) * abi;
            p[1] = (axp * a[1] + bxp * b[1]) * abi;
            p[2] = (axp * a[2] + bxp * b[2]) * abi;
            abi*=0.5;
            dx.eval (abi, p[0] - a[0], p[0] - b[0], sp.lv1 + 1, sp.lv2 + 1);
            dy.eval (abi, p[1] - a[1], p[1] - b[1], sp.lv1 + 1, sp.lv2 + 1);
            dz.eval (abi, p[2] - a[2], p[2] - b[2], sp.lv1 + 1, sp.lv2 + 1);
            for (kc = 0; kc < sp.len; ++kc) {
                unsigned short key=sp.lstates[kc];
                ls2= key&0xF;
                ls1= key>>4;
                l1 = aux.Lxyz (sp.lv1, ls1, 0);
                m1 = aux.Lxyz (sp.lv1, ls1, 1);
                n1 = aux.Lxyz (sp.lv1, ls1, 2);
                l2 = aux.Lxyz (sp.lv2, ls2, 0);
                m2 = aux.Lxyz (sp.lv2, ls2, 1);
                n2 = aux.Lxyz (sp.lv2, ls2, 2);
                nfact = aux.normalization_factor (sp.lv1, ls1) *
                        aux.normalization_factor (sp.lv2, ls2) * c12;
                sov = nfact * s12;
                ovl = sov * dx.getValue (l1, l2, 0) *
                      dy.getValue (m1, m2, 0) * dz.getValue (n1, n2, 0);
                pvx = sov * dy.getValue (m1, m2, 0) * dz.getValue (n1, n2, 0);
                pvy = sov * dx.getValue (l1, l2, 0) * dz.getValue (n1, n2, 0);
                pvz = sov * dy.getValue (m1, m2, 0) * dx.getValue (l1, l2, 0);
                //////////////////// dipole moments ////////////////////////////////
                (mvals + kc)->dx +=
                    (ovl * p[0] + dx.getValue (l1, l2, 1) * pvx);
                (mvals + kc)->dy +=
                    (ovl * p[1] + dy.getValue (m1, m2, 1) * pvy);
                (mvals + kc)->dz +=
                    (ovl * p[2] + dz.getValue (n1, n2, 1) * pvz);
                //////////////////// quadrupole moments ////////////////////////////
                (mvals + kc)->qxx +=
                    (ovl * (p[0] * p[0] + abi) +
                     2.0 * pvx * ( p[0]*dx.getValue(l1,l2,1)+
                                   dx.getValue(l1,l2,2) ) );
                (mvals + kc)->qxy +=
                    (sov * dz.getValue (n1, n2, 0) *
                     (dx.getValue (l1, l2, 1) * dy.getValue (m1, m2, 1) +
                      dx.getValue (l1, l2, 1) * dy.getValue (m1, m2, 0) * p[1] +
                      dx.getValue (l1, l2, 0) * dy.getValue (m1, m2, 1) * p[0])+
                     ovl * p[0] * p[1]);
                (mvals + kc)->qxz +=
                    (sov * dy.getValue (m1, m2, 0) *
                     (dx.getValue (l1, l2, 1) * dz.getValue (n1, n2, 1) +
                      dx.getValue (l1, l2, 1) * dz.getValue (n1, n2,0) * p[2] +
                      dx.getValue (l1, l2, 0) * dz.getValue (n1, n2,1) * p[0]) +
                     ovl * p[0] * p[2]);
                (mvals + kc)->qyy +=
                    (ovl * (p[1] * p[1] + abi) +
                     2.0 * pvy * (p[0] * dy.getValue (m1, m2, 1) +
                                  dy.getValue (m1, m2, 2)));
                (mvals + kc)->qyz +=
                    (sov * dx.getValue (l1, l2, 0) *
                     (dy.getValue (m1, m2, 1) * dz.getValue (n1, n2, 1) +
                      dy.getValue (m1, m2, 1) * dz.getValue (n1, n2,0) * p[2] +
                      dy.getValue (m1, m2, 0) * dz.getValue (n1, n2,1) * p[1]) +
                     ovl * p[1] * p[2]);
                (mvals + kc)->qzz +=
                    (ovl * (p[2] * p[2] + abi) +
                     2.0 * pvz * (p[2] * dz.getValue (n1, n2, 1) +
                                  dz.getValue (n1, n2, 2)));
            }
        }
    }
}


void MomentInts (const Basis & basis ) {
    int ish, icen, jsh, jcen;

    int nshell = basis.number_of_shells ();
    int lmax = basis.maxLvalue ();
    const Shell *shell = basis.shell_ptr ();
    const Center *center = basis.center_ptr ();
    const AuxFunctions& aux(*(basis.auxfun_ptr()));
    int maxlst = aux.number_of_lstates (lmax);
    int maxlst2 = maxlst * maxlst;
    ShellPairData sp;
    sp.lstates = new unsigned short[maxlst2];
    MomInts *mvals = new MomInts[maxlst2];
    double *factors = new double[maxlst2];
    MD_Dfunction dx (lmax + 2);
    MD_Dfunction dy (lmax + 2);
    MD_Dfunction dz (lmax + 2);
    FILE* out=create_file("MOMINTS.DAT");
    FILE* rout=create_file("RMOM.DAT");
    if (!out) {
        fatal_error("could not open MOMINTS.DAT for writing");
    }
    int ir0 = 0;
    int nls1 = 0;
    for (ish = 0; ish < nshell; ++ish) {
        ir0 = ir0 + nls1;
        sp.npr1 = (shell + ish)->number_of_prims ();
        sp.lv1 = (shell + ish)->Lvalue ();
        icen = (shell + ish)->center ();
        sp.al1 = (shell + ish)->alf_ptr ();
        sp.co1 = (shell + ish)->cof_ptr ();
        sp.a = (center + icen)->r_vec ();
        nls1 = aux.number_of_lstates (sp.lv1);
        int jr0 = 0;
        int nls2 = 0;
        for (jsh = 0; jsh <= ish; ++jsh) {
            jr0 = jr0 + nls2;
            sp.npr2 = (shell + jsh)->number_of_prims ();
            sp.lv2 = (shell + jsh)->Lvalue ();
            jcen = (shell + jsh)->center ();
            sp.al2 = (shell + jsh)->alf_ptr ();
            sp.co2 = (shell + jsh)->cof_ptr ();
            nls2 = aux.number_of_lstates (sp.lv2);
            sp.b = (center + jcen)->r_vec ();
            sp.ab2=dist_sqr(sp.a,sp.b);
            int knt = 0;
            int ir = ir0;
            for (int ils = 0; ils < nls1; ++ir, ++ils) {
                int ijr = ir * (ir + 1) / 2 + jr0;
                int jr = jr0;
                int jend = nls2;
                if (ish == jsh)
                    jend = ils + 1;
                for (int jls = 0; jls < jend; ++ijr, ++jr, ++jls) {
                    (mvals + knt)->dx = 0.0;
                    (mvals + knt)->dy = 0.0;
                    (mvals + knt)->dz = 0.0;
                    (mvals + knt)->qxx = 0.0;
                    (mvals + knt)->qxy = 0.0;
                    (mvals + knt)->qxz = 0.0;
                    (mvals + knt)->qyy = 0.0;
                    (mvals + knt)->qyz = 0.0;
                    (mvals + knt)->qzz = 0.0;
                    (mvals + knt)->ijr = ijr;
                    (sp.lstates)[knt] = (unsigned short)((ils<<4)+jls);
                    factors[knt] = 4.0;
                    if (ir == jr) {
                        factors[knt] = 2.0;
                    }
                    ++knt;
                }
            }
            sp.len = knt;
            if (!knt)
                continue;
            calc_moments (mvals, sp, aux, dx, dy, dz);
            Fwrite(mvals,sizeof(MomInts),knt,rout);
            for (int ils = 0; ils < sp.len; ++ils) {
                double fact = factors[ils];
                (((mvals + ils)->dx)) *= fact;
                ((mvals + ils)->dy) *= fact;
                ((mvals + ils)->dz) *= fact;
                ((mvals + ils)->qxx) *= fact;
                ((mvals + ils)->qxy) *= fact;
                ((mvals + ils)->qxz) *= fact;
                ((mvals + ils)->qyy) *= fact;
                ((mvals + ils)->qyz) *= fact;
                ((mvals + ils)->qzz) *= fact;
                Fwrite((mvals+ils),sizeof(MomInts),1,out);
            }
        }
    }
    fclose(rout);
    fclose(out);
    delete [] factors;
    delete [] mvals;
    delete [] sp.lstates;
}


void AnalyzeMoments (const double *Pmat, const Center * center, int ncen, int no2) {
    MomInts mint;
    MomInts eM, nM, tM;

    nM.dx = 0.0;
    nM.dy = 0.0;
    nM.dz = 0.0;
    nM.qxx = 0.0;
    nM.qxy = 0.0;
    nM.qxz = 0.0;
    nM.qyy = 0.0;
    nM.qyz = 0.0;
    nM.qzz = 0.0;
    for (int ic = 0; ic < ncen; ++ic) {
        double qi = (center + ic)->charge ();
        const double *ri = (center + ic)->r_vec ();
        nM.dx += qi * ri[0];
        nM.dy += qi * ri[1];
        nM.dz += qi * ri[2];
        nM.qxx += qi * ri[0] * ri[0];
        nM.qxy += qi * ri[0] * ri[1];
        nM.qxz += qi * ri[0] * ri[2];
        nM.qyy += qi * ri[1] * ri[1];
        nM.qyz += qi * ri[1] * ri[2];
        nM.qzz += qi * ri[2] * ri[2];
    }
    eM.dx = 0.0;
    eM.dy = 0.0;
    eM.dz = 0.0;
    eM.qxx = 0.0;
    eM.qxy = 0.0;
    eM.qxz = 0.0;
    eM.qyy = 0.0;
    eM.qyz = 0.0;
    eM.qzz = 0.0;
    FILE* in=open_file("MOMINTS.DAT");
    if (!in) {
        fatal_error("could not open MOMINTS.DAT for reading");
    }
    for (int i = 0; i < no2; ++i) {
        Fread(&mint,sizeof(MomInts),1,in);
        double pij = Pmat[mint.ijr];
        eM.dx += pij * (mint.dx);
        eM.dy += pij * (mint.dy);
        eM.dz += pij * (mint.dz);
        eM.qxx += pij * (mint.qxx);
        eM.qxy += pij * (mint.qxy);
        eM.qxz += pij * (mint.qxz);
        eM.qyy += pij * (mint.qyy);
        eM.qyz += pij * (mint.qyz);
        eM.qzz += pij * (mint.qzz);
    }
    fclose(in);
    tM.dx = nM.dx - eM.dx;
    tM.dy = nM.dy - eM.dy;
    tM.dz = nM.dz - eM.dz;
    tM.qxx = nM.qxx - eM.qxx;
    tM.qxy = nM.qxy - eM.qxy;
    tM.qxz = nM.qxz - eM.qxz;
    tM.qyy = nM.qyy - eM.qyy;
    tM.qyz = nM.qyz - eM.qyz;
    tM.qzz = nM.qzz - eM.qzz;
    double q00 = (tM.qzz - 0.5 * (tM.qxx + tM.qyy));
    FILE* out=create_file("moments.dat");
    fprintf(out," MULTIPOLE MOMENT ANALYSIS \n");
    fprintf(out," units in bohr - hartree atomic units \n\n");
    fprintf(out," DIPOLE MOMENTS \n\n");
    fprintf(out,"          Total         Electronic        Nuclear \n");
    fprintf(out," x  %15.7le  %15.7le  %15.7le \n", tM.dx, eM.dx, nM.dx);
    fprintf(out," y  %15.7le  %15.7le  %15.7le \n", tM.dy, eM.dy, nM.dy);
    fprintf(out," z  %15.7le  %15.7le  %15.7le \n", tM.dz, eM.dz, nM.dz);
    fprintf(out,"\n");
    fprintf(out," QUADRUPOLE MOMENTS \n\n");
    fprintf(out,"          Total         Electronic        Nuclear \n");
    fprintf(out," xx  %15.7le  %15.7le  %15.7le \n", tM.qxx, eM.qxx, nM.qxx);
    fprintf(out," xy  %15.7le  %15.7le  %15.7le \n", tM.qxy, eM.qxy, nM.qxy);
    fprintf(out," xz  %15.7le  %15.7le  %15.7le \n", tM.qxz, eM.qxz, nM.qxz);
    fprintf(out," yy  %15.7le  %15.7le  %15.7le \n", tM.qyy, eM.qyy, nM.qyy);
    fprintf(out," yz  %15.7le  %15.7le  %15.7le \n", tM.qyz, eM.qyz, nM.qyz);
    fprintf(out," zz  %15.7le  %15.7le  %15.7le \n", tM.qzz, eM.qzz, nM.qzz);
    fprintf(out,"\n");
    fprintf(out," Dipole moment     = %25.15le \n", tM.dz);
    fprintf(out," Quadrupole moment = %25.15le \n", q00);
    fclose ( out );
}


void AnalyzeMoments (const double *PmatA, const double* PmatB,
                     const Center * center, int ncen, int no2) {
    MomInts mint;
    MomInts eM, nM, tM;

    nM.dx = 0.0;
    nM.dy = 0.0;
    nM.dz = 0.0;
    nM.qxx = 0.0;
    nM.qxy = 0.0;
    nM.qxz = 0.0;
    nM.qyy = 0.0;
    nM.qyz = 0.0;
    nM.qzz = 0.0;
    for (int ic = 0; ic < ncen; ++ic) {
        double qi = (center + ic)->charge ();
        const double *ri = (center + ic)->r_vec ();
        nM.dx += qi * ri[0];
        nM.dy += qi * ri[1];
        nM.dz += qi * ri[2];
        nM.qxx += qi * ri[0] * ri[0];
        nM.qxy += qi * ri[0] * ri[1];
        nM.qxz += qi * ri[0] * ri[2];
        nM.qyy += qi * ri[1] * ri[1];
        nM.qyz += qi * ri[1] * ri[2];
        nM.qzz += qi * ri[2] * ri[2];
    }
    eM.dx = 0.0;
    eM.dy = 0.0;
    eM.dz = 0.0;
    eM.qxx = 0.0;
    eM.qxy = 0.0;
    eM.qxz = 0.0;
    eM.qyy = 0.0;
    eM.qyz = 0.0;
    eM.qzz = 0.0;
    FILE* in=open_file("MOMINTS.DAT");
    if (!in) {
        fatal_error("could not open MOMINTS.DAT for reading");
    }
    for (int i = 0; i < no2; ++i) {
        Fread(&mint,sizeof(MomInts),1,in);
        int ijr=mint.ijr;
        double pij = 0.5*(PmatA[ijr]+PmatB[ijr]);
        eM.dx += pij * (mint.dx);
        eM.dy += pij * (mint.dy);
        eM.dz += pij * (mint.dz);
        eM.qxx += pij * (mint.qxx);
        eM.qxy += pij * (mint.qxy);
        eM.qxz += pij * (mint.qxz);
        eM.qyy += pij * (mint.qyy);
        eM.qyz += pij * (mint.qyz);
        eM.qzz += pij * (mint.qzz);
    }
    fclose(in);
    tM.dx = nM.dx - eM.dx;
    tM.dy = nM.dy - eM.dy;
    tM.dz = nM.dz - eM.dz;
    tM.qxx = nM.qxx - eM.qxx;
    tM.qxy = nM.qxy - eM.qxy;
    tM.qxz = nM.qxz - eM.qxz;
    tM.qyy = nM.qyy - eM.qyy;
    tM.qyz = nM.qyz - eM.qyz;
    tM.qzz = nM.qzz - eM.qzz;
    double q00 = (tM.qzz - 0.5 * (tM.qxx + tM.qyy));
    FILE* out = create_file ("moments.dat");
    fprintf(out," MULTIPOLE MOMENT ANALYSIS \n");
    fprintf(out," units in bohr - hartree atomic units \n");
    fprintf(out,"\n");
    fprintf(out," DIPOLE MOMENTS \n\n");
    fprintf(out, "          Total         Electronic        Nuclear \n");
    fprintf(out," x  %15.7le  %15.7le  %15.7le \n", tM.dx, eM.dx, nM.dx);
    fprintf(out," y  %15.7le  %15.7le  %15.7le \n", tM.dy, eM.dy, nM.dy);
    fprintf(out," z  %15.7le  %15.7le  %15.7le \n", tM.dz, eM.dz, nM.dz);
    fprintf(out,"\n");
    fprintf(out," QUADRUPOLE MOMENTS \n \n");
    fprintf(out,"          Total         Electronic        Nuclear \n");
    fprintf(out," xx  %15.7le  %15.7le  %15.7le \n", tM.qxx, eM.qxx, nM.qxx);
    fprintf(out," xy  %15.7le  %15.7le  %15.7le \n", tM.qxy, eM.qxy, nM.qxy);
    fprintf(out," xz  %15.7le  %15.7le  %15.7le \n", tM.qxz, eM.qxz, nM.qxz);
    fprintf(out," yy  %15.7le  %15.7le  %15.7le \n", tM.qyy, eM.qyy, nM.qyy);
    fprintf(out," yz  %15.7le  %15.7le  %15.7le \n", tM.qyz, eM.qyz, nM.qyz);
    fprintf(out," zz  %15.7le  %15.7le  %15.7le \n", tM.qzz, eM.qzz, nM.qzz);
    fprintf(out,"\n");
    fprintf(out," Dipole moment     = %25.15le \n", tM.dz);
    fprintf(out," Quadrupole moment = %25.15le \n", q00);
    fclose (out);
}

void AnalyzeMOMoments (double *Cmat,int no2, int norb) {
    MomInts mint;

    FILE* in=open_file("RMOM.DAT");
    double *dx=new double[no2];
    double *dy=new double[no2];
    double *dz=new double[no2];
    double *wk=new double[norb*norb];
    for (int i = 0; i < no2; ++i) {
        dx[i]=dy[i]=dz[i]=0.0;
    }
    for (int i = 0; i < no2; ++i) {
        Fread(&mint,sizeof(MomInts),1,in);
        int ijr=mint.ijr;
        dx[ijr] = mint.dx;
        dy[ijr] = mint.dy;
        dz[ijr] = mint.dz;
    }
    fclose(in);
    SymmPack::sp_trans(norb,dx,Cmat,wk);
    SymmPack::sp_trans(norb,dy,Cmat,wk);
    SymmPack::sp_trans(norb,dz,Cmat,wk);
    FILE* out=create_file("mol_dipmom.dat");
    fprintf(out," MULTIPOLE MOMENT ANALYSIS \n");
    fprintf(out," units in bohr - hartree atomic units \n\n");
    fprintf(out," MO TRANSISITION DIPOLE MOMENTS \n\n");
    fprintf(out," orbital 1       orbital2      dx-dy-dz\n");
    int ij=0;
    for (int i=0; i<norb; ++i) {
        for (int j=0; j<=i; ++j,++ij) {
            fprintf(out," %12d %12d %15.6le %15.6le %15.6le\n",i,j,
                    dx[ij],dy[ij],dz[ij]);
        }
    }
    fclose ( out );
    delete [] wk;
    delete [] dz;
    delete [] dy;
    delete [] dx;
}

void
AnalyzeMOMoments (double *CmatA,double *CmatB,int no2, int norb) {
    MomInts mint;

    FILE* in=open_file("RMOM.DAT");
    double *dx=new double[no2];
    double *dy=new double[no2];
    double *dz=new double[no2];
    double *wk=new double[norb*norb];
    for (int i = 0; i < no2; ++i) {
        dx[i]=dy[i]=dz[i]=0.0;
    }
    for (int i = 0; i < no2; ++i) {
        Fread(&mint,sizeof(MomInts),1,in);
        int ijr=mint.ijr;
        dx[ijr] = mint.dx;
        dy[ijr] = mint.dy;
        dz[ijr] = mint.dz;
    }
    SymmPack::sp_trans(norb,dx,CmatA,wk);
    SymmPack::sp_trans(norb,dy,CmatA,wk);
    SymmPack::sp_trans(norb,dz,CmatA,wk);
    FILE* out=create_file("mol_dipmom.dat");
    fprintf(out," MULTIPOLE MOMENT ANALYSIS \n");
    fprintf(out," units in bohr - hartree atomic units \n\n");
    fprintf(out," ALPHA MO TRANSISITION DIPOLE MOMENTS \n\n");
    fprintf(out," orbital 1       orbital2      dx-dy-dz\n");
    int ij=0;
    for (int i=0; i<norb; ++i) {
        for (int j=0; j<=i; ++j,++ij) {
            fprintf(out," %12d %12d %15.6le %15.6le %15.6le\n",i,j,
                    dx[ij],dy[ij],dz[ij]);
        }
    }
    rewind(in);
    for (int i = 0; i < no2; ++i) {
        dx[i]=dy[i]=dz[i]=0.0;
    }
    for (int i = 0; i < no2; ++i) {
        Fread(&mint,sizeof(MomInts),1,in);
        int ijr=mint.ijr;
        dx[ijr] = mint.dx;
        dy[ijr] = mint.dy;
        dz[ijr] = mint.dz;
    }
    fclose(in);
    SymmPack::sp_trans(norb,dx,CmatB,wk);
    SymmPack::sp_trans(norb,dy,CmatB,wk);
    SymmPack::sp_trans(norb,dz,CmatB,wk);
    fprintf(out," BETA MO TRANSISITION DIPOLE MOMENTS \n\n");
    fprintf(out," orbital 1       orbital2      dx-dy-dz\n");
    ij=0;
    for (int i=0; i<norb; ++i) {
        for (int j=0; j<=i; ++j,++ij) {
            fprintf(out," %12d %12d %15.6le %15.6le %15.6le\n",i,j,
                    dx[ij],dy[ij],dz[ij]);
        }
    }
    fclose ( out );
    delete [] wk;
    delete [] dz;
    delete [] dy;
    delete [] dx;
}

} // end namespace

