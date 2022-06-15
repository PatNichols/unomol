#include "TwoElectronInts.hpp"

namespace unomol {

#define UNO_MASK 0xF
#define UNO_SHIFT 4U
#define UNO_SHIFT2 8U


//#define UNOMOL_MD_INTS


struct PrimBlock {
    double p[3];
    double pa[3];
    double pexp;
    double sab;
    double c12;

    PrimBlock() = delete;

    PrimBlock(
        const double& a1,const double& a2,
        const double& c1,const double& c2,
        const double * r1,const double * r2,const double& rsq)
    {
        pexp = a1 + a2;
        double abi = 1./pexp;
        c12 = c1 * c2;
        sab = std::exp(-a1*a2*rsq*abi);
        for (int k=0;k<3;++k) p[k] = (a1*r1[k]+a2*r2[k])*abi;
        for (int k=0;k<3;++k) pa[k] = p[k] - r1[k];
    }
    
};

struct ShellBlock {
    std::vector< PrimBlock > prims;
    double r12[3];
    double r12sq;
    int l1,l2,lt;
    int off1,off2;

    ShellBlock() = delete;
    
    explicit ShellBlock(const Shell& sh1, const Shell& sh2,
        const Center * centers,
        int off1_,int off2_)
    {
        const double * r1 = (centers+sh1.center())->r_vec();
        const double * r2 = (centers+sh2.center())->r_vec();
        off1 = off1_;
        off2 = off2_;
        l1 = sh1.Lvalue();
        l2 = sh2.Lvalue();
        lt = l1 + l2;
        for (int k=0;k<3;++k) r12[k] = r1[k] - r2[k];
        r12sq = 0.0;
        for (int k=0;k<3;++k) r12sq += r12[k] * r12[k];
        int np1 = sh1.number_of_prims();
        int np2 = sh2.number_of_prims();
        prims.clear();
        prims.reserve(np1*np2);
        for (int ip1=0;ip1<np1;++ip1) {
            for (int ip2=0;ip2<np2;++ip2) {
                prims.push_back( 
                    PrimBlock(sh1.alf(ip1),sh2.alf(ip2),
                        sh1.cof(ip1),sh2.cof(ip2),
                        r1,r2,r12sq) );
            }
        }
    }
};


void TwoElectronInts::calculate(const Basis& basis) {
    const double SRterm= 34.9868366552497250;
    const double threshold=1.e-14;
    const double ints_eps = 1.e-14;
    const Shell* shells(basis.shell_ptr());
    const Center* centers(basis.center_ptr());
    const AuxFunctions& aux(*basis.auxfun_ptr());
    const int nshell=basis.number_of_shells();
    const int maxl=basis.maxLvalue();
    const int maxlst=aux.maxLstates();
    const int maxl4 = maxlst * maxlst * maxlst * maxlst;
    Rys rys(maxl);
    std::vector< ShellBlock > shell_blocks;
    std::vector< double > norms(maxl4);
    std::vector< std::uint32_t > lstates(maxl4);
    std::vector< TwoInts > sints(maxl4);
    cache.open_for_writing();
    putils::Stopwatch timer;
    timer.start();
    /// Create ShellBlocks
    shell_blocks.clear();
    for (int ish=0; ish<nshell; ++ish) {
        for (int jsh=0; jsh<=ish; ++jsh) {
            int off1 = basis.offset(ish);
            int off2 = basis.offset(jsh);
            const Shell& sh1 = shells[ish];
            const Shell& sh2 = shells[jsh];
            shell_blocks.push_back(
                ShellBlock(sh1,sh2,centers,
                    off1,off2) );
        }
    }
    // compute all erints over pairs of ShellBlocks
    int nblks = shell_blocks.size();
    std::cerr << " # of blocks = " << nblks << "\n";
    int st_block = start * ( start + 1 ) / 2;
    int ncalc = 0;
    std::cerr << " start block = " << st_block << "\n";
    int end_blks = (start+1) * ( start + 2) / 2;
    for (int ib1=st_block; ib1 < nblks; ++ib1)  {
        for (int ib2=0; ib2< nblks; ++ib2) {
            const ShellBlock& sb1 = shell_blocks[ib1];
            const ShellBlock& sb2 = shell_blocks[ib2];
            int ir0 = sb1.off1;
            int jr0 = sb1.off2;
            int kr0 = sb2.off1;
            int lr0 = sb2.off2;
            if ( jr0 > ir0 ) break;
            if ( kr0 > ir0 ) break;
            int nls1 = aux.number_of_lstates(sb1.l1);
            int nls2 = aux.number_of_lstates(sb1.l2);
            int nls3 = aux.number_of_lstates(sb2.l1);
            int nls4 = aux.number_of_lstates(sb2.l2);
                    int knt=0;
                    for (int ils=0; ils<nls1;++ils) {
                        int ir = ir0 + ils;
                        for (int jls=0; jls<nls2;++jls) {
                            int jr = jr0 + jls;
                            if ( jr > ir ) break;
                            for (int kls=0; kls<nls3;++kls) {
                                int kr = kr0 + kls;
                                if ( kr > ir) break;
                                for (int lls=0; lls<nls4;++lls) {
                                    int lr = lr0 + lls;
                                    if ( lr > kr || ( ir == kr && lr > jr) ) break;
                                    sints[knt].val=0.0;
                                    sints[knt].i=(unsigned int)ir;
                                    sints[knt].j=(unsigned int)jr;
                                    sints[knt].k=(unsigned int)kr;
                                    sints[knt].l=(unsigned int)lr;
                                    unsigned int l12 = (ils<<UNO_SHIFT) + jls;
                                    unsigned int l34=(kls<<UNO_SHIFT)+lls;
                                    lstates[knt]=(l12<<UNO_SHIFT2)+l34;
                                    norms[knt] =
                                        aux.normalization_factor(sb1.l1,ils)*
                                        aux.normalization_factor(sb1.l2,jls)*
                                        aux.normalization_factor(sb2.l1,kls)*
                                        aux.normalization_factor(sb2.l2,lls);
                                    ++knt;
                                }
                            }
                        }
                    }
                    if (!knt) continue;
                    ncalc += knt;
            int nroots = (sb1.lt + sb2.lt)/2 + 1;
            // end loop to find int l state data
            int np1 = sb1.prims.size();
            int np2 = sb2.prims.size();
            for (int ip1=0; ip1<np1; ++ip1) {
                for (int ip2=0; ip2<np2; ++ip2) {
                    const PrimBlock& sp1 = sb1.prims[ip1];
                    const PrimBlock& sp2 = sb2.prims[ip2];
                    double txp = sp1.pexp + sp2.pexp;
                    double sr = SRterm * sp1.sab * sp2.sab /
                                (sp1.pexp * sp2.pexp * sqrt(txp));
                    if (sr < 1.e-12) continue;
                    sr *= sp1.c12 * sp2.c12;        
                    rys.Recur(sp1.p ,sp2.p,
                              sp1.pa,sp2.pa,
                              sp1.pexp,
                              sp2.pexp,
                              txp,
                              sb1.lt,
                              sb2.lt,
                              nroots);
                    for (int kc=0; kc<knt; ++kc) {
                        unsigned int key=lstates[kc];
                        int lls=key&UNO_MASK;
                        key>>=UNO_SHIFT;
                        int kls=key&UNO_MASK;
                        key>>=UNO_SHIFT;
                        int jls=key&UNO_MASK;
                        key>>=UNO_SHIFT;
                        int ils=key&UNO_MASK;
                        const int *lv1 = aux.l_vector(sb1.l1,ils);
                        const int *lv2 = aux.l_vector(sb1.l2,jls);
                        const int *lv3 = aux.l_vector(sb2.l1,kls);
                        const int *lv4 = aux.l_vector(sb2.l2,lls);
                        double sum=
                            rys.Shift(sb1.r12,sb2.r12,lv1,lv2,lv3,lv4,nroots);
                        sints[kc].val += sum * sr * norms[kc];
                    } // end loop over l states combos
                }  // end loop over prim shell block 2
            } // end loop over prim shell block 1
            for (int kc=0; kc<knt; ++kc) {
                if ( fabs(sints[kc].val) >= 1.e-14) {
                    cache.write(&sints[kc],1);
                }
            }
        }  // end loop over shell block 2
    } // end loop over shell block 1
    timer.stop();
    cache.close();
    std::cerr << "Time for Two Electrons Integrals = " << timer.elapsed_time() << " seconds\n";
    size_t nb = cache.total_size()/sizeof(TwoInts);
    std::cerr << " # of integrals written  = " << nb << "\n";
    std::cerr << " # if integrals calc.    = " << ncalc << "\n";
}


void
TwoElectronInts::formGmatrix(const double* Pmat,double *Gmat) {
    constexpr const int BINSIZE=8196;
    TwoInts sints[BINSIZE];

    cache.open_for_reading();
    size_t nints = cache.total_size()/sizeof(TwoInts);
    size_t nbin = nints/BINSIZE;
    size_t nextra = nints%BINSIZE;
    for (int ibin=0; ibin<nbin; ++ibin) {
        cache.read(sints,BINSIZE);
        for (int ix=0; ix<BINSIZE; ++ix) {
            double val=(sints+ix)->val;
            int i=(sints+ix)->i;
            int j=(sints+ix)->j;
            int k=(sints+ix)->k;
            int l=(sints+ix)->l;
            int ii=i*(i+1)/2;
            int ij=ii+j;
            int ik=ii+k;
            int il=ii+l;
            int jk,jl;
            int kk = ( k * ( k + 1 ) ) / 2;
            int kl = kk + l;
            if (j >= k) {
                int jj = (j * ( j + 1 ) ) / 2;
                jk = jj + k;
                jl = jj + l;
            } else {
                jk = kk + j;
                if (j>l) {
                    jl=j*(j+1)/2+l;
                } else {
                    jl=l*(l+1)/2+j;
                }
            }
            double da=val*2.0*Pmat[ij];
            double db=val*2.0*Pmat[kl];
            double sjl=val*Pmat[ik];
            double sjk=val*Pmat[il];
            double sik=val*Pmat[jl];
            double sil=val*Pmat[jk];
            if (k!=l) {
                db=db+db;
                Gmat[ik]-=sik;
                if (i!=j && j>=k) Gmat[jk]-=sjk;
            }
            Gmat[il]-=sil;
            Gmat[ij]+=db;
            if (i!=j && j>=l) Gmat[jl]-=sjl;
            if (ij!=kl) {
                if (i!=j) da=da+da;
                if (j<=k) {
                    Gmat[jk]-=sjk;
                    if (i!=j && i<=k) Gmat[ik]-=sik;
                    if (k!=l && j<=l) Gmat[jl]-=sjl;
                }
                Gmat[kl]+=da;
            }
        }
    }
    if (nextra) {
        cache.read(sints,nextra);
        for (int ix=0; ix<nextra; ++ix) {
            double val=(sints+ix)->val;
            int i=(sints+ix)->i;
            int j=(sints+ix)->j;
            int k=(sints+ix)->k;
            int l=(sints+ix)->l;
            int ii=i*(i+1)/2;
            int ij=ii+j;
            int ik=ii+k;
            int il=ii+l;
            int jk,jl;
            int kk = ( k * ( k + 1 ) ) / 2;
            int kl = kk + l;
            if (j >= k) {
                int jj = (j * ( j + 1 ) ) / 2;
                jk = jj + k;
                jl = jj + l;
            } else {
                jk = kk + j;
                if (j>l) {
                    jl=j*(j+1)/2+l;
                } else {
                    jl=l*(l+1)/2+j;
                }
            }
            double da=val*2.0*Pmat[ij];
            double db=val*2.0*Pmat[kl];
            double sjl=val*Pmat[ik];
            double sjk=val*Pmat[il];
            double sik=val*Pmat[jl];
            double sil=val*Pmat[jk];
            if (k!=l) {
                db=db+db;
                Gmat[ik]-=sik;
                if (i!=j && j>=k) Gmat[jk]-=sjk;
            }
            Gmat[il]-=sil;
            Gmat[ij]+=db;
            if (i!=j && j>=l) Gmat[jl]-=sjl;
            if (ij!=kl) {
                if (i!=j) da=da+da;
                if (j<=k) {
                    Gmat[jk]-=sjk;
                    if (i!=j && i<=k) Gmat[ik]-=sik;
                    if (k!=l && j<=l) Gmat[jl]-=sjl;
                }
                Gmat[kl]+=da;
            }
        }
    }
    cache.close();
}

void
TwoElectronInts::formGmatrix(const double* PmatA,const double *PmatB,
                             double *GmatA,double *GmatB) {
    constexpr const int BINSIZE=8196;
    TwoInts sints[BINSIZE];

    cache.open_for_reading();
    size_t nints = cache.total_size()/sizeof(TwoInts);
    size_t nbin = nints/BINSIZE;
    size_t nextra = nints%BINSIZE;
    for (size_t ibin=0; ibin<nbin; ++ibin) {
        cache.read(sints,BINSIZE);
        for (int ix=0; ix<BINSIZE; ++ix) {
            double val=(sints+ix)->val;
            int i=(sints+ix)->i;
            int j=(sints+ix)->j;
            int k=(sints+ix)->k;
            int l=(sints+ix)->l;
            int ii=i*(i+1)/2;
            int ij=ii+j;
            int ik=ii+k;
            int il=ii+l;
            int jk,jl;
            int kl=k*(k+1)/2+l;
            if (j>k) {
                jk=j*(j+1)/2+k;
            } else {
                jk=k*(k+1)/2+j;
            }
            if (j>l) {
                jl=j*(j+1)/2+l;
            } else {
                jl=l*(l+1)/2+j;
            }
            double da=val*(PmatA[ij]+PmatB[ij]);
            double db=val*(PmatA[kl]+PmatB[kl]);
            double sjlA=val*PmatA[ik];
            double sjkA=val*PmatA[il];
            double sikA=val*PmatA[jl];
            double silA=val*PmatA[jk];
            double sjlB=val*PmatB[ik];
            double sjkB=val*PmatB[il];
            double sikB=val*PmatB[jl];
            double silB=val*PmatB[jk];
            if (k!=l) {
                db=db+db;
                GmatA[ik]-=sikA;
                GmatB[ik]-=sikB;
                if (i!=j && j>=k) {
                    GmatA[jk]-=sjkA;
                    GmatB[jk]-=sjkB;
                }
            }
            GmatA[il]-=silA;
            GmatA[ij]+=db;
            GmatB[il]-=silB;
            GmatB[ij]+=db;
            if (i!=j && j>=l) {
                GmatA[jl]-=sjlA;
                GmatB[jl]-=sjlB;
            }
            if (ij!=kl) {
                if (i!=j) da=da+da;
                if (j<=k) {
                    GmatA[jk]-=sjkA;
                    GmatB[jk]-=sjkB;
                    if (i!=j && i<=k) {
                        GmatA[ik]-=sikA;
                        GmatB[ik]-=sikB;
                    }
                    if (k!=l && j<=l) {
                        GmatA[jl]-=sjlA;
                        GmatB[jl]-=sjlB;
                    }
                }
                GmatA[kl]+=da;
                GmatB[kl]+=da;
            }
        }
    }
    if (nextra!=0) {
        cache.read(sints,nextra);
        for (int ix=0; ix<nextra; ++ix) {
            double val=(sints+ix)->val;
            int i=(sints+ix)->i;
            int j=(sints+ix)->j;
            int k=(sints+ix)->k;
            int l=(sints+ix)->l;
            int ii=i*(i+1)/2;
            int ij=ii+j;
            int ik=ii+k;
            int il=ii+l;
            int jk,jl;
            int kl=k*(k+1)/2+l;
            if (j>k) {
                jk=j*(j+1)/2+k;
            } else {
                jk=k*(k+1)/2+j;
            }
            if (j>l) {
                jl=j*(j+1)/2+l;
            } else {
                jl=l*(l+1)/2+j;
            }
            double da=val*(PmatA[ij]+PmatB[ij]);
            double db=val*(PmatA[kl]+PmatB[kl]);
            double sjlA=val*PmatA[ik];
            double sjkA=val*PmatA[il];
            double sikA=val*PmatA[jl];
            double silA=val*PmatA[jk];
            double sjlB=val*PmatB[ik];
            double sjkB=val*PmatB[il];
            double sikB=val*PmatB[jl];
            double silB=val*PmatB[jk];
            if (k!=l) {
                db=db+db;
                GmatA[ik]-=sikA;
                GmatB[ik]-=sikB;
                if (i!=j && j>=k) {
                    GmatA[jk]-=sjkA;
                    GmatB[jk]-=sjkB;
                }
            }
            GmatA[il]-=silA;
            GmatA[ij]+=db;
            GmatB[il]-=silB;
            GmatB[ij]+=db;
            if (i!=j && j>=l) {
                GmatA[jl]-=sjlA;
                GmatB[jl]-=sjlB;
            }
            if (ij!=kl) {
                if (i!=j) da=da+da;
                if (j<=k) {
                    GmatA[jk]-=sjkA;
                    GmatB[jk]-=sjkB;
                    if (i!=j && i<=k) {
                        GmatA[ik]-=sikA;
                        GmatB[ik]-=sikB;
                    }
                    if (k!=l && j<=l) {
                        GmatA[jl]-=sjlA;
                        GmatB[jl]-=sjlB;
                    }
                }
                GmatA[kl]+=da;
                GmatB[kl]+=da;
            }
        }
    }
    cache.close();
}

}
