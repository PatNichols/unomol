///////////////////////////////////////////////////////////////////////
//                   ATSCAT PROGRAM
///////////////////////////////////////////////////////////////////////
//   This program calculates the phase shifts,integrated cross sections,
// differential cross sections,and if wanted the physical solution for
// elastic positron atom scattering. The energies for these calculations
// are supplied by the user in the file scatin.dat. The potential is
// supplied by the user in a file vspolin.dat.
//   This program compiles under the GNU g++ compiler version 2.96 on
// a Linux box. With SGI you might try the -lang-std option.
// I really hope it won't compile with visual c++. Should you
// want to try to do so, you might try changing the names of the
// included header files to match those of your system.
///////////////////////////////////////////////////////////////////////
//  Added Born Closure to converge differential cross-sections
//  -------------------------------------------------------------------
//  Essentially One use the First Born Approximation
//  to calculate the scattering amplitude for _all_ values of L.
//  Then one removes the contribution from all partial waves
//  with L<=lmax. This leaves the contribution to the scattering
//  amplitude for all partial waves with L>lmax which is added
//  on to the scattering amplitude we have calculated with the
//  usual method. This should give us a _very_ accurate
//  approximation to the scattering amplitude. This could fail
//  for very high energies. So compare dxsec.XX.out with dxsec_bc.XX.dat!
///////////////////////////////////////////////////////////////////////
// copyright 2001 Patrick Jay Nichols and Texas Tech University
///////////////////////////////////////////////////////////////////////
// This program is released under the GNU Public License
///////////////////////////////////////////////////////////////////////
// send bug reports,comments,etc. to
//  ripjn@spudhammer.phys.ttu.edu
///////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <string>
#include <cmath>
#include <cfloat>
#include <iomanip>

// the value of a Rydberg in eV from NIST
#define RYDBERG 13.60569171
// small constant very approx equal to DBL_EPSILON
#define eps 1.e-15
#define SCAT_HUGE 1.e100 // default potential when r is about zero

class radial_mesh {
private:
    size_t npts;
    size_t nreg;
    size_t *nstp;
    double *rstt;
    double *delr;

    std::istream& read_from_stream(std::istream& in)
    {
        if (nstp) delete [] nstp;
        if (delr) delete [] delr;
        if (rstt) delete [] rstt;
        in >> nreg;
        rstt=new double[nreg+1];
        delr=new double[nreg];
        nstp=new size_t[nreg];
        double rs,re;
        in >> rs;
        npts = 0;
        for (int i=0; i<nreg; i++) {
            in >> re;
            in >> delr[i];
            if ( re < rs) {
                std::cerr<<"Your mesh has a hole in it!"<<std::endl;
                std::cerr<< " rst = " << rs << "\n";
                std::cerr<< "rend = " << re << "\n";
                exit(EXIT_FAILURE);
            }
            rstt[i] = rs;
            nstp[i] = (size_t) ( ( re - rs ) / delr[i] );
            rs = re;
            npts += nstp[i];
        }
        rstt[nreg]=rs;
        ++npts;
        return in;
    }
public:
    radial_mesh():npts(0),nreg(0),nstp(nullptr),rstt(nullptr),delr(nullptr) {}

    ~radial_mesh()
    {
        delete [] nstp;
        delete [] delr;
        delete [] rstt;
    }

    double rstart(size_t ir) const noexcept {
        return rstt[ir];
    }

    const double& deltar(size_t ir) const noexcept {
        return delr[ir];
    }

    const size_t& nsteps(size_t ir) const noexcept {
        return nstp[ir];
    }

    const size_t& num_regions() const noexcept {
        return nreg;
    }

    const size_t& num_points() const noexcept {
        return npts;
    }

    friend std::istream& operator >> ( std::istream& is, radial_mesh& mesh) {
        return mesh.read_from_stream(is);
    }

    void write_mesh(std::ostream& os)
    {
        os << "Radial Mesh\n";
        os << "# of regions = " << nreg << "\n";
        os << "# of points = " << npts << "\n";
        os << "rstart         ";
        os << "end         ";
        os << "step size   ";
        os << " # of steps\n";
        for (int i=0; i<nreg; ++i) {
            os << std::setw(15) << std::setprecision(8) << std::fixed << rstt[i] << " ";
            os << std::setw(15) << std::setprecision(8) << std::fixed << rstt[i+1] << " ";
            os << std::setw(15) << std::setprecision(8) << std::fixed << delr[i] << " ";
            os << std::setw(6)  << nstp[i] << "\n";
        }
    }
};


double si_int(double x) {
    if (x>=20.0) {
/////////////////////////////////////////////////////////////
// SEE ABRAMOWITZ AND STEGUN PG. 243 FOR THIS ONE
/////////////////////////////////////////////////////////////
        double f,g,cx,sx;
        f=g=1.0e0;
        double x2=x*x;
        double x2n=x2;
        double s=-1.0e0;
        double fact1=2.0e0;
        double fact2=6.0e0;
        for (size_t i=1; i<5; i++)  {
            f+=s*fact1/x2n;
            g+=s*fact2/x2n;
            s=-s;
            x2n*=x2;
            fact1*=(2*i+2)*(2*i+1);
            fact2*=(2*i+3)*(2*i+2);
        }
        sx = sin(x);
        cx = cos(x);
        fact1=-cx/x;
        fact2=-sx/x2;
        return (fact1*f+fact2*g);
    }
    if (x<=1.0) {
        if (x==0.0) return (DBL_MAX-M_PI_2);
////////////////////////////////////////////////////////
// USE INTEGRATION OF THE POWER SERIES OF THE INTEGRAND
////////////////////////////////////////////////////////
        double sum,x2,xfac,denom;
        sum=x;
        x2=x*x;
        xfac=x;
        denom=1.0;
        for (size_t i=3; i<=17; i+=2) {
            xfac*=(-x2);
            denom*=static_cast<double>(i*(i-1));
            sum+=xfac/denom/i;
        }
        return (sum-M_PI_2);
    }
/////////////////////////////////////////////////////
// USE 401 pt. 1/3 Simpson's Rule INTEGRATION
/////////////////////////////////////////////////////
    const int nst = 500;
    double dr = x/nst;
    double sum4 = 0.0;
    for (size_t i=0; i<nst; i+=2) {
        double r = (i+1)*dr;
        sum4+=sin(r)/r;
    }
    double sum2 = 0.0;
    for (size_t i=1; i<nst; i+=2) {
        double r = (i+1)*dr;
        sum2+=sin(r)/r;
    }
    double endpts=1.0+sin(x)/x;
    return ((dr/3.0)*(endpts+sum2*2.+sum4*4.)-M_PI_2);
}

////////////////////////////////////////////////////////////////////////
// SPLINE FIT CLASS
//
// written by Patrick J. Nichols,Texas Tech University 2001
//  send comments to patrick.j.nichols@ttu.edu
////////////////////////////////////////////////////////////////////////
//  this class borrows heavily from the Fortran77 programs
//   splint and splinb written by Joe Eccles 5/12/78 and
//   additions to these programs by Vasil Babamov (5/1/80)
////////////////////////////////////////////////////////////////////////
//
//  One needs to allocate two arrays.
//  One holds the funtion values for the points in question , the
//   other holds the domain values . These are passed to this class
//   along with the derivatives at the endpoints if you have them.
//  Note: the pointers are copied in a shallow manner so do not
//   free them until this class is destroyed.
////////////////////////////////////////////////////////////////////////
class spline_fit {
private:
    double *_c1,*_c2,*_c3;
    double *fk,*xk;
    size_t np;

/////////////////////////////////////////////////////////////////
// USE DIVIDED DIFFERENCES TO FIND THE DERIVATIVES AT THE
// END POINTS IF THE DERIVATIVES ARE NOT SUPPLIED BY THE
// USER
/////////////////////////////////////////////////////////////////
    void find_derivatives(double& d0,double& dn) {
        double f0,f1;
        f0=(fk[1]-fk[0])/(xk[1]-xk[0]);
        f1=(fk[2]-fk[1])/(xk[2]-xk[1]);
        d0=0.5*(f0-f1);
        f0=(fk[np-1]-fk[np-2])/(xk[np-1]-xk[np-2]);
        f1=(fk[np-2]-fk[np-3])/(xk[np-2]-xk[np-3]);
        dn=0.5*(f0-f1);
    }

/////////////////////////////////////////////////////////////////
// FINDS THE COEFFICIENTS OF THE SPLINE FIT
// BY FINDING THE DERIVATIVES AT THE KNOTS
/////////////////////////////////////////////////////////////////
    void find_coefficients(const double& d0,const double& dn) {
        double tst=1.0e4;
        size_t nm1=np-1;
        double dxl=xk[1]-xk[0];
        double dfl=(fk[1]-fk[1])*3.0/dxl;
        _c1[0]=d0;
        _c1[1]=d0;
        _c2[0]=0.0;
        _c2[1]=1.0;
        int ip1=2;
        for (int i=1; i<nm1; ip1++,i++) {
            double dxh=xk[ip1]-xk[i];
            double dfh=(fk[ip1]-fk[i])*3.0/dxh;
            _c1[ip1]=dfh-2.0*_c1[i]+(dfl-2.0*_c1[i]-_c1[i-1])*dxh/dxl;
            _c2[ip1]=-(2.0*_c2[i]+_c2[i-1])*dxh/dxl-2.0*_c2[i];
            dxl=dxh;
            if (fabs(_c2[ip1])<tst) {
                dfl=dfh;
                continue;
            }
            double fp=(dfl+dfh)*0.16666666666666666666667;
            double a=(fp-_c1[ip1])/_c2[ip1];
            double b=1.0/_c2[ip1];
            for (int j=0; j<=i; j++) {
                _c1[j]+=a*_c2[j];
                _c2[j]*=b;
            }
            _c1[ip1]=fp;
            _c2[ip1]=1.0;
            dfl=dfh;
        }
        double a=(dn-_c1[np-1])/_c2[np-1];
        _c1[0]+=a*_c2[0];
        for (int j=1; j<np; j++) {
            _c1[j]+=a*_c2[j];
            double cdiff=1.0/(xk[j]-xk[j-1]);
            double tp1=(fk[j]-fk[j-1])*cdiff;
            double tp2=(_c1[j]+_c1[j-1]);
            _c3[j-1]=(-2.0*tp1+tp2)*cdiff*cdiff;
            _c2[j-1]=(3.0*tp1-tp2-_c1[j-1])*cdiff;
        }
        _c1[np-1]=dn;
    }
public:
    spline_fit():_c1(0),_c2(0),_c3(0),fk(0),xk(0),np(0) {}

    spline_fit(double *f,double *x,size_t n):
        _c1(new double[n]),_c2(new double[n]),_c3(new double[n]),
        fk(0),xk(0),np(n) {
        fk = new double[n];
        xk = new double[n];
        memcpy(fk,f,sizeof(double)*n);
        memcpy(xk,x,sizeof(double)*n);
        double d0,d1;
        find_derivatives(d0,d1);
        find_coefficients(d0,d1);
    }

    spline_fit(double *f,double *x,const double& d0,
               const double& dn,size_t n):
        _c1(new double[n]),_c2(new double[n]),_c3(new double[n]),
        fk(0),xk(0),np(n) {
        fk = new double[n];
        xk = new double[n];
        memcpy(fk,f,sizeof(double)*n);
        memcpy(xk,x,sizeof(double)*n);
        find_coefficients(d0,dn);
    }

/// dtor not we do not delete xk or fk
    ~spline_fit() {
        delete [] xk;
        delete [] fk;
        delete [] _c3;
        delete [] _c2;
        delete [] _c1;
    }


/////////////////////////////////////////////////////////////////
// FINDS THE VALUE OF THE FUNCTION AT A POINT X BASED UPON THE
// SPLINE FIT
/////////////////////////////////////////////////////////////////
    double value(const double& x) const {
        double a;
        const size_t npm1=np-1;
        if (x==xk[npm1]) return fk[npm1];
        if (x==xk[0]) return fk[0];
        if (x>=xk[npm1]) {
            a=x-xk[npm1];
            return (((_c3[npm1]*a+_c2[npm1])*a+_c1[npm1])*a+fk[npm1]);
        }

        if (x<=xk[0]) {
            if (x==xk[0]) return fk[0];
            a=x-xk[0];
            return (((_c3[0]*a+_c2[0])*a+_c1[0])*a+fk[0]);
        }
// use binary search to find proper knot
        size_t nh=npm1;
        size_t nl=0;
        size_t i;
        do {
            i=(nl+nh)>>1;
            double diff=x-xk[i];
            if (diff==0.0) return fk[i];
            if (diff<0.0) nh=i;
            else nl=i;
        } while ((nh-nl)>1);
// interpolate value from this knot
        if (x==xk[nl]) return fk[nl];
        a=x-xk[nl];
        return (((_c3[nl]*a+_c2[nl])*a+_c1[nl])*a+fk[nl]);
    }
};


class sph_potential
{
public:
    sph_potential():sfit(nullptr),maxr(0),minr(0),alf(0),npts(0) {}

    ~sph_potential() {
        if (sfit) delete sfit;
    }

    double rmax() const noexcept { return maxr;}
    double rmin() const noexcept { return minr;}
    double alpha() const noexcept { return alf;} 
    int npoints() const noexcept {return npts;}

///////////////////////////////////////////////////
// Asymptotic potential for r>rmax
// spline fit for r<rmax
// you can add your own potentials here as well.
///////////////////////////////////////////////////
    double operator()(const double& r) const noexcept 
    {
        double v;
        if (r<=maxr) {
            if ( r>=minr) {
                return sfit->value(r);
            }        
            return SCAT_HUGE;
        }
        double r2=r*r;
        return (-alf/r2/r2);
    }

    friend std::istream& operator >> (std::istream& in,sph_potential& pot)
    {
        in >> pot.npts;
        in >> pot.maxr;
        in >> pot.alf;
        double *vp=new double[pot.npts];
        double *rp=new double[pot.npts];
        for (int i=0; i<pot.npts; i++) {
            in >> rp[i];
            in >> vp[i];
            vp[i]*=2.0;
        }
        pot.minr = rp[0];
        if (pot.sfit) delete pot.sfit;
        pot.sfit = new spline_fit(vp,rp,pot.npts);
        delete [] rp;
        delete [] vp;
        return in;
    }
private:
    spline_fit * sfit;
    double maxr;
    double minr;
    double alf;
    int npts;
};

constexpr double ipow(double x,int n) noexcept {
    if (!n) return 1.0;
    if (n<0) {
        n=-n;
        x=1.0/x;
    }
    double y=x;
    while (--n) y*=x;
    return y;
}


inline void sphbess(double& sphj,double& sphn,const double& x,int l) {
    const double fourpi = 12.566370614359172; // 4*pi
    int i, n;
    int nstop, nstrt, lp1, lmax;
    double rfact1, rfact2, rfact3, factr, ltest;
    double xi, xsq, hxsq, series, sj0, sj1, sn0, sn1, fact, sjl, sjlm1,
           sjlp1, sav, scale;

    if (x <= DBL_EPSILON) {
        sphj=0.0;
        if (!l) {
            sphj = 1.0;
        }
        sphn = -SCAT_HUGE;
        return;
    }
    if (x <= 0.1 && l > 1) {
        /**********************************************************
         * Ascending Bessel series
         **********************************************************/
        rfact1 = 1.0;
        nstop = 2 * l + 1;
        for (i = 1; i <= nstop; i += 2) {
            rfact1 *= ((double) i);
        }
        hxsq = 0.5 * x * x;
        rfact2 = 2.0 * l + 3.0;
        rfact3 = 2.0 * (2 * l + 5) * rfact2;
        series = 1.0 - hxsq / rfact2 + hxsq * hxsq / rfact3;
        sphj = series * ipow (x, l) / rfact1;
        /**********************************************************
         * ascending Neumann series
         **********************************************************/
        rfact1 = 1.0;
        lp1 = l + 1;
        nstop = 2 * l - 1;
        for (i = 1; i <= nstop; i += 2) {
            rfact1 *= (double) i;
        }
        rfact2 = 1.0 - 2.0 * l;
        rfact3 = 2.0 * (3.0 - 2 * l) * rfact1;
        series = 1.0 - hxsq / rfact2 - hxsq * hxsq / rfact3;
        sphn = -rfact1 * series / pow (x, lp1);
        return;
    }
    double csx = cos(x);
    double snx = sin(x);
    xi = 1.0 / x;
    sj0 = snx * xi;
    sj1 = (snx * xi - csx) * xi;
    sn0 = -csx * xi;
    sn1 = -(csx * xi + snx) * xi;
    if (l == 0) {
        sphj = sj0;
        sphn = sn0;
        return;
    } else {
        if (l == 1) {
            sphj = sj1;
            sphn = sn1;
            return;
        } else {
            /***************************************************************
             * Start the upward recurrence for the neumann function and
             *  for the Bessel functions of x*x|>(l*(l+1))
             ***************************************************************/
            for (i = 2; i <= l; i++) {
                fact = 2.0 * i - 1.0;
                sphn = fact * sn1 * xi - sn0;
                sn0 = sn1;
                sn1 = sphn;
            }
            xsq = x * x;
            ltest = (double) l *(l + 1);
            if (xsq >= ltest) {
                for (i = 2; i <= l; i++) {
                    fact = 2.0 * i - 1.0;
                    sphj = fact * sj1 * xi - sj0;
                    sj0 = sj1;
                    sj1 = sphj;
                }
                return;
            }
            /***************************************************************
             * start the downward recurrence for the Bessel function
             ***************************************************************/
            lmax = l + 20;
            sjlp1 = 0.0;
            sjl = 1.e-34;
            nstrt = lmax - 1;
            nstop = l + 1;
            for (n = nstrt; n >= nstop; n--) {
                factr = 2.0 * n + 1.0;
                sjlm1 = factr * sjl * xi - sjlp1;
                sjlp1 = sjl;
                sjl = sjlm1;
            }
            sav = sjl;
            nstrt = l;
            nstop = 1;
            for (n = nstrt; n >= nstop; n--) {
                factr = 2.0 * n + 1.0;
                sjlm1 = factr * sjl * xi - sjlp1;
                sjlp1 = sjl;
                sjl = sjlm1;
            }
            /****************************************************************
             * scale the Bessel Function
             ***************************************************************/
            scale = sj0 / sjl;
            sphj = scale * sav;
        }
        return;
    }
}

inline double si(double x) {
    if (x>=20.0) {
/////////////////////////////////////////////////////////////
// SEE ABRAMOWITZ AND STEGUN PG. 243 FOR THIS ONE
/////////////////////////////////////////////////////////////
        double f,g,cx,sx;
        f=g=1.0e0;
        double x2=x*x;
        double x2n=x2;
        double s=-1.0e0;
        double fact1=2.0e0;
        double fact2=6.0e0;
        for (size_t i=1; i<5; i++)  {
            f+=s*fact1/x2n;
            g+=s*fact2/x2n;
            s=-s;
            x2n*=x2;
            fact1*=(2*i+2)*(2*i+1);
            fact2*=(2*i+3)*(2*i+2);
        }
        sx = sin(x);
        cx = cos(x);
        fact1=-cx/x;
        fact2=-sx/x2;
        return (fact1*f+fact2*g);
    }
    if (x<=1.0) {
        if (x==0.0) return (DBL_MAX-M_PI_2);
////////////////////////////////////////////////////////
// USE INTEGRATION OF THE POWER SERIES OF THE INTEGRAND
////////////////////////////////////////////////////////
        double sum,x2,xfac,denom;
        sum=x;
        x2=x*x;
        xfac=x;
        denom=1.0;
        for (size_t i=3; i<=17; i+=2) {
            xfac*=(-x2);
            denom*=static_cast<double>(i*(i-1));
            sum+=xfac/denom/i;
        }
        return (sum-M_PI_2);
    }
/////////////////////////////////////////////////////
// USE 401 pt. 1/3 Simpson's Rule INTEGRATION
/////////////////////////////////////////////////////
    double qsum=0.0;
    double rstep=x/200.0;
    double endpts=1.0+sin(x)/x;
    double r=rstep*0.5;
    for (size_t i=2; i<=400; i+=2) {
        qsum+=sin(r)/r;
        r+=rstep;
    }
    double q4=qsum*4.0;
    qsum=0.0;
    r=rstep;
    for (size_t i=3; i<400; i+=2) {
        qsum+=sin(r)/r;
        r+=rstep;
    }
    double q2=2.0*qsum;
    return ((rstep/6.0)*(endpts+q2+q4)-M_PI_2);
}


///////////////////////////////////////////////////
// this propagates the solution out until
//  we reach convergence or we run out of points.
// We use trapeziodal integration here so keep
//  the step sizes small.
//
///////////////////////////////////////////////////
void do_calculation(double energy,const sph_potential& scat_pot, const radial_mesh& mesh,
                   double& shift,double& xsec,
                    const double& tolerance,int l)
{
    double rstop,atal,sphj,sphn;
    const double pi4 = 12.566370614359172; // 4*pi

    energy/=RYDBERG;
    double k=sqrt(energy);
    bool endit=false;
    size_t ncheck=0;
    double sum1 = 0.0;
    double sum2 = 0.0;
    double lfact=l*(l+1);
    double rs = mesh.rstart(0);
    sphbess(sphj,sphn,k*rs,l);
    double u = sphj * k;
    double v = scat_pot(rs);
    double dr = mesh.deltar(0);;
    double prod = u * v * rs * dr  * 0.5;
    double trap1 = prod * sphj * k;
    double trap2 = prod * sphn;
    size_t nreg = mesh.num_regions();
    for (int ir=0; ir<nreg; ir++) {
        rs=mesh.rstart(ir);
        dr=mesh.deltar(ir);
        size_t ns=mesh.nsteps(ir);
        prod = rs * u * v * dr * 0.5;
        trap1 = prod * sphj * k;
        trap2 = prod * sphn;
        for (int is=1; is<=ns; is++) {
            double r = rs + is * dr;
            sphbess(sphj,sphn,k*r,l);
            u=k*r*sphj*(1.0-sum2)+r*sphn*sum1;
            sum1 += trap1;
            sum2 += trap2;
            v=scat_pot(r);
            prod=r*u*v*dr*0.5;
            trap1=k*prod*sphj;
            trap2=sphn*prod;
            sum1 += trap1;
            sum2 += trap2;
/// test to see if we are in the asymptotic region //
            double tst=fabs(energy-lfact/r/r);
            if (tst<=eps || fabs(v)>tst || (fabs(tst-v)/tst)<0.999998 ) {
                continue;
            }
            ++ncheck;
            double afact=1.0-sum2;
            if (fabs(afact)<=eps) continue;
            double atalp1=atan(-sum1/afact/k);
            if (ncheck==1) atal=atalp1;
            if (ncheck<=19 || fabs(atalp1)<=eps ) continue;
            if (fabs((atalp1-atal)/atalp1)<tolerance) {
                endit=true;
                rstop = r;
                break;
            } 
            ncheck=0;
        }
        if (endit) break;
    }
    if (!endit) {
        std::cout<<"Covergence was not reached !!"<<std::endl;
    }
    fprintf(stdout,"L= %3u Energy = %15.10lf eV\n",l,energy*RYDBERG);
    fprintf(stdout," Rstop = %15.7le\n",rstop);
    fprintf(stdout,"sphj= %12.5le sphn=%12.5le uhomo= %12.5le v=%12.5le\n",
            sphj,sphn,u,v);
    fprintf(stdout,"t1= %12.5le t2= %12.5le sum1= %12.5le sum2= %12.5le\n",
            trap1,trap2,sum1,sum2);
    double factor=1.0-sum2;
    double tandl=-sum1/factor/k;
    shift=atan(tandl);
    double snsh = sin(shift);
    xsec=pi4*(2*l+1)*snsh*snsh/energy;
}

///////////////////////////////////////////////////
// this propagates the solution out until
//  we reach convergence or we run out of points.
// We use trapeziodal integration here so keep
//  the step sizes small.
// This version stores the physical solution
///////////////////////////////////////////////////
void do_calculation(double energy,const sph_potential& scat_pot,const radial_mesh& mesh,
                    double& shift,double& xsec,
                    double* uk,double* rk,const double& tolerance,int& unpts,int l) {
    double rstop,atal,sphj,sphn;
    const double pi4 = 12.566370614359172; // 4*pi
    energy/=RYDBERG;
    double k=sqrt(energy);
    bool endit=false;
    size_t kk = 0;
    size_t ncheck=0;
    double sum1 = 0.0;
    double sum2 = 0.0;
    double lfact=l*(l+1);
    double rs = mesh.rstart(0);
    sphbess(sphj,sphn,k*rs,l);
    double u = sphj * k;
    uk[kk] = u;
    rk[kk] = rs;
    ++kk;
    double v = scat_pot(rs);
    double dr = mesh.deltar(0);
    double prod = u * v * rs * dr  * 0.5;
    double trap1 = prod * sphj * k;
    double trap2 = prod * sphn;
    size_t nreg = mesh.num_regions();
    for (int ir=0; ir<nreg; ir++) {
        rs=mesh.rstart(ir);
        dr=mesh.deltar(ir);
        size_t ns=mesh.nsteps(ir);
        prod = rs * u * v * dr * 0.5;
        trap1 = prod * sphj * k;
        trap2 = prod * sphn;
        for (int is=1; is<=ns; is++) {
            double r = rs + is * dr;
            sphbess(sphj,sphn,k*r,l);
            u=k*r*sphj*(1.0-sum2)+r*sphn*sum1;
            uk[kk] = u;
            rk[kk] = r;
            ++kk;
            sum1 += trap1;
            sum2 += trap2;
            v=scat_pot(r);
            prod=r*u*v*dr*0.5;
            trap1=k*prod*sphj;
            trap2=sphn*prod;
            sum1 += trap1;
            sum2 += trap2;
/// test to see if we are in the asymptotic region //
            double tst=fabs(energy-lfact/r/r);
            if (tst<=eps || fabs(v)>tst || (fabs(tst-v)/tst)<0.999998 ) {
                continue;
            }
            ++ncheck;
            double afact=1.0-sum2;
            if (fabs(afact)<=eps) continue;
            double atalp1=atan(-sum1/afact/k);
            if (ncheck==1) atal=atalp1;
            if (ncheck<=19 || fabs(atalp1)<=eps ) continue;
            if (fabs((atalp1-atal)/atalp1)<tolerance) {
                endit=true;
                rstop=r;
                break;
            }
            ncheck=0;
        }
        if (endit) break;
    }
    if (!endit) {
        std::cout<<"Covergence was not reached !!"<<std::endl;
    }
    fprintf(stdout,"L= %3u Energy = %15.10lf eV\n",l,energy*RYDBERG);
    fprintf(stdout," Rstop = %15.7le\n",rstop);
    fprintf(stdout,"sphj= %12.5le sphn=%12.5le uhomo= %12.5le v=%12.5le\n",
            sphj,sphn,u,v);
    fprintf(stdout,"t1= %12.5le t2= %12.5le sum1= %12.5le sum2= %12.5le\n",
            trap1,trap2,sum1,sum2);
    double factor=1.0-sum2;
    double tandl=-sum1/factor/k;
    shift=atan(tandl);
    double snsh = sin(shift);
    xsec=pi4*(2*l+1)*snsh*snsh/energy;
    unpts = kk;
    double ufact = cos(shift)/factor;
    for (size_t m=0; m<kk; ++m) uk[m] *= ufact;
}


inline void write_u(double* u,double* r,int np,int ie,int l) {
    char nums[32];
    sprintf(nums,"uhomo.%03d.%03d.dat",ie,l);
    FILE * out = fopen(nums,"w");
    for (int i=0; i< np; i++) {
        fprintf(out,"%15.10lf %20.10le\n",r[i],u[i]);
    }
    fclose(out);
}

/////////////////////////////////////
// The first Born approximation to the
// total scattering amplitude.
///////////////////////////////////////////////////////////////////////
// let q=2*k*sin(theta/2)
// and rm=r where v(r) approximately equal alpha/r^4
// for us rm=rmax
//  fb(k,theta)= -1/q* integral(r=0:rmax) of sin(q*r)*r*V(r) dr +
//     alpha*0.5*(sin(q*rm)/rm/q/rm+cos(q*rm)/rm+q*si(rm))
// where si(x)= Si(x)- 0.5*pi =-integral(rm:infinity dt sin(t)/t)
// when theta->0
//   sin(q*r)/q*r->1
//  fb(k,0)= -integral(r=0:rmax) of r*r*V(r) dr + alpha/rm
///////////////////////////////////////////////////////////////////////
double bornf(const double& energy,const double& theta,
             const sph_potential& scat_pot,int nrad=500) {
    double rest,qint,end,wt,sum2,sum4,q,qrm,cqrm,sqrm;

    const double kappa=sqrt(energy/RYDBERG);
    const double alpha = scat_pot.alpha();
    const double rmax = scat_pot.rmax();
    double rstep=rmax/nrad;
    double trstep=2*rstep;
    wt=-rstep/3.0;
    double sum=0.0;
    double r=rstep;
    if (fabs(theta)<DBL_EPSILON) {
        for (size_t i=2; i<=nrad; i+=2) {
            sum+=r*r*scat_pot(r);
            r+=trstep;
        }
        sum4=4.0*sum;
        sum=0.0;
        r=trstep;
        for (size_t i=3; i<nrad; i+=2) {
            sum+=r*r*scat_pot(r);
            r+=trstep;
        }
        sum2=2.0*sum;
        end=-alpha/rmax/rmax;
        return (wt*(sum2+sum4+end)+alpha/rmax);
    }
    q=2.0*kappa*sin(theta*0.5);
    for (size_t i=2; i<=nrad; i+=2) {
        sum+=sin(q*r)*r*scat_pot(r);
        r+=trstep;
    }
    sum4=4.0*sum;
    sum=0.0;
    r=trstep;
    for (size_t i=3; i<nrad; i+=2) {
        sum+=sin(q*r)*r*scat_pot(r);
        r+=trstep;
    }
    sum2=2.0*sum;
    qrm=q*rmax;
    sqrm = sin(qrm);
    cqrm = cos(qrm);
    end=-sqrm*alpha/(rmax*rmax*rmax);
    qint=wt*(sum2+sum4+end)/q;
    rest=alpha*0.5*(sqrm/rmax/qrm+cqrm/rmax+q*si(qrm));
    return (qint+rest);
}

void psborn(const double& energy,const sph_potential& scat_pot,
            double* fbps,double* fbf,int lmax,int nrad=500) {
    double qsum2,qsum4,end;
    const double kappa=sqrt(energy/RYDBERG);
    const double alpha = scat_pot.alpha();
    const double rmax = scat_pot.rmax();
    double k2=kappa*kappa;
    double pika=M_PI*kappa*alpha;
    double krm=kappa*rmax;
    double tkrm=2.0*krm;
    double rstep=rmax/nrad;
    double trstep=2*rstep;
    double wt=rstep/3.0;
    double qsum=0.0;
    double r=rstep;
//////////////////////
// l=0 term
//////////////////////
    double skr=sin(krm);
    end=-skr*skr*alpha/rmax/rmax/rmax/rmax;
    for (size_t i=2; i<=nrad; i+=2) {
        skr=sin(kappa*r);
        qsum+=skr*skr*scat_pot(r);
        r+=trstep;
    }
    qsum4=4.0*qsum;
    qsum=0.0;
    r=trstep;
    for (size_t i=3; i<nrad; i+=2) {
        skr=sin(kappa*r);
        qsum+=skr*skr*scat_pot(r);
        r+=trstep;
    }
    qsum2=2.0*qsum;
    double qint=wt*(end+qsum2+qsum4)/k2;
    double qint2=1+(2.0*krm*krm-1.0)*cos(tkrm);
    qint2+=krm*sin(tkrm);
    qint2+=4.0*krm*krm*krm*si(tkrm);
    double fact=12.0*k2*rmax*rmax*rmax;
    qint2*=-alpha/fact;
    qint+=qint2;
    fbps[0]=atan(-kappa*qint);
    fbf[0]=-qint;
/////////////////////
// loop over l terms l!=0
/////////////////////
    double sj,sn;
    for (size_t l=1; l<=lmax; l++) {
        qsum=0.0;
        r=rstep;
        end=0.0;
        if (l==1) end=(alpha*k2)/9.0;
        for (size_t i=2; i<=nrad; i+=2) {
            sphbess(sj,sn,kappa*r,l);
            double r2=r*r;
            double r4=r2*r2;
            qsum+=r2*(scat_pot(r)+alpha/r4)*sj*sj;
            r+=trstep;
        }
        qsum4=4.0*qsum;
        qsum=0.0;
        r=trstep;
        for (size_t i=3; i<nrad; i+=2) {
            sphbess(sj,sn,kappa*r,l);
            double r2=r*r;
            double r4=r2*r2;
            qsum+=r2*(scat_pot(r)+alpha/r4)*sj*sj;
            r+=trstep;
        }
        qsum2=2.0*qsum;
        qint=wt*(end+qsum2+qsum4);
        fact=(2*l-1)*(2*l+1)*(2*l+3);
        qint2=-pika/fact;
        qint+=qint2;
        fbps[l]=atan(-kappa*qint);
        fbf[l]=(-(2.0*l+1.0)*qint);
    }
}

void born_proj(const double& energy,const sph_potential& scat_pot,double* fbl,
int lmax) {
    const size_t ntheta = 500;
    double pl[3];
    double tstep=M_PI/ntheta;
    double tstep2=2.0*tstep;
    double sum=0.0;
    double *lsum2=new double[lmax+1];
    double *lsum4=new double[lmax+1];
    for (size_t i=0; i<=lmax; i++) {
        lsum2[i]=0.0;
        lsum4[i]=0.0;
    }
    double t=tstep;
    for (size_t i=2; i<=ntheta; i+=2) {
        double fbk=bornf(energy,t,scat_pot);
        double cost = cos(t);
        double sint = sin(t);
        pl[0]=1.0;
        pl[1]=cost;
        lsum4[0]+=fbk*pl[0]*sint;
        lsum4[1]+=fbk*pl[1]*sint;
        for (size_t l=2; l<=lmax; l++) {
            pl[2]=(cost*(2*l-1)*pl[1]-(l-1)*pl[0])/l;
            lsum4[l]+=fbk*pl[2]*sint;
            pl[0]=pl[1];
            pl[1]=pl[2];
        }
        t+=tstep2;
    }
    t=tstep2;
    for (size_t i=3; i<ntheta; i+=2) {
        double fbk=bornf(energy,t,scat_pot);
        double cost = cos(t);
        double sint = sin(t);
        pl[0]=1.0;
        pl[1]=cost;
        lsum2[0]+=fbk*pl[0]*sint;
        lsum2[1]+=fbk*pl[1]*sint;
        for (size_t l=2; l<=lmax; l++) {
            pl[2]=(cost*(2*l-1)*pl[1]-(l-1)*pl[0])/l;
            lsum2[l]+=fbk*pl[2]*sint;
            pl[0]=pl[1];
            pl[1]=pl[2];
        }
        t+=tstep2;
    }
    double wt=tstep/6.0;
    for (size_t l=0; l<=lmax; l++) {
        fbl[l]=wt*(2.0*lsum2[l]+4.0*lsum4[l])*(2.0*l+1.0);
    }
    delete [] lsum2;
    delete [] lsum4;
}

void dxsec(double energy,const sph_potential& sfit,const double *pshift,
           std::ofstream& sout,int ie,int lmax,int nrad=500)
{
    char suffix_s[32];
    sprintf(suffix_s,".%0d.out",ie);
    std::string suffix(suffix_s);
    sprintf(suffix_s,".%0d.dat",ie);
    std::string suffix2(suffix_s);
////// open dxsec file
    std::string title="dxsec" + suffix2;
    std::ofstream out(title.c_str());
    title="dxsec_bc" + suffix2;
    std::ofstream bcout(title.c_str());
    title="bpshft" + suffix;
    std::ofstream pout(title.c_str());
//// here's the good stuff
    double k=sqrt(energy/RYDBERG);
    double kinv=1.0/k;
    double *r_amp=new double[lmax+1];
    double *i_amp=new double[lmax+1];
    double *b_amp_proj=new double[lmax+1];
    double *b_pshift=new double[lmax+1];
    double *b_amp=new double[lmax+1];
    born_proj(energy,sfit,b_amp_proj,lmax);
    psborn(energy,sfit,b_pshift,b_amp,lmax);
    for (size_t l=0; l<=lmax; l++) {
        double cs,ss;
        ss = sin(pshift[l]);
        cs = cos(pshift[l]);
        r_amp[l]=kinv*(2.0*l+1)*cs*ss;
        i_amp[l]=kinv*(2.0*l+1)*ss*ss;
        pout << std::setw(3) << l << " ";
        pout << std::setw(20) << std::setprecision(10) << std::scientific << b_pshift[l] << " ";
        pout << std::setw(20) << std::setprecision(10) << std::scientific << b_amp[l] << " ";
        pout << std::setw(20) << std::setprecision(10) << std::scientific << b_amp_proj[l] << "\n";
    }
    pout.close();
    double theta,cthet,dstep,dtheta,deg,pl[3],rsum,isum,b_sum_proj,b_sum;
    theta=0.0;
    dstep=5.0;
    size_t nstep=1+static_cast<size_t>(rint(180.0/dstep));
    dtheta=dstep*M_PI/180.0;
    deg=0.0;
    sout << "Angle(degrees)     DCS no BornClosure       DCS w/ Born Closure ";
    sout << "     DCS                   ";
    sout << "     DCS first Born Approx.\n";
    sout << "                                            (calc amps)         ";
    sout << "     (proj amps)   )      \n";
    for (size_t i=0; i<nstep; i++) {
        double cthet=cos(theta);
        double fbk=bornf(energy,theta,sfit);
        pl[0]=1.0;
        pl[1]=cthet;
        rsum=r_amp[0]+r_amp[1]*pl[1];
        isum=i_amp[0]+i_amp[1]*pl[1];
        b_sum_proj=b_amp_proj[0]+b_amp_proj[1]*pl[1];
        b_sum=b_amp[0]+b_amp[1]*pl[1];
        for (size_t l=2; l<=lmax; l++) {
            pl[2]=(cthet*(2.0*l-1)*pl[1]-(l-1.0)*pl[0])/l;
            rsum+=r_amp[l]*pl[2];
            isum+=i_amp[l]*pl[2];
            b_sum_proj+=b_amp_proj[l]*pl[2];
            b_sum+=b_amp[l]*pl[2];
            pl[0]=pl[1];
            pl[1]=pl[2];
        }
        double dxsec0=(rsum*rsum+isum*isum);
// find R-closure correction
        double deltaf=fbk-b_sum_proj;
        double deltaf2=fbk-b_sum;
        double dxsec1=dxsec0+2.0*rsum*deltaf+deltaf*deltaf;
        double dxsec2=fbk*fbk;
        double dxsec3=dxsec0+2.0*rsum*deltaf2+deltaf2*deltaf2;
        sout << std::setw(15) << std::setprecision(10) << std::fixed << deg << " ";
        sout << std::setw(24) << std::setprecision(15) << std::scientific << dxsec0 << " ";
        sout << std::setw(24) << std::setprecision(15) << std::scientific << dxsec3 << " ";
        sout << std::setw(24) << std::setprecision(15) << std::scientific << dxsec1 << " ";
        sout << std::setw(24) << std::setprecision(15) << std::scientific << dxsec2 << "\n";
        out << std::setw(15) << std::setprecision(10) << std::fixed << deg << " ";
        out << std::setw(24) << std::setprecision(16) << std::scientific << dxsec0 << "\n";
        bcout << std::setw(15) << std::setprecision(10) << std::fixed << deg << " ";
        bcout << std::setw(24) << std::setprecision(16) << std::scientific << dxsec1 << "\n";
        deg+=dstep;
        theta+=dtheta;
    }
    sout<<"\n";
    out.close();
    bcout.close();
////////////////////////////////////////////////////////
// compute integrated cross section from
// integrating DCS
////////////////////////////////////////////////////////
    double sthet;
    double sum0=0.0;
    double sum1=0.0;
    double sum2=0.0;
    double sum3=0.0;
    double dt=M_PI/180;
    double t=dt;
    double wt=2.0*M_PI*dt;
    for (size_t i=1; i<180; i++) {
        double fbk=bornf(energy,t,sfit);
        sthet = sin(t);
        cthet = cos(t);
        pl[0]=1.0;
        pl[1]=cthet;
        rsum=r_amp[0]+r_amp[1]*pl[1];
        isum=i_amp[0]+i_amp[1]*pl[1];
        b_sum_proj=b_amp_proj[0]+b_amp_proj[1]*pl[1];
        b_sum=b_amp[0]+b_amp[1]*pl[1];
        for (size_t l=2; l<=lmax; l++) {
            pl[2]=(cthet*(2.0*l-1)*pl[1]-(l-1.0)*pl[0])/l;
            rsum+=r_amp[l]*pl[2];
            isum+=i_amp[l]*pl[2];
            b_sum_proj+=b_amp_proj[l]*pl[2];
            b_sum+=b_amp[l]*pl[2];
            pl[0]=pl[1];
            pl[1]=pl[2];
        }
        double dxsec0=(rsum*rsum+isum*isum);
// find R-closure correction
        double deltaf=fbk-b_sum_proj;
        double deltaf2=fbk-b_sum;
        double dxsec1=dxsec0+2.0*rsum*deltaf+deltaf*deltaf;
        double dxsec2=fbk*fbk;
        double dxsec3=dxsec0+2.0*rsum*deltaf2+deltaf2*deltaf2;
        sum0+=dxsec0*sthet;
        sum1+=dxsec1*sthet;
        sum2+=dxsec2*sthet;
        sum3+=dxsec3*sthet;
        t+=dt;
    }
    double aixsec0=wt*sum0;
    double aixsec1=wt*sum1;
    double aixsec2=wt*sum2;
    double aixsec3=wt*sum3;
    sout<<"Integrated Cross Section from Phase Shifts = "<<std::endl;
    sout<<"           ";
    sout << std::setw(24) << std::setprecision(15) << std::scientific << aixsec0 << "\n";

    sout<<"Integrated Cross Section with Born Closure (Calc Proj)= "<<std::endl;
    sout<<"           ";
    sout << std::setw(24) << std::setprecision(15) << std::scientific << aixsec3 << "\n";

    sout<<"Integrated Cross Section with Born Closure (Proj Amps)= "<<std::endl;
    sout<<"           ";
    sout << std::setw(24) << std::setprecision(15) << std::scientific << aixsec1 << "\n";

    sout<<"Integrated Cross Section from First Born DCS amps = "<<std::endl;
    sout<<"           ";
    sout << std::setw(24) << std::setprecision(15) << std::scientific << aixsec2 << "\n";

    isum=0.0;
    pl[0]=1.0;
    pl[1]=1.0;
    isum=i_amp[0]+i_amp[1]*pl[1];
    for (size_t l=2; l<=lmax; l++) {
        pl[2]=((2.0*l-1)*pl[1]-(l-1.0)*pl[0])/l;
        isum+=i_amp[l]*pl[2];
        pl[0]=pl[1];
        pl[1]=pl[2];
    }
    isum*=4.0*M_PI*kinv;

    sout<<"Integated Cross Section from Optical Theorem="<<std::endl;
    sout<<"           ";
    sout << std::setw(24) << std::setprecision(15) << std::scientific << isum << "\n";
//////// free allocated arrays //////////////////
    delete [] b_amp;
    delete [] b_pshift;
    delete [] b_amp_proj;
    delete [] i_amp;
    delete [] r_amp;
}

int main() {
    int npts,prtopt;
    int nnrg;
    double *u,*r;
    int lmax = 10;
// open the file scatin.dat and read the parameters for this run
    std::ifstream in("scatin.dat");
    if (!in) {
        std::cerr<<"Could not open the file scatin.dat! does it exist ?"<<std::endl;
        exit(EXIT_FAILURE);
    }
    in >> lmax;
    in >> prtopt;
    if (lmax>50) {
        std::cerr<<"max l value is > 50 "<<std::endl;
        std::cerr<<"This will cause catastrophic failure in "<<std::endl
                 <<"the spherical Bessel function."<<std::endl;
        std::cerr<<"You do not want that"<<std::endl;
        in.close();
        exit(EXIT_FAILURE);
    }
    double tolerance;
    in >> tolerance;
    radial_mesh rmesh;
    in >> rmesh;
    in >> nnrg;
    double *enrg=new double[nnrg];
    for (int i=0; i<nnrg; i++) {
        in >> enrg[i];
    }
    in.close();
// read potential in and spline fit it
    sph_potential spot;
    std::ifstream vin("vsp.dat");
    if (!vin) {
        std::cerr<<"Could not open the file vsp.dat! does it exist ?"<<std::endl;
        exit(EXIT_FAILURE);
    }
    vin >> spot;
    vin.close();
////////////////////////////////////////
// write input parameters to file
    std::ofstream out("scat.out");
    out << " Elastic Atomic Scattering Program for Patmol\n\n";
    out<<" The max L value is "<<lmax  << "\n";
    out<<" The asymptotic dipole polarizablity is: ";
    out<< std::setw(15) << std::setprecision(10) << std::fixed << spot.alpha() << "\n";
    out<<" The asymptotic radial value is : ";
    out<< std::setw(15) << std::setprecision(10) << std::fixed << spot.rmax() << "\n";
    out<<" The tolerance is : ";
    out<< std::setw(15) << std::setprecision(10) << std::scientific << tolerance << "\n";
    out<< "\n";
    rmesh.write_mesh(out);
    out << "\n";
    out<<" There are "<<nnrg<<" energies \n";
    out<<std::endl;
    out<<" The energies are :"<<std::endl;
    for (int i=0; i<nnrg; i++) {
        out<< std::setw(20) << std::setprecision(10) << std::fixed << enrg[i] << "\n";
    }
    out<<"\n";
///////////////////////////////////
// now we do the calculations
    double *pshift=new double[lmax+1];
    double *xsec=new double[lmax+1];
    std::ofstream pout("pshift.out");
    std::ofstream xout("xsec.out");
    if (prtopt) {
        size_t tnpts=rmesh.num_points();
        u=new double[tnpts];
        r=new double[tnpts];
    }    
    for (int ie=0; ie<nnrg; ie++) {
        double e=enrg[ie];
        double txsec=0.0;
        double tpshift=0.0;
        if (!prtopt) {
            for (int l=0; l<=lmax; l++) {
                double psh,xs;
                do_calculation(e,spot,rmesh,psh,xs,tolerance,l);
                pshift[l] = psh;
                tpshift += psh;
                xsec[l] = xs;
                txsec += xs;
            }
        } else {
            for (int l=0; l<=lmax; l++) {
                double xs,psh;
                int np;
                do_calculation(e,spot,rmesh,psh,xs,u,r,tolerance,np,l);
                write_u(u,r,np,ie,l);
                pshift[l] = psh;
                tpshift += psh;
                xsec[l] = xs;
                txsec += xs;
            }
        }
// write output to scat.out
        out<<"Energy = ";
        out<< std::setw(15) << std::setprecision(10) << std::fixed << (e/RYDBERG) << " Rydbergs\n";
        out<< std::setw(15) << std::setprecision(10) << std::fixed << (e) << " eV\n";
        out << "L  Phase Shift          Cross-Section\n";
        out << "--------------------------------------------\n";
        for (int l=0; l<=lmax; l++) {
            out << std::setw(2) << l << " ";
            out<< std::setw(20) << std::setprecision(10) << std::scientific << (pshift[l]) << " ";
            out<< std::setw(20) << std::setprecision(10) << std::scientific << (xsec[l]) << "\n";
        }
        out << "   ";
        out<< std::setw(20) << std::setprecision(10) << std::scientific << (tpshift) << " ";
        out<< std::setw(20) << std::setprecision(10) << std::scientific << (txsec) << "\n\n";
// write output to pshift.out
        for (int l=0; l<=lmax; l++) {
            pout << std::setw(3) << l << " ";
            pout<< std::setw(20) << std::setprecision(10) << std::scientific << (e) << " ";
            pout<< std::setw(20) << std::setprecision(10) << std::scientific << (pshift[l]) << "\n";
        }
// write output to xsec.out
        xout<< std::setw(20) << std::setprecision(10) << std::fixed << (e) << " ";
        xout<< std::setw(20) << std::setprecision(10) << std::scientific << (txsec) << "\n\n";
// write out diff xsections
        dxsec(e,spot,pshift,out,ie,lmax);
        out<<"Integrated Cross Section from Integration = "<<std::endl;
        out<<"           ";
        out << std::setw(24) << std::setprecision(15) << std::scientific << txsec << "\n\n";
    }
    if (prtopt) {
        delete [] r;
        delete [] u;
    }
    delete [] xsec;
    delete [] pshift;
    delete [] enrg;
    out.close();
    pout.close();
    xout.close();
}

