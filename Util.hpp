#ifndef UNOMOL_UTIL_HPP
#define UNOMOL_UTIL_HPP
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iomanip>
using namespace std;

namespace unomol {

void fatal_error(const char* message);

template <class T> inline T** new_matrix(const int nrows,const int ncols) {
    T** m=new T*[nrows];
    for (int i=0; i<nrows; ++i) m[i]=new T[ncols];
    return m;
}

template <class T> inline void delete_matrix( T** m,int nrows) {
    int i=nrows;
    for (; i--;) delete [] m[i];
    delete [] m;
}

template <class T> inline
T*** new_tensor3(const int nrows,const int ncols,
                 const int npages) {
    int i,j;
    T*** t=new T**[nrows];
    for (i=0; i<nrows; ++i) {
        t[i]=new T*[ncols];
        for (j=0; j<ncols; ++j) t[i][j]=new T[npages];
    }
    return t;
}

template <class T> inline
void delete_tensor3( T*** t,int nrows,int ncols) {
    int i=nrows;
    for (; i--;) {
        for (int j=ncols; j--;) delete [] t[i][j];
        delete [] t[i];
    }
    delete [] t;
}

template <class T> inline
T**** new_tensor4(const int n1,const int n2,
                  const int n3,const int n4) {
    int i,j,k;
    T**** t=new T***[n1];
    for (int i=0; i<n1; ++i) {
        t[i]=new T**[n2];
        for (int j=0; j<n2; ++j) {
            t[i][j]=new T*[n3];
            for (k=0; k<n3; ++k) {
                t[i][j][k]=new T[n4];
            }
        }
    }
    return t;
}

template <class T> inline
void delete_tensor4( T**** t,int n1,int n2,int n3) {
    int i;
    for (int i=n1; i--;) {
        for (int j=n2; j--;) {
            for (int k=n3; k--;) delete [] t[i][j][k];
            delete [] t[i][j];
        }
        delete [] t[i];
    }
    delete [] t;
}

FILE* open_file(const char* name);

FILE* create_file(const char* name);

void make_filename(int fileno,int rank,const char* base,
                   char* filename);

void copy_identity (const int nr, const int nc, double *__restrict__ z, const int ldz);

void copy_trans (const int n, double* __restrict__ z);

constexpr double dist_sqr(const double *a,const double *b) noexcept {
    double tx=a[0]-b[0];
    double ty=a[1]-b[1];
    double tz=a[2]-b[2];
    return (tx*tx+ty*ty+tz*tz);
}

}
#endif
