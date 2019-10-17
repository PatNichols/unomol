#ifndef _UTIL_CC_
#define _UTIL_CC_
#include <cstdio>
#include <cstdlib>
#include <cstring>
#ifdef _MPI_ABI_
#include <mpi.h>
#endif
#include "pconfig.h"
using namespace std;

inline void fatal_error(const char* message)
{
    fprintf(stderr,"%s\n",message);
#ifdef _MPI_API_
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
#else
    exit(EXIT_FAILURE);
#endif
}

template <class T>
inline T** new_matrix(const int nrows,const int ncols)
{
    register int i;
    T** m=new T*[nrows];
    for (i=0;i<nrows;++i) m[i]=new T[ncols];
    return m;
}

template <class T>
inline void delete_matrix( T** m,int nrows)
{
    register int i=nrows;
    for (;i--;) delete [] m[i];
    delete [] m;
}

template <class T>
inline T*** new_tensor3(const int nrows,const int ncols,
                        const int npages)
{
    register int i,j;
    T*** t=new T**[nrows];
    for (i=0;i<nrows;++i)
    {
        t[i]=new T*[ncols];
        for (j=0;j<ncols;++j) t[i][j]=new T[npages];
    }
    return t;
}

template <class T>
inline void delete_tensor3( T*** t,int nrows,int ncols)
{
    register int i=nrows;
    for (;i--;)
    {
        for (register int j=ncols;j--;) delete [] t[i][j];
        delete [] t[i];
    }
    delete [] t;
}

template <class T>
inline T**** new_tensor4(const int n1,const int n2,
                         const int n3,const int n4)
{
    register int i,j,k;
    T**** t=new T***[n1];
    for (i=0;i<n1;++i)
    {
        t[i]=new T**[n2];
        for (j=0;j<n2;++j)
        {
            t[i][j]=new T*[n3];
            for (k=0;k<n3;++k)
            {
                t[i][j][k]=new T[n4];
            }
        }
    }
    return t;
}

template <class T>
inline void delete_tensor4( T**** t,int n1,int n2,int n3)
{
    register int i;
    for (register int i=n1;i--;)
    {
        for (register int j=n2;j--;)
        {
            for (register int k=n3;k--;) delete [] t[i][j][k];
            delete [] t[i][j];
        }
        delete [] t[i];
    }
    delete [] t;
}

inline FILE* open_file(const char* name)
{
    FILE* fp=fopen(name,"r");
    if (!fp) {
        fprintf(stderr,"Could not find %s\n",name);
        fatal_error("Does this file exist?");
    }
    return fp;
}

inline FILE* create_file(const char* name)
{
    FILE* fp=fopen(name,"w");
    if (!fp) {
        fprintf(stderr,"Could not open %s\n",name);
        fatal_error("Too many open files?");
    }
    return fp;
}

inline void make_filename(int fileno,int rank,const char* base,
                          char* filename)
{
    char prefix[9];
    char to_text[]={ '0','1','2','3','4','5','6','7','8','9' };
#ifdef _USE_TMP_
    const char tmpdir[]="/tmp/";
    int tmpsize=6;
#elif defined _USE_SCRATCH_
    const char tmpdir[]="/scratch/";
    int tmpsize=10;
#else
    const char tmpdir[]="./";
    int tmpsize=3;
#endif
    size_t baselen = strlen(base) + 1; 
    if (fileno>9999) {
        fatal_error("filenumber>9999 in make_filename\n");
    }
    if (rank>99) {
        fatal_error("rank>99 in make_filename\n");
    }
    strncpy(filename,tmpdir,tmpsize);
    prefix[0]=to_text[fileno/1000];
    prefix[1]=to_text[fileno/100%10];
    prefix[2]=to_text[fileno/10%10];
    prefix[3]=to_text[fileno%10];
    prefix[4]='.';
    prefix[5]=to_text[rank/10];
    prefix[6]=to_text[rank%10];
    prefix[7]='.';
    prefix[8]='\0';
    strncat(filename,prefix,9);
    strncat(filename,base,baselen);
}

inline double dist_sqr(const double* a,const double* b)
{
    register double tx,ty,tz;
    tx=a[0]-b[0];
    ty=a[1]-b[1];
    tz=a[2]-b[2];
    return (tx*tx+ty*ty+tz*tz);
}


inline void
copy_identity (const int nr, const int nc, double *__restrict__ z, const int ldz)
{
    register int i, j;
    double *zp;

    for (i = 0; i < nr; i++)
    {
        zp = z + i * ldz;
        for (j = 0; j < nc; j++)
        {
            (*zp++) = 0.0;
        }
        (*(z + i * ldz + i)) = 1.0;
    }
}

inline void
copy_trans (const int n, double* __restrict__ z)
{
    int i, j;
    double *row, *col, tmp;

    for (i = 0; i < n; i++)
    {
        row = z + i * n;
        col = z + i;
        for (j = 0; j < i; j++)
        {
            tmp = (*row);
            (*row) = (*col);
            (*col) = tmp;
            ++row;
            col += n;
        }
    }
}

#endif
