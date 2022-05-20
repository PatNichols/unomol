
#include "Util.hpp"

namespace unomol {

void fatal_error(const char* message) {
    fprintf(stderr,"%s\n",message);
    exit(EXIT_FAILURE);
}

FILE* open_file(const char* name) {
    FILE* fp=fopen(name,"r");
    if (!fp) {
        fprintf(stderr,"Could not find %s\n",name);
        fatal_error("Does this file exist?");
    }
    return fp;
}

FILE* create_file(const char* name) {
    FILE* fp=fopen(name,"w");
    if (!fp) {
        fprintf(stderr,"Could not open %s\n",name);
        fatal_error("Too many open files?");
    }
    return fp;
}

void make_filename(int fileno,int rank,const char* base,
                   char* filename) {
    char prefix[12];
    sprintf(filename,"%04d_%03d_%s",fileno,rank,base);
}

void
copy_identity (const int nr, const int nc, double *__restrict__ z, const int ldz) {
    int i, j;
    double *zp;

    memset(z,0x0,sizeof(double)*nr*nc);
    for (i = 0; i < nr; i++) {
        (*(z + i * ldz + i)) = 1.0;
    }
}

void
copy_trans (const int n, double* __restrict__ z) {
    int i, j;
    double *row, *col, tmp;

    for (i = 0; i < n; i++) {
        row = z + i * n;
        col = z + i;
        for (j = 0; j < i; j++) {
            tmp = (*row);
            (*row) = (*col);
            (*col) = tmp;
            ++row;
            col += n;
        }
    }
}

}

