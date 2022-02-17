
#include "Util.hpp"

namespace unomol {

void fatal_error(const char* message) {
    fprintf(stderr,"%s\n",message);
#ifdef _MPI_API_
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
#else
    exit(EXIT_FAILURE);
#endif
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
    char to_text[]= { '0','1','2','3','4','5','6','7','8','9' };
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
    prefix[4]='_';
    prefix[5]=to_text[rank/100];
    prefix[6]=to_text[rank/10%10];
    prefix[7]=to_text[rank%10];
    prefix[8]='_';
    prefix[9]='\0';
    strncat(filename,prefix,10);
    strncat(filename,base,baselen);
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

