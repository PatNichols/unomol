#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>

typedef struct
{
    double v;
    unsigned short int i,j,k,l;
} teststruct;

void maxfilesize()
{

    size_t ts=sizeof(teststruct);
    size_t pgsz=get_pagesize();
}



int main()
{
    FILE *fp,*fpd;
    DIR *dirp;

    fp=fopen("pconfig.h","w");
    if (!fp) {
        fprintf(stderr,"Error could not open pconfig.h\n");
        fprintf(stderr,"Check Read Write permissions in this directory\n");
        fprintf(stderr,"Check for too many open files\n");
        exit(EXIT_SUCCESS);
    }
    fprintf(fp,"/***** Configuration FILE for UNOMOL  *****/");
    fprintf(fp,"#ifndef _PCONFIG_H_\n");
    fprintf(fp,"#define _PCONFIG_H_\n");
    int lsize=sizeof(unsigned long);
    if (lsize!=4) {
        fprintf(fp,"#define _64BIT_API_\n");
    }
    dirp=opendir("/tmp");
    if (dirp!=0) {
        fpd=fopen("/tmp/somefile","w");
        if (fpd!=0) {
            fclose(fpd);
            fprintf(fp,"#endif\n");
            fclose(fp);
            closedir(dirp);
            return EXIT_SUCCESS;
        }
    }
    dirp=opendir("/scratch");
    if (dirp!=0) {
        fpd=fopen("/scratch/somefile","w");
        if (fpd!=0) {
            fprintf(fp,"#define _USE_SCRATCH\n");
            fclose(fpd);
            fprintf(fp,"#endif\n");
            fclose(fp);
            closedir(dirp);
            return EXIT_SUCCESS;
        }
    }
    fprintf(fp,"#define _USE_PWD\n");
    fprintf(fp,"#endif\n");


    fclose(fp);
    return EXIT_SUCCESS;
}
