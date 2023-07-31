#ifndef PUTILS_C_H
#define PUTILS_C_H
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <errno.h>
#include <fcntl.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

void fatal_error();

void * Malloc(size_t n);
void * Grow(void **ptr,size_t old_size,size_t new_size);
void * Calloc(size_t n);
void * AlignedAlloc(size_t al,size_t n);

#define MALLOC_PTR(type) ( type * ) Malloc( sizeof ( type ) )
#define MALLOC_ARRAY(type,n) (type *)Malloc(sizeof(type)*(n))
#define CHAR_ALLOC(n) (char*)Malloc(n+1)

int Open(const char *name,int flgs);

FILE * Fdopen(int fdes,const char *mode);

FILE * Fopen(const char *name,const char *mode);

FILE * Fmemopen(char *b,size_t bsize,const char *mode);

FILE * Popen(const char *cmd,const char *mode);

void Write(int fd,const void *buff,size_t sz);

void Read(int fd,void *buff,size_t sz);

void BlockingRead(int fd,void *buff,size_t sz);

void Fwrite(const void *p,size_t osize,size_t cnt,FILE *fp);

void Fread(void *p,size_t osize,size_t cnt,FILE *fp);

void Fseek(FILE *fp,long pos,int whence);

void Pipe(int *fds);

int Fork();

#ifdef __cplusplus
}
#endif

#endif
