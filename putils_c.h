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

#ifdef __cplusplus
extern "C" {
#endif

#define BUFF_SIZE 64

void * Malloc(size_t n);
void * Grow(void **ptr,size_t old_size,size_t new_size);
void * Calloc(size_t n);

#define MALLOC_PTR(type) ( type * ) Malloc( sizeof ( type ) )
#define MALLOC_ARRAY(type,n) (type *)Malloc(sizeof(type)*(n))
#define CHAR_ALLOC(n) (char*)Malloc(n)

int Open(const char *name,int flgs);

FILE * Fdopen(int fdes,const char *mode);
FILE * Fopen(const char *name,const char *mode);
FILE * Fmemopen(char *b,size_t bsize,const char *mode);

ssize_t Write(int fd,const void *buff,size_t sz);
ssize_t Read(int fd,void *buff,size_t sz);
ssize_t BlockingRead(int fd,void *buff,size_t sz);

size_t Fwrite(const void *p,size_t osize,size_t cnt,FILE *fp);
size_t Fread(void *p,size_t osize,size_t cnt,FILE *fp);
size_t BlockingFread(void *ptr,size_t osize,size_t cnt,FILE *fp);

void Fseek(FILE *fp,long pos,int whence);

void Pipe(int *fds);
int Fork();

#ifdef __cplusplus
}
#endif

#endif
