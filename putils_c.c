#include "putils_c.h"

void default_exit_function() { exit(EXIT_FAILURE);}

void (*exit_fun)() = default_exit_function;

#ifdef __cplusplus
extern "C" {
#endif

void set_exit_function(void (*f)()) {
    exit_fun = f;
}

void fatal_error() { 
    (*exit_fun)();
}

void * Malloc(size_t n) {
    void * p = malloc(n);
    if (p) return p;
    fprintf(stderr,"malloc failed for size %lu\n",n);
    fatal_error();
    return 0x0;
}

void * Calloc(size_t n) {
    void * p = Malloc(n);
    memset(p,0x0,n);
    return p;
}

void * Resize(void **p,size_t old_size,size_t new_size)
{
    size_t cpy_size;
    void *tmp;

    if (old_size==new_size) return p;
    tmp = Malloc(new_size);
    cpy_size = (new_size > old_size) ? old_size:new_size;
    memcpy(tmp,*p,cpy_size);
    free(*p);
    *p = tmp;
    return *p;
}

FILE * Fopen(const char *name,const char *mode) {
    FILE *fp = fopen(name,mode);
    if ( !fp ) {
        fprintf(stderr,"could not open file %s in mode %s\n",name,mode);
        fprintf(stderr,"%s\n",strerror(errno));
        fatal_error();
    }
    return fp;
}

int Open(const char *name,int flgs)
{
    int f = open(name,flgs);
    if (f>=0) return f;
    fprintf(stderr,"error in opening %s : %s\n",name,strerror(errno));
    fatal_error();
    return -1;
}

FILE *Fdopen(int desc,const char *mode)
{
    FILE *fp = fdopen(desc,mode);
    if (fp) return fp;
    fprintf(stderr,"fdopen failed for desc %d in mode %s error = %s\n",
            desc,mode,strerror(errno));
    fatal_error();
    return 0x0;
}

FILE * Fmemopen(char *b,size_t bsize,const char *mode)
{
    FILE * fp = fmemopen(b,bsize,mode);
    if (fp) return fp;
    fprintf(stderr,"error in fmemopen %s\n",strerror(errno));
    fatal_error();
    return 0x0;
}

void Read(int fd,void *buff,size_t sz)
{
    ssize_t c = 0;
    if (sz==0) {
        fprintf(stderr,"warning size = 0 in read!\n");
        return;
    }
    c = read(fd,buff,sz);
    if (c < sz) {
        fprintf(stderr,"read failed %s\n",strerror(errno));
        fatal_error();
    }
    return;
}


void Write(int fd,const void *buff,size_t sz)
{
    ssize_t e = 0;
    ssize_t c = 0;
    ssize_t sz0 = sz;
    ssize_t rem = sz0;
    const char * bp = buff;
    if (sz==0) {
        fprintf(stderr,"warning sz = 0 in write!\n");
        return;
    }
    while (1) {
        c = write(fd,bp,rem);
        if (c < 0) {
            fprintf(stderr,"block write failed %s\n",strerror(errno));
            fatal_error();
        }
        rem -= c;
        if (rem==0) break;
        bp += c;
    }
    return;
}

void BlockingRead(int fd,void *buff,size_t sz)
{
    struct timespec tdurr;
    struct timespec twait;
    double acc = 0.0;
    double tmax = 20.0;
    double dwait = 100000*1.e-9;
    twait.tv_sec = 0;
    twait.tv_nsec = 100000;
    ssize_t rem = sz;
    char * bp = buff;
    if (sz==0) {
        fprintf(stderr,"warning size = 0 in read!\n");
        return;
    }
    while (1) {
        ssize_t c = read(fd,(void*)bp,rem);
        fprintf(stderr,"c readblock read size %ld\n",c);
        if (c<0) {
            fprintf(stderr,"blocking read failed : %s\n",strerror(errno));
            fatal_error();
        }
        if (c==0) {
            int e = nanosleep(&twait,&tdurr);
            if (e < 0) {
                if (errno == EINTR) continue;
                fprintf(stderr,"error in nanosleep %s \n",strerror(errno));
                fatal_error();
            }
            acc += dwait;
            if ( acc > tmax) {
                fprintf(stderr,"block read timed out!\n");
                fatal_error();
            }
        } else {
            rem -= c;
            if (rem==0) break;
            bp += c;
        }
    }
    return;
}


void Fwrite(const void *p,size_t osize,size_t cnt,FILE *fp)
{
    int64_t e = fwrite(p,osize,cnt,fp);
    if (e==cnt) return;
    e = ferror(fp);
    fprintf(stderr,"Fwrite failed ferror = %s\n",strerror(e));
    fatal_error();
}

void Fread(void *p,size_t osize,size_t cnt,FILE *fp)
{
    int64_t e = fread(p,osize,cnt,fp);
    if (e==cnt) return;
    e = ferror(fp);
    fprintf(stderr,"Fread failed ferror = %s\n",strerror(e));
    fatal_error();
}

void Fseek(FILE *fp,long pos,int whence)
{
    int e = fseek(fp,pos,whence);
    if (!e) return;
    fprintf(stderr,"error in fseek %s\n",strerror(errno));
    fatal_error();
}

int Fork()
{
    int p = fork();
    if ( p < 0) {
        fprintf(stderr,"fork failed %s!\n",strerror(errno));
        fatal_error();
    }
    return p;
}

void Pipe(int *fds)
{
    int e = pipe(fds);
    if (e) {
        fprintf(stderr,"opening pipe failed %s!\n",strerror(errno));
        fatal_error();
    }
}

FILE * Popen(const char *cmd,const char *mode)
{
    FILE * fp = popen(cmd,mode);
    if ( !fp ) {
        fprintf(stderr,"popen failed %s mode %s : %s\n",
            cmd,mode,strerror(errno));
        fatal_error();
    }
    return fp;
}


#ifdef __cplusplus
}
#endif
