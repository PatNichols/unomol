#include "putils_c.h"

#ifdef __cplusplus
extern "C" {
#endif

    void * Malloc(size_t n) {
        void * p = malloc(n);
        if (p) return p;
        fprintf(stderr,"malloc failed for size %lu\n",n);
        exit(EXIT_FAILURE);
    }

    void * Calloc(size_t n) {
        void * p = Malloc(n);
        memset(p,0x0,n);
        return p;
    }

    void * Grow(void **p,size_t old_size,size_t new_size)
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
        if (fp) return fp;
        fprintf(stderr,"could not open file %s in mode %s\n",name,mode);
        exit(EXIT_FAILURE);
    }

    int Open(const char *name,int flgs)
    {
        int f = open(name,flgs);
        if (f>=0) return f;
        fprintf(stderr,"error in opening %s : %s\n",name,strerror(errno));
        exit(EXIT_FAILURE);
    }

    FILE *Fdopen(int desc,const char *mode)
    {
        FILE *fp = fdopen(desc,mode);
        if (fp) return fp;
        fprintf(stderr,"fdopen failed for desc %d in mode %s error = %s\n",
                desc,mode,strerror(errno));
        exit(EXIT_FAILURE);
    }

    FILE * Fmemopen(char *b,size_t bsize,const char *mode)
    {
        FILE * fp = fmemopen(b,bsize,mode);
        if (fp) return fp;
        fprintf(stderr,"error in fmemopen %s\n",strerror(errno));
        exit(EXIT_FAILURE);
    }

    void Write(int fd,const void *buff,size_t sz)
    {
        size_t e = 0;
        ssize_t c = 0;
        ssize_t rem = sz;
        const void * bp = buff;
        if (sz==0) {
            fprintf(stderr,"size = 0 in write!\n");
            return;
        }
        while (1) {
            c = write(fd,bp,rem);
            if (c < 0) {
                fprintf(stderr,"write failed  %s\n",strerror(errno));
                exit(EXIT_FAILURE);
            }
            rem -= c;
            if (rem==0) break;
            bp += c;
        }
        return;
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
            exit(EXIT_FAILURE);
        }
        return;
    }


    void BlockingWrite(int fd,const void *buff,size_t sz)
    {
        ssize_t e = 0;
        ssize_t c = 0;
        ssize_t sz0 = sz;
        ssize_t rem = sz0;
        const void * bp = buff;
        if (sz==0) {
            fprintf(stderr,"sz = 0 in write!\n");
            return;
        }
        while (1) {
            c = write(fd,bp,rem);
            if (c < 0) {
                fprintf(stderr,"block write failed %s\n",strerror(errno));
                exit(EXIT_FAILURE);
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

        twait.tv_sec = 0;
        twait.tv_nsec = 100000;
        ssize_t rem = sz;
        void * bp = buff;
        if (sz==0) {
            fprintf(stderr,"size = 0 in read!\n");
            return;
        }
        while (1) {
            ssize_t c = read(fd,bp,rem);
            fprintf(stderr,"c readblock read size %ld\n",c);
            if (c<0) {
                fprintf(stderr,"blocking read failed : %s\n",strerror(errno));
                exit(EXIT_FAILURE);
            }
            if (c==0) {
                tdurr.tv_sec = 0;
                tdurr.tv_nsec = 100000;
                int e = nanosleep(&twait,&tdurr);
                if (e < 0) {
                    if (errno == EINTR) return;
                    fprintf(stderr,"error in nanosleep %s \n",strerror(errno));
                    exit(EXIT_FAILURE);
                }
                fprintf(stderr,"c wait\n");
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
        exit(EXIT_FAILURE);
    }

    void Fread(void *p,size_t osize,size_t cnt,FILE *fp)
    {
        int64_t e = fread(p,osize,cnt,fp);
        if (e==cnt) return;
        e = ferror(fp);
        fprintf(stderr,"Fread failed ferror = %s\n",strerror(e));
        exit(EXIT_FAILURE);
    }

    /***
     *  This loops waiting for the required size to output to file
     *  note::: this can be dangerous!
     ***/
    void BlockingFread(void *buff,size_t osize,size_t cnt,FILE *fp)
    {
        struct timespec tdurr;
        struct timespec twait;

        twait.tv_sec = 0;
        twait.tv_nsec = 100000;
        tdurr.tv_sec = 0;
        tdurr.tv_nsec = 0;
        size_t e = 0;
        size_t c = 0;
        size_t rem = cnt;
        void * bp = buff;
        if (cnt==0) {
            fprintf(stderr,"WARNING : cnt = 0 in Fread!\n");
            return;
        }
        while (1) {
            if (feof(fp)) {
                e = nanosleep(&twait,&tdurr);
            }
            c = fread(bp,osize,rem,fp);
            if (c!=rem) {
                tdurr.tv_nsec = 0;
                e = nanosleep(&twait,&tdurr);
            } else {
                rem -= c;
                if (rem==0) break;
                bp += c;
            }
        }
        return;
    }


    void Fseek(FILE *fp,long pos,int whence)
    {
        int e = fseek(fp,pos,whence);
        if (!e) return;
        fprintf(stderr,"error in fseek %s\n",strerror(errno));
        exit(EXIT_FAILURE);
    }

    int Fork()
    {
        int p = fork();
        if (p>=0) return p;
        fprintf(stderr,"fork failed!\n");
        exit(EXIT_FAILURE);
    }

    void Pipe(int *fds)
    {
        int e = pipe(fds);
        if (!e) return;
        fprintf(stderr,"opening pipe failed %s!\n",strerror(errno));
        exit(EXIT_FAILURE);
    }

#ifdef __cplusplus
}
#endif
