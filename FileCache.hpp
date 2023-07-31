#ifndef FILECACHE_HPP
#define FILECACHE_HPP
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <filesystem>

#include "putils_c.h"

namespace putils {

#define DEF_CACHE_SIZE 1048576LL * 64LL

class FileCache
{
public:
    typedef std::int64_t size_type;

    FileCache():
        fp(nullptr),
        m_file_name("./tmp.dat"),
        cache(nullptr),
        m_cache_size(DEF_CACHE_SIZE),
        curr_pos(0),
        curr_file(0),
        m_file_size(0),
        m_mem_size(0),
        mode(0),
        save_me(0)
    {
        cache = new char[m_cache_size+1];
    }

    explicit FileCache(const std::string& fname,size_type cache_size_):
        fp(nullptr),
        m_file_name(fname),
        cache(nullptr),
        m_cache_size(cache_size_),
        curr_pos(0),
        curr_file(0),
        m_file_size(0),
        m_mem_size(0),
        mode(0),
        save_me(0)
    {
        cache = new char[m_cache_size+1];
    }

    explicit FileCache(const std::string& fname):
        fp(nullptr),
        m_file_name(fname),
        cache(nullptr),
        m_cache_size(DEF_CACHE_SIZE),
        curr_pos(0),
        curr_file(0),
        m_file_size(0),
        m_mem_size(0),
        mode(0),
        save_me(0)
    {
            std::string lfname = m_file_name + ".meta";
            std::ifstream in(lfname.c_str());
            in.read((char*)&m_file_size,sizeof(size_type));
            in.read((char*)&m_cache_size,sizeof(size_type));
            in.read((char*)&m_mem_size,sizeof(size_type));
            cache = new char[m_cache_size+1];
            in.read(cache,m_mem_size);
            in.close();
            save_me = true;
    }

    explicit FileCache(size_type cache_size_):
        fp(nullptr),
        m_file_name("./tmp.dat"),
        cache(nullptr),
        m_cache_size(cache_size_),
        curr_pos(0),
        curr_file(0),
        m_file_size(0),
        m_mem_size(0),
        mode(0),
        save_me(0)
    {
        cache = new char[m_cache_size+1];
    }

    FileCache(const FileCache&) = delete;

    ~FileCache() {
        close();
        if (save_me) {
            save();
        } else {
            remove(m_file_name.c_str());
        }
        delete [] cache;
        cache = nullptr;
    }


    std::string get_file_name() const noexcept {
        return std::string(m_file_name);
    }

    void set_save(bool b=true) {
        save_me = b;
    }

    size_type file_size(size_type i) const noexcept {
        return m_file_size;
    }

    size_type total_size() const noexcept {
        return m_mem_size + m_file_size;
    }

    constexpr size_type max_size_of_cache() const noexcept {
        return m_cache_size;
    }
    
    constexpr size_type memory_size() const noexcept {
        return m_mem_size;
    }

    void save() {
        close();
        std::string fname = m_file_name + ".meta";
        std::ofstream out(fname.c_str());
        out.write((const char*)&m_file_size,sizeof(size_type));
        out.write((const char*)&m_cache_size,sizeof(size_type));
        out.write((const char*)&m_mem_size,sizeof(size_type));
        out.write(cache,m_mem_size);
        out.close();
    }


    void load() {
        close();
        if (cache) delete [] cache;
        std::string fname = m_file_name + ".meta";
        std::ifstream in(fname.c_str());
        in.read((char*)&m_file_size,sizeof(size_type));
        in.read((char*)&m_cache_size,sizeof(size_type));
        in.read((char*)&m_mem_size,sizeof(size_type));
        cache = new char[m_cache_size+1];
        in.read(cache,m_mem_size);
        in.close();
        save_me = true;
    }

    void close() {
        if (mode==0) return;
        if (mode==1) {
            if ( curr_file == 1 ) {
                m_file_size = curr_pos;
            } else {
                m_mem_size = curr_pos;
            }
        }
        if (fp) fclose(fp);
        mode = 0;
        fp = nullptr;
    }

    void open(const char *mode) {
        if (mode[0] == 'r') return open_for_reading();
        if (mode[0] == 'w') return open_for_writing();
        if (mode[0] == 'a') return open_for_appending();
        std::cerr << "illegal mode in open " << mode << "\n";
        exit(EXIT_FAILURE);
    }

    void open_for_reading() {
        close();
        fp = Fmemopen(cache,m_cache_size,"r");
        curr_file = 0;
        curr_pos = 0;
        mode = 2;
    }

    void open_for_writing() {
        close();
        fp = Fmemopen(cache,m_cache_size,"w");
        curr_file = 0;
        curr_pos = 0;
        mode = 1;
        m_mem_size = 0;
        m_file_size = 0;
    }

    void open_for_appending() {
        close();
        mode = 1;
        if ( m_file_size ) {
            fp = Fopen(m_file_name.c_str(),"a");
            curr_file = 1;
            curr_pos = m_file_size;
            return;
        }
        if ( m_mem_size == 0) {
            open_for_writing();
            return;
        }
        fp = Fmemopen(cache,m_cache_size,"w");
        Fseek(fp,m_mem_size,SEEK_SET);
        curr_file = 0;
        curr_pos = m_mem_size;
    }

    void seek( size_type pos) {
        size_type cmode = mode;
        close();
        mode = cmode;
        if (mode == 0) {
            std::cerr << "seek on closed file system not possible\n";
            exit(EXIT_FAILURE);
        } 
        if (pos > (m_file_size+m_mem_size)) {
            std::cerr << "seek past end of data\n";
            std::cerr << "pos = " << pos << " size = " << (m_mem_size + m_file_size) << "\n";
            exit(EXIT_FAILURE);            
        }
        if (pos==0) {
            if (mode==1) open_for_reading();
            if (mode==2) open_for_writing();
            return;
        }
        if (pos <= m_mem_size) {
            curr_file = 0;
            curr_pos = pos;
            if (mode==1) {
                fp = Fmemopen(cache,m_cache_size,"w");
                m_mem_size = pos;
                m_file_size = 0;            
            }else{
                fp = Fmemopen(cache,m_cache_size,"r");
            }  
            Fseek(fp,pos,SEEK_SET);
            return;
        }else{
            curr_file = 1;
            curr_pos = pos;
            if (mode==1) {
                fp = Fopen(m_file_name.c_str(),"w");
                m_file_size = pos;
            }else{
                fp = Fopen(m_file_name.c_str(),"r");
            }  
            Fseek(fp,pos,SEEK_SET);
            return;
        }
    }

    void rewind() {
        int cmode = mode;
        close();
        mode = cmode;
        if (mode==1) {
            open_for_writing();
        } else {
            open_for_reading();
        }
    }

    template < class Tp >
    void read(Tp *ptr,size_t n) {
        if (n==0) return;
        if (ptr==nullptr) {
            std::cerr << "null pointer in read!\n";
            exit(EXIT_FAILURE);
        }
        if (curr_file == 1) {
            Fread((void*)ptr,sizeof(Tp),n,fp);
            curr_pos += n * sizeof(Tp);
            return;
        }
        size_type left =  (m_mem_size - curr_pos)/sizeof(Tp);
        if (left >= n ) {
            Fread((void*)ptr,sizeof(Tp),n,fp);
            curr_pos += n * sizeof(Tp);
            return;
        }
        if ( left ) {
            Fread((void*)ptr,sizeof(Tp),size_t(left),fp);
        }
        fclose(fp);
        curr_file = 1;
        curr_pos = 0;
        ptr += left;
        n -= left;
        fp = Fopen(m_file_name.c_str(),"r");
        Fread((void*)ptr,sizeof(Tp),n,fp);
    }
    
    template < class Tp > 
    void write(const Tp *ptr,size_type n) {
        if (n==0) return;
        if (ptr==nullptr) {
           std::cerr << "null pointer in write!\n";
           exit(EXIT_FAILURE);
       }        
        if (curr_file == 1)
        {
            Fwrite((const void*)ptr,sizeof(Tp),n,fp);
            curr_pos += n * sizeof(Tp);
            return;
        }
        size_type left = (m_cache_size - curr_pos)/sizeof(Tp);
         if (left >= n ) {
            Fwrite((const void*)ptr,sizeof(Tp),n,fp);
            curr_pos += n * sizeof(Tp);
            return;
        }
        if ( left ) {
            Fwrite((const void*)ptr,sizeof(Tp),size_t(left),fp);
            curr_pos += left * sizeof(Tp);
        }
        m_mem_size = curr_pos;
        curr_file = 1;
        curr_pos = 0;
        fclose(fp);
        fp = Fopen(m_file_name.c_str(),"w");
        ptr += left;
        n -= left;
        Fwrite((const void*)ptr,sizeof(Tp),n,fp);
    }         

    void write(const void *ptr,size_t n) {
        return write(reinterpret_cast<const char*>(ptr),n);
    }

    void read(void *ptr,size_t n) {
        return reade(reinterpret_cast<char*>(ptr),n);
    }

private:
    FILE * fp;
    std::string m_file_name;
    char * cache; // the in-memory cache
    size_type m_cache_size; // max number of bytes that fit in memory cache
    size_type curr_pos; 
    size_type curr_file;
    size_type m_file_size; // actual number of bytes in file
    size_type m_mem_size;  // actual number of bytes in memory cache
    int mode;
    bool save_me;
};

}
#endif
