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
#include <functional>
#include "putils_c.h"

namespace putils {

#define DEF_CACHE_SIZE 1048576LL * 32LL

class Cache
{
public:
    typedef std::int64_t size_type;

    Cache() {
        fname = "./tmp.dat";
        m_cache_size = DEF_CACHE_SIZE;
        fp = nullptr;
        curr_pos = 0;
        curr_file = 0;
        mode = 0;
        m_file_size = 0;
        m_mem_size = 0;
        cache = new char[m_cache_size+1];
        cache[m_cache_size] = '\0';
        save_me = false;
    }

    explicit Cache(const std::string& dir_str,const std::string& name_str) {
        fname = dir_str + std::string("/") + name_str;
        m_cache_size = DEF_CACHE_SIZE;
        fp = nullptr;
        curr_pos = 0;
        curr_file = 0;
        mode = 0;
        m_file_size = 0;
        m_mem_size = 0;
        cache = new char[m_cache_size+1];
        cache[m_cache_size] = '\0';
        save_me = false;
        // create directory
        {
            namespace fs = std::filesystem;
            fs::path p(dir_str);
            if (!fs::exists(p)) {
                fs::create_directory(p);
            }
        }
    }

    explicit Cache(const std::string& dir_str,const std::string& name_str,
        size_type cache_size_) {
        fname = dir_str + std::string("/") + name_str;
        m_cache_size = cache_size_;
        fp = nullptr;
        curr_pos = 0;
        curr_file = 0;
        mode = 0;
        m_mem_size = 0;
        m_file_size = 0;
        cache = new char[m_cache_size+1];
        cache[m_cache_size] = '\0';
        save_me = false;
        // create directory
        {
            namespace fs = std::filesystem;
            fs::path p(dir_str);
            if (!fs::exists(p)) {
                fs::create_directory(p);
            }
        }
    }

    explicit Cache(const std::string& name_str,
        size_type cache_size_) {
        fname = name_str;
        m_cache_size = cache_size_;
        fp = nullptr;
        curr_pos = 0;
        curr_file = 0;
        mode = 0;
        m_mem_size = 0;
        m_file_size = 0;
        cache = new char[m_cache_size+1];
        cache[m_cache_size] = '\0';
        save_me = false;
    }

    explicit Cache(size_type m_cache_size_) {
        fname = "tmp.dat";
        m_cache_size = m_cache_size_;
        fp = nullptr;
        curr_pos = 0;
        curr_file = 0;
        mode = 0;
        m_file_size = 0;
        m_mem_size = 0;
        cache = new char[m_cache_size+1];
        cache[m_cache_size] = '\0';
        save_me = false;
    }

    Cache(const Cache&) = delete;

    ~Cache() {
        close();
        if (save_me) {
            std::string mname = fname + ".meta";
            std::ofstream out(mname.c_str());
            out.write((const char*)&m_file_size,sizeof(int64_t));
            out.write((const char*)&m_mem_size,sizeof(int64_t));
            out.write((const char*)&m_cache_size,sizeof(int64_t));            
            out.write(cache,m_cache_size);
            out.close();
        } else {
            remove(fname.c_str());
        }
        delete [] cache;
        cache = nullptr;
    }
    
    size_t hash_cache() {
        return std::hash<std::string>{}(std::string(cache));
    }

    void load()
    {
            close();
            if ( cache ) {
                delete [] cache;
            } 
            std::string mname = fname + ".meta";
            std::ifstream in(mname.c_str());
            in.read((char*)&m_file_size,sizeof(int64_t));
            in.read((char*)&m_mem_size,sizeof(int64_t));
            in.read((char*)&m_cache_size,sizeof(int64_t));            
            cache = new char[m_cache_size+1];
            cache[m_cache_size]='\0';
            in.read(cache,m_cache_size);
            in.close();
            std::cout << "loaded\n";
            std::cout << "file size  = " << m_file_size << "\n";
            std::cout << "mem  size  = " << m_mem_size << "\n";
            std::cout << "cache size = " << m_cache_size << "\n";
            std::cout << "hash       = " << std::hash<std::string>{}(std::string(cache)) << "\n";
            std::cout << "\n";
            save_me = true;
            curr_file = 0;
            curr_pos = 0;
            mode = 0;
    }

    void info()
    {
            std::cout << "file size  = " << m_file_size << "\n";
            std::cout << "mem  size  = " << m_mem_size << "\n";
            std::cout << "cache size = " << m_cache_size << "\n";    
            std::cout << "hash       = " << std::hash<std::string>{}(std::string(cache)) << "\n";
    }

    std::string get_file_name() const noexcept {
        return std::string(fname);
    }

    void set_save(bool b = true) {
        save_me = b;
    }

    void set_fname(const std::string& dir_str,const std::string& filename_str=std::string("tmp"))
    {
        fname = dir_str + std::string("/") + filename_str;
    }

    size_type mem_size() const noexcept { 
        return m_mem_size;
    }
    
    size_type file_size() const noexcept {
        return m_file_size;
    }

    size_type total_size() const noexcept {
        return (m_file_size + m_mem_size);
    }

    constexpr size_type cache_size() const noexcept {
        return m_cache_size;
    }

    void close() {
        if (mode==0) return;
        if (fp) {
            fclose(fp);
            fp = nullptr;
        }
        if (curr_pos && mode==1) {
            if ( curr_file == 1 ) m_file_size = curr_pos;
            else m_mem_size = curr_pos;
        }
        mode = 0;
    }

    void open(const char *mode) {
        if (mode[0] == 'r') return open_for_reading();
        if (mode[0] == 'w') return open_for_writing();
        if (mode[0] == 'a') return open_for_appending();
        std::cerr << "illegal mode in open " << mode << "\n";
        exit(EXIT_FAILURE);
    }

    void open_for_reading() {
        if (mode!=0) close();
        fp = Fmemopen(cache,m_cache_size,"r");
        curr_file = 0;
        curr_pos = 0;
        mode = 2;
    }

    void open_for_writing() {
        if (mode!=0) close();
        fp = Fmemopen(cache,m_cache_size,"w");
        curr_file = 0;
        curr_pos = 0;
        m_file_size = 0;
        m_mem_size = 0;
        mode = 1;
    }

    void open_for_appending() {
        if (mode!=0) close();
        constexpr size_type edge = 256LL;
        if ( m_file_size!=0 ) {
            fp = Fopen(fname.c_str(),"a");
            Fseek(fp,m_file_size,SEEK_SET);
            curr_file = 1;
            curr_pos = m_file_size;
            mode = 1;
            m_file_size = 0;
        } else {
            fp = Fmemopen(cache,m_cache_size,"a");
            Fseek(fp,m_mem_size,SEEK_SET);
            curr_file = 0;
            curr_pos = m_mem_size;
            mode = 1;
            m_mem_size = 0; 
            m_file_size = 0;       
        }
    }

    template < class Tp > void read(Tp *ptr,size_type n) {
        if (n==0) return;
        if (ptr==nullptr) {
            std::cerr << "null pointer in read!\n";
            exit(EXIT_FAILURE);
        }
        if ( curr_file == 1 ) {
            if ( ( n * sizeof(Tp)) > m_file_size) {
                std::cerr << "Attempt to read past end of file\n";
                exit(-1);
            }
            Fread(ptr,sizeof(Tp),n,fp);
            curr_pos += n * sizeof(Tp);
            return;            
        }
        size_type n_left = (m_mem_size - curr_pos)/sizeof(Tp);
        if (n_left >= n ) {
            Fread(ptr,sizeof(Tp),n,fp);
            curr_pos += sizeof(Tp) * n;
            return;
        }
        if (n_left) {
            Fread(ptr,sizeof(Tp),n_left,fp);
            curr_pos += sizeof(Tp) * n_left;    
        }
        fclose(fp);
        fp = Fopen(fname.c_str(),"r");
        curr_pos = 0;
        curr_file = 1;
        return read((ptr+n_left),(n-n_left));
    }

    void write(const void *ptr,size_type n) {
        return write<char>(reinterpret_cast<const char*>(ptr),n);
    }

    void read(void *ptr,size_type n) {
        return read<char>(reinterpret_cast<char*>(ptr),n);
    }

    template < class Tp > void write(const Tp *ptr,size_type n)
    {
        if (n==0) return;
        if (ptr==nullptr) {
            std::cerr << "null pointer in write!\n";
            exit(EXIT_FAILURE);
        }
        if ( curr_file == 1) {
            Fwrite(ptr,sizeof(Tp),n,fp);
            curr_pos += n * sizeof(Tp);
            return;        
        }
        size_type n_left = (m_cache_size-curr_pos)/sizeof(Tp);
        if (n_left >= n ) {
            Fwrite(ptr,sizeof(Tp),n,fp);
            curr_pos += n * sizeof(Tp);
            return;
        }
        if ( n_left) 
        {
            Fwrite(ptr,sizeof(Tp),n_left,fp);
            curr_pos += n_left * sizeof(Tp);
        }
        m_mem_size = curr_pos;
        fclose(fp);        
        fp = Fopen(fname.c_str(),"w");
        curr_file = 1;
        curr_pos = 0;
        return write((ptr+n_left),(n-n_left));
    }

private:
    char * cache;
    size_type m_cache_size;
    size_type curr_pos;
    size_type curr_file;
    size_type m_mem_size;
    size_type m_file_size;
    std::string fname;
    FILE *fp;
    int mode;
    bool save_me;
};

}
#endif
