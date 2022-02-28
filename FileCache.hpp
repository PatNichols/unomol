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

#define DEF_CACHE_SIZE 1048576LL * 32LL
#define DEF_MAX_FSZ 1048576LL * 512LL

class FileCache
{
public:
    typedef std::int64_t size_type;

    FileCache() {
        prefix = "./tmp/tmp.";
        cache_size = DEF_CACHE_SIZE;
        max_file_size = DEF_MAX_FSZ;
        fp = nullptr;
        curr_pos = 0;
        curr_file = 0;
        mode = 0;
        file_sz.clear();
        cache = new char[cache_size];
        save_me = false;
        std::string dir = "./tmp";
        // create directory
        {
            namespace fs = std::filesystem;
            std::string dir("./tmp");
            fs::path p(dir);
            if (!fs::exists(p)) {
                fs::create_directory(p);
            }
        }
    }

    explicit FileCache(const std::string& dir_str,const std::string& name_str) {
        prefix = dir_str + std::string("/") + name_str + std::string(".");
        cache_size = DEF_CACHE_SIZE;
        max_file_size = DEF_MAX_FSZ;
        fp = nullptr;
        curr_pos = 0;
        curr_file = 0;
        mode = 0;
        file_sz.clear();
        cache = new char[cache_size];
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

    explicit FileCache(const std::string& dir_str,const std::string& name_str,size_type cache_size_,size_type file_size) {
        prefix = dir_str + std::string("/") + name_str + std::string(".");
        cache_size = cache_size_;
        max_file_size = file_size;
        fp = nullptr;
        curr_pos = 0;
        curr_file = 0;
        mode = 0;
        file_sz.clear();
        cache = new char[cache_size];
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

    explicit FileCache(size_type cache_size_,size_type file_size) {
        prefix = "./tmp/tmp.";
        cache_size = cache_size_;
        max_file_size = file_size;
        fp = nullptr;
        curr_pos = 0;
        curr_file = 0;
        mode = 0;
        file_sz.clear();
        cache = new char[cache_size];
        save_me = false;
        // create directory
        {
            namespace fs = std::filesystem;
            fs::path p(std::string("./tmp"));
            if (!fs::exists(p)) {
                fs::create_directory(p);
            }
        }
    }

    FileCache(const FileCache&) = delete;

    ~FileCache() {
        close();
        if (save_me) save();
        else {
            size_t nfile = file_sz.size();
            for (size_t k=1; k<nfile; ++k) {
                std::string fname = prefix + std::to_string(k);
                remove(fname.c_str());
            }
        }
        delete [] cache;
        cache = nullptr;
    }


    std::string get_prefix() const noexcept {
        return std::string(prefix);
    }

    void set_save() {
        save_me = true;
    }

    void set_prefix(const std::string& dir_str,const std::string& filename_str=std::string("tmp"))
    {
        prefix = dir_str + std::string("/") + filename_str + std::string(".");
    }

    void set_max_file_size(size_type max_size) {
        max_file_size = max_size;
    }

    size_type file_size(size_type i) const noexcept {
        return file_sz[i];
    }

    size_type num_files() const noexcept {
        return file_sz.size();
    }

    size_type total_size() const noexcept {
        size_type s = 0;
        for ( size_type k : file_sz ) {
            s += k;
        }
        return s;
    }

    constexpr size_type size_of_cache() const noexcept {
        return cache_size;
    }

    void save() {
        close();
        std::string fname = prefix + "meta";
        std::ofstream out(fname.c_str());
        out << file_sz.size() << "\n";
        for ( size_type s : file_sz ) {
            out << s << "\n";
        }
        out << cache_size << "\n";
        out << max_file_size << "\n";
        out.close();
        fname = prefix + "0";
        std::ofstream xout(fname.c_str());
        xout.write(cache,file_sz[0]);
        xout.close();
    }


    void load() {
        close();
        std::string fname = prefix + "meta";
        file_sz.clear();
        std::ifstream in(fname.c_str());
        size_type nf = 0;
        size_type s;
        in >> nf;
        for ( size_type k=0; k<nf; ++k ) {
            in >> s;
            file_sz.push_back(s);
        }
        size_type cache_size_in;
        in >>  cache_size_in;
        in >>  max_file_size;
        in.close();
        if (cache_size_in!=cache_size) {
            delete [] cache;
            cache_size = cache_size_in;
            cache = new char[cache_size];
        }
        fname = prefix + "0";
        std::ifstream xin(fname.c_str());
        xin.read(cache,file_sz[0]);
        xin.close();
        std::cout << " max file size = " << max_file_size << "\n";
        std::cout << " cache size    = " << cache_size << "\n";
        std::cout << " size 0        = " << file_sz[0] << "\n";
        mode = 0;
    }

    void close() {
        if (mode==0) return;
        if (curr_pos && mode==1) {
            file_sz.push_back(curr_pos);
        }
        if (fp) fclose(fp);
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
        fp = Fmemopen(cache,cache_size,"r");
        curr_file = 0;
        curr_pos = 0;
        mode = 2;
    }

    void open_for_writing() {
        if (mode!=0) close();
        file_sz.clear();
        fp = Fmemopen(cache,cache_size,"w");
        curr_file = 0;
        curr_pos = 0;
        mode = 1;
    }

    void open_for_appending() {
        if (mode!=0) close();
        constexpr size_type edge = 256LL;
        mode = 1;
        // find the last file that can be appended
        if ((file_sz[0]+edge)<cache_size) {
            fp = Fmemopen(cache,cache_size,"w");
            Fseek(fp,file_sz[0],SEEK_SET);
            return;
        }
        curr_file = 1;
        size_type cpos = cache_size;
        size_type nfiles = file_sz.size();
        for (size_type k=1; k<nfiles; ++k) {
            if ( (file_sz[k]+edge) < max_file_size) {
                curr_file = k;
                curr_pos = file_sz[k];
                break;
            }
        }
        std::string fname(prefix);
        fname += std::to_string(curr_file);
        fp = Fopen(fname.c_str(),"a");
    }

    void seek( size_type pos) {
        size_type cmode = mode;
        close();
        mode = cmode;
        if (pos < file_sz[0]) {
            if (mode==2) {
                fp = Fmemopen(cache,cache_size,"r");
                Fseek(fp,pos,SEEK_SET);
                curr_pos = pos;
                curr_file = 0;
                return;
            }
            if (mode==1) {
                fp = Fmemopen(cache,cache_size,"w");
                Fseek(fp,pos,SEEK_SET);
                curr_pos = pos;
                curr_file = 0;
                file_sz.clear();
            }
        }
        curr_file = 0;
        size_type cpos = file_sz[0];
        for (size_type k=1; k<file_sz.size(); ++k) {
            if ( (cpos + file_sz[k]) > pos ) {
                curr_pos = pos - cpos;
                curr_file = k;
                break;
            }
            cpos += file_sz[k];
        }
        if (curr_file==0) {
            std::cerr << "seek past end of data\n";
            std::cerr << "pos = " << pos << " size = " << cpos << "\n";
            exit(EXIT_FAILURE);
        }
        std::string fname(prefix);
        if (mode == 2) {
            fname += std::to_string(curr_file);
            fp = Fopen(fname.c_str(),"r");
            Fseek(fp,curr_pos,SEEK_SET);
        }
        if (mode == 1) {
            fname += std::to_string(curr_file);
            fp = Fopen(fname.c_str(),"w");
            Fseek(fp,curr_pos,SEEK_SET);
            file_sz.erase(file_sz.begin()+curr_file-1,file_sz.end());
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

    void read(void *ptr,size_type n) {
        if (n==0) return;
        if (ptr==nullptr) {
            std::cerr << "null pointer in read!\n";
            exit(EXIT_FAILURE);
        }
        size_type left = file_sz[curr_file] - curr_pos;
        if (left >= n ) {
            Fread(ptr,1,n,fp);
            curr_pos += n;
            return;
        }
        {
            size_type nr = Fread(ptr,1,left,fp);
            ++curr_file;
            fclose(fp);
            std::string fname(prefix);
            fname += std::to_string(curr_file);
            fp = Fopen(fname.c_str(),"r");
            curr_pos = 0;
            char * cptr = (char*)ptr;
            read((cptr+nr),(n-nr));
        }
    }

    template < class Tp > void read(Tp *ptr,size_type n) {
        if (n==0) return;
        if (ptr==nullptr) {
            std::cerr << "null pointer in read!\n";
            exit(EXIT_FAILURE);
        }
        size_type bytes_left = file_sz[curr_file] - curr_pos;
        size_type n_left = bytes_left/sizeof(Tp);
        if (n_left >= n ) {
            Fread(ptr,sizeof(Tp),n,fp);
            curr_pos += n * sizeof(Tp);
            return;
        }
        {
//            std::cerr << "opening next file \n";
//            std::cerr << "reading " << n_left << "\n";
//            std::cerr << "need    " << n << "\n";
            size_type nr = Fread(ptr,sizeof(Tp),n_left,fp);
            ++curr_file;
            fclose(fp);
            std::string fname(prefix);
            fname += std::to_string(curr_file);
            fp = Fopen(fname.c_str(),"r");
            curr_pos = 0;
            read((ptr+nr),(n-nr));
        }
    }

    void write(const void *ptr,size_type n) {
        if (n==0) return;
        if (ptr==nullptr) {
            std::cerr << "null pointer in write!\n";
            exit(EXIT_FAILURE);
        }
        size_type left = max_file_size - curr_pos;
        if (curr_file ==  0) left = cache_size - curr_pos;
        if (left >= n ) {
            Fwrite(ptr,1,n,fp);
            curr_pos += n;
            return;
        }
        {
            size_type nw = Fwrite(ptr,1,left,fp);
            file_sz.push_back((curr_pos+nw));
            ++curr_file;
            fclose(fp);
            std::string fname(prefix);
            fname += std::to_string(curr_file);
            fp = Fopen(fname.c_str(),"w");
            curr_pos = 0;
            const char * cptr = (char*)ptr;
            write((cptr+nw),(n-nw));
        }
    }

    template < class Tp > void write(const Tp *ptr,size_type n)
    {
        if (n==0) return;
        if (ptr==nullptr) {
            std::cerr << "null pointer in write!\n";
            exit(EXIT_FAILURE);
        }
        size_type bytes_left = max_file_size - curr_pos;
        if (curr_file ==  0) bytes_left = cache_size - curr_pos;
        size_type n_left = bytes_left/sizeof(Tp);
        if (n_left >= n ) {
            Fwrite(ptr,sizeof(Tp),n,fp);
            curr_pos += n * sizeof(Tp);
            return;
        }
        {
            size_type nw = Fwrite(ptr,sizeof(Tp),n_left,fp);
            file_sz.push_back((curr_pos+nw*sizeof(Tp)));
            ++curr_file;
            fclose(fp);
            std::string fname(prefix);
            fname += std::to_string(curr_file);
            fp = Fopen(fname.c_str(),"w");
            curr_pos = 0;
            write((ptr+nw),(n-nw));
        }
    }
private:
    char * cache;
    size_type cache_size;
    size_type curr_pos;
    size_type curr_file;
    size_type max_file_size;
    std::vector< size_type > file_sz;
    std::string prefix;
    FILE *fp;
    int mode;
    bool save_me;
};

}
#endif
