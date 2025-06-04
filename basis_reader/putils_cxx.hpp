#pragma once
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <stdexcept>
namespace putils
{
struct FileOpenError
{
    std::string msg;
    FileOpenError() = delete;
    FileOpenError(const char *where,const char *name,const char *mode):
        msg(where)
    {
        msg += " could not open the file ";
        msg += name;
        msg += " in mode ";
        msg += mode;
    }
    virtual const char * what() const noexcept
    {
        return msg.c_str();
    }
};

struct FormatError
{
    std::string msg;
    FormatError():msg(" unstated format error ") {}
    FormatError(const char *where,const char *name,const char *what_):
        msg(where)
    {
        msg += " format error in ";
        msg += name;
        msg += " ";
        msg += what_;
    }
    FormatError(const char *where,const char *what_):
        msg(where)
    {
        msg += " format error ";
        msg += what_;
    }
    virtual const char * what() const noexcept
    {
        return msg.c_str();
    }
};

struct FatalError
{
    std::string msg;
    FatalError():
        msg(" fatal error ") 
    {
    }
    explicit FatalError(const char *where, const char *what_):
        msg(where)
    {
        msg += " fatal error ";
        msg += what_;
        msg += " ";
        msg += strerror(errno);
    }

    explicit FatalError(const char *where):
        msg(where)
    {
        msg += " fatal error ";
        msg += strerror(errno);
    }
    virtual const char * what() const noexcept
    {
        return msg.c_str();
    }
};

struct  ReadError
{
    std::string msg;
    ReadError():msg(" istream read error ") {}

    explicit ReadError(const char *where,const char *name):
        msg(where)
    {
        msg = where;
        msg += " istream read error in ";
        msg += name;
    }
    
    explicit ReadError(const char *where):
        msg(where)
    {
        msg += " istream read error ";
    }
    virtual const char * what() const noexcept
    {
        return msg.c_str();
    }
};

struct OpenFileError: public std::exception
{
    std::string msg;
    explicit OpenFileError( const char * filename, const char *mode):
        msg("could not open the file ")
    {
        msg += filename;
        msg += " in mode ";
        msg += mode;
    }
    explicit OpenFileError( const char * where, const char * filename, const char *mode):
        msg(where)
    {
        msg += "could not open the file ";
        msg += filename;
        msg += " in mode ";
        msg += mode;
    }

    virtual const char * what() const noexcept
    {
        return msg.c_str();
    }
};

std::size_t tokenize_string(const std::string& str, const std::string& delims, std::vector<std::string>& tokens);
std::string string_toupper(const std::string& str);
std::string string_tolower(const std::string& str);
std::string strings_join(int argc,char **argv);
std::string strings_join( const std::vector<std::string>& strs);
bool string_compare_no_case( const std::string& s1, const std::string& s2);
template < class Tp > Tp string_to_type(const std::string& s);
void GetLine(std::istream& is,std::string str);
}
