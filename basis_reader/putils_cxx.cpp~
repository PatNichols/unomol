#include "putils_cxx.hpp"
namespace putils {

std::size_t tokenize_string(const std::string& str, const std::string& delims, std::vector<std::string>& tokens)
{
    std::size_t f;
    std::size_t s;
    tokens.clear();
    s = str.find_first_not_of(delims,0);
    f = str.find_first_of(delims,s);
    while (s!=std::string::npos) {
        tokens.push_back(str.substr(s,f-s));
        s = str.find_first_not_of(delims,f);
        f = str.find_first_of(delims,s);
    }
    return tokens.size();
}

std::string string_tolower(const std::string& str)
{
    std::string new_str(str);
    std::size_t sz = str.size();
    for (auto k=0;k<sz;++k) new_str[k] = tolower(new_str[k]);
    return new_str;    
}

std::string string_toupper(const std::string& str)
{
    std::string new_str(str);
    std::size_t sz = str.size();
    for (auto k=0;k<sz;++k) new_str[k] = toupper(new_str[k]);
    return new_str;    
}

std::string strings_join(int argc,char **argv)
{
    std::string jstr;
    if (!argc) return jstr;
    if ( argv == 0x0) {
        throw FormatError(__FUNCTION__," NULL STRING");
    }
    if ( argv[0] == 0x0) {
        throw FormatError(__FUNCTION__," NULL STRING");
    }
    jstr = argv[0];
    for ( auto i=1;i<argc;++i) {
        if ( argv[i] == 0x0) {
            throw FormatError(__FUNCTION__," NULL STRING");
        }
        jstr += " ";
        jstr += argv[i];
    }
    return jstr;
}

std::string strings_join( const std::vector<std::string>& strs)
{
    std::string jstr;
    std::size_t sz = strs.size();
    if (sz == 0) return jstr;
    jstr = strs[0];
    for ( auto i=1;i<sz;++i) {
        jstr += " ";
        jstr += strs[i];
    }
    return jstr;
}

bool string_compare_no_case( const std::string& s1, const std::string& s2)
{
    std::size_t sz = s1.size();
    if (sz != s2.size() ) return false;
    for (auto i=0;i<sz;++i) {
        if ( toupper(s1[i]) != toupper(s2[i]) ) return false;
    }
    return true;
}

template < class Tp > Tp string_to_type(const std::string& str)
{
    Tp x;
    try {
        std::istringstream in(str);
        in >> x;
    } catch (...) {
        std::cerr << " cannot convert " << str << " to specified type\n";
        throw FatalError(__FUNCTION__);
    }
    return x;
}

template <> bool string_to_type<bool>(const std::string& s)
{
    if ( s.size() == 0) throw FatalError(__FUNCTION__,"null string");
    if (s[0] == 'F' || s[0] == 'f' || s[0] == '0') return false;
    return true;
}

template <> std::string string_to_type<std::string>(const std::string& s)
{
    return std::string(s);
}

#define PUTILS_FUN(TYPE,FUN)\
template <> TYPE string_to_type<TYPE>(const std::string& s) {		\
    TYPE x;								\
    try {								\
        x = FUN/**/(s);							\
    } catch (...) {							\
        std::cerr << " cannot convert " << s << " to specified type\n";	\
        throw FatalError(__FUNCTION__);					\
    }									\
    return x;								\
}

PUTILS_FUN(long double,std::stold)
PUTILS_FUN(double,std::stod)
PUTILS_FUN(float,std::stof)
PUTILS_FUN(unsigned long long,std::stoull)
PUTILS_FUN(unsigned long,std::stoul)
PUTILS_FUN(long long,std::stoll)
PUTILS_FUN(long,std::stol)
PUTILS_FUN(int,std::stoi)
#undef PUTILS_FUN

#define PUTILS_FUN(TYPE,FUN)\
template <> TYPE string_to_type<TYPE>(const std::string& s) {		\
    TYPE x;								\
    try {								\
        x = static_cast</**/TYPE/**/>(/**/FUN/**/(s));			\
    } catch (...) {							\
        std::cerr << " cannot convert " << s << " to specified type\n";	\
        throw FatalError(__FUNCTION__);					\
    }									\
    return x;								\
}

PUTILS_FUN(unsigned int,std::stoul)
PUTILS_FUN(unsigned short,std::stoul)
PUTILS_FUN(unsigned char,std::stoul)
PUTILS_FUN(short,std::stoi)
PUTILS_FUN(char,std::stoi)
#undef PUTILS_FUN
} // end namespace std
