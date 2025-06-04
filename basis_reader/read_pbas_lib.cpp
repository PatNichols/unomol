
struct pShell
{
    int npr,lsh;
    std::vector<double> alf;
    std::vector<double> cof;

    pShell(int n,int l,const double *al,const double *co):
        npr(n),
        lsh(l),
        alf(al,n),
        cof(co,n)
    {}
};

struct BasisSet
{
    std::string atom;
    std::string name;
    std::vector<pShell> shells;
};

void read_basis_set(std::istream& is,
    BasisSet& b)
{
    b.name.clear();
    b.atom.clear();
    b.shells.clear();
    in >> str;
    size_t p = str.find(":");
    if ( p == std::string::npos) {
        format_error();
    }
    b.atom = str.substr(0,p);
    b.name = str.substr(p+1);
    in >> b.nsh;
    alf.clear();
    cof.clear();
    b.shells.resize(b.nsh);
    std::vector<double> alf;
    std::vector<double> cof;
    int npr,lsh;
    for (int j=0;j<b.nsh;++j)
    { 
        in >> npr;
        in >> lsh;
        alf.resize(npr);
        cof.resize(npr);
        for (int k=0;k<npr;++k) {
            in >> alf[k];
            in >> cof[k];
        }
        b.shells[j] = std::move(Shell(npr,lsh,alf.data(),cof.data());
    }
    return b;
}

bool read_atomic_sets(std::istream& in,
    BasisSet& bout, 
    const std::string& target_atom,
    const std::string& target_basis_name)
{    
    in >> atom;
    in >> nsets;
    if (atom.compare(target_atom))
    {
        for (int i=0;i<nsets;++i)
        {
            BasisSet b;
            read_basis_set(in,b);
        }
    }else{
        for (int i=0;i<nsets;++i)
        {
            BasisSet b;
            read_basis_set(in,b);
            if (!b.name.compare(target_basis_name))
            {
                bout = std::move(b);
                return true;
            } 
        }    
    }
    return false;
}


void ReadCenter()
{
}

void ReadInput()
{
    

}

