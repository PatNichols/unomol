#ifndef _STRUCTS_CC_
#define _STRUCTS_CC_

struct ShellPairData
{
    double ab2;
    const double *a,*b;
    int npr1,lv1;
    const double *al1,*co1;
    int npr2,lv2;
    const double *al2,*co2;
    unsigned short *lstates;
    unsigned int len;
};

struct MomInts
{
    double dx,dy,dz,qxx,qxy,qxz,qyy,qyz,qzz;
    unsigned int ijr;
};

struct Moments
{
    double dx,dy,dz,qxx,qxy,qxz,qyy,qyz,qzz;
    unsigned int ijr;
};

#endif

