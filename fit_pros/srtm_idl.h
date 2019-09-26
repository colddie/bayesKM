#ifndef SRTM_IDL_H
#define SRTM_IDL_H

// #define EXPORT __declspec(dllexport)

class SRTM_IDL
{
private :

    const int parNr=3;
    int fitframeNr;
    double *t, *cr, *ct, *tis, *w; /* These are pointers, not allocated */
    double pmin[MAX_PARAMS], pmax[MAX_PARAMS];
    double wss_wo_penalty;
    /* Local functions */


public :
    extern "C"  int srtm_idl(int argc, char **argv);
    double srtmFunc(int parNr, double *p, void*);

};



#endif


