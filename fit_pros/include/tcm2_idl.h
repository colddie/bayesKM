#ifndef TCM2_IDL_H
#define TCM2_IDL_H

// #define EXPORT __declspec(dllexport)

class TCM2_IDL
{
private :
    const int parNr=4;
    DFT input, data;
    double *petmeas, *petsim;
    double fVb;
    double pmin[MAX_PARAMETERS], pmax[MAX_PARAMETERS];
    double fk1k2;
    int fitframeNr;
    double wss_wo_penalty;

public :
    extern "C" int tcm2_idl(int argc, char **argv);
    double cm3Func(int parNr, double *p, void*);

};



#endif


