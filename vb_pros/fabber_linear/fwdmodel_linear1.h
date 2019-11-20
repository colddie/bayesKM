// fwdmodel_sine.h - A simple sine curve fitting model
// First we will create the interface fwdmodel_exp.h file which shows the methods we will need to implement:
// We have not made our methods virtual, so nobody will be able to create a subclass of our model. If we wanted 
// this to be the case all the non-static methods would need to be virtual, and we would need to add a virtual destructor. 
// This is sometimes useful when you want to create variations on a basic model.
// Most of the code above is completely generic to any model. The only parts which are specific to our exp-function model are:
// The name ExpFwdModel
// The private variables m_num (the number of exponentials in our sum) and m_dt (the time between data points).

#ifndef FWDMODEL_LINEAR1_H
#define FWDMODEL_LINEAR1_H

#include "fabber_core/fwdmodel.h"

#include "newmat.h"
// #include "tpccm.h"
// #include "sim1cm.c"

#include <string>

class PetFwdModel : public FwdModel {
public:
    static FwdModel* NewInstance();

    PetFwdModel()
        : m_include_offset(false)
    {
    }

    void GetOptions(std::vector<OptionSpec>& opts) const;
    std::string GetDescription() const;
    std::string ModelVersion() const;

    void Initialize(FabberRunData& args);
    int NumParams() const;
    void NameParams(std::vector<std::string>& names) const;
    void HardcodedInitialDists(MVNDist& prior, MVNDist& posterior) const;
    void Evaluate(const NEWMAT::ColumnVector& params, NEWMAT::ColumnVector& result) const;

private:
    bool m_include_offset;
    static FactoryRegistration<FwdModelFactory,PetFwdModel> registration;

protected:
    // NEWMAT::ColumnVector plasma_c;
    // NEWMAT::ColumnVector plasma_t;
    // double *plasma_c;
    double *plasma_t;
    int     nsample;
    double  pmean1;
    double  pmean2;
    double  pprecision1;
    double  pprecision2;    

    // bool  usepriorimg;
    //NEWMAT::RowVector pmeanimg1;
    //NEWMAT::RowVector pmeanimg2;
};

#endif
