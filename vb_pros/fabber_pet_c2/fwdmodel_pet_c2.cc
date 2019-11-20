//  fwdmodel_sine.cc - Implements a simple sine curve fitting model


#include "fwdmodel_pet_c2.h"

#include "fabber_core/fwdmodel.h"

#include "dist_mvn.h"
#include "easylog.h"
#include "rundata.h"
#include "tools.h"
#include "version.h"

#include <newmatio.h>
#include <math.h>

using namespace std;
using namespace NEWMAT;

// We need to implement a couple of methods to ensure that our model is visible to the Fabber system:
// The first line here registers our model so that it is known to Fabber by the name exp The second line 
// is a Factory method used so that Fabber can create a new instance of our model when its name appears on the command line:
FactoryRegistration<FwdModelFactory, PetFwdModel> PetFwdModel::registration("pet");

FwdModel *PetFwdModel::NewInstance()
{
    return new PetFwdModel();
}

string PetFwdModel::GetDescription() const
{
    return "Example model which uses a activity function from parametric parameters";
}

string PetFwdModel::ModelVersion() const
{
    return "1.0, add support to bowsher prior";
}

// We’ve given our model a version number, if we update it at some later stage we should change the number returned 
// so anybody using the model will know it has changed and what version they have. There’s also a brief description
//  which fabber will return when the user requests help on the model:
// Each option is listed in the OPTIONS array which ends with an empty option (important!)
//
// An option is described by:
// It’s name which generally should not include underscores (hyphen is OK as in this case). The name translates into a 
// command line option e.g. --num-exps.
// An option type. Possibilities are:
// OPT_BOOL for a Yes/No boolean
// OPT_FLOAT for a decimal number
// OPT_INT for a whole number
// OPT_STR for text
// OPT_MATRIX for a small matrix (specified by giving the filename of a text file which contains the matrix data in
//  tab-separated form)
// OPT_IMAGE for a 3D image specified as a Nifti file
// OPT_TIMESERIES for a 4D image specified as a Nifti file
// OPT_FILE for a generic filename
// A brief description of the option. This will be displayed when --help is requested for the model
// OPT_NONREQ if the option is not mandatory (does not need to be specified) or OPT_REQ if the option must be provided 
// by the user.
// An indication of the default value. This value is not actually used to initialize anything but is shown in --help to
//  explain to the user what the default is if the option is not given. So it can contain any text (e.g. "0.7 for PASL, 
//  1.3 for pCASL". You should not specify a default for a mandatory option (OPT_REQ)
static OptionSpec OPTIONS[] = {
    // { "use-offset", OPT_BOOL, "If True, allow an additional constant offset parameter", OPT_NONREQ, "false" },
    { "plasma_t", OPT_MATRIX, "AIF time stamps", OPT_REQ, "false" },
    { "plasma_c", OPT_MATRIX, "AIF value samples in ml-1", OPT_REQ, "false" },
    { "nsample", OPT_INT, "number of frames", OPT_REQ, "false" },
    { "pmean1", OPT_INT, "prior means1", OPT_NONREQ, "0.0" },
    { "pmean2", OPT_INT, "prior means2", OPT_NONREQ, "0.0" },
    { "pmean3", OPT_INT, "prior means3", OPT_NONREQ, "0.0" },
    { "pmean4", OPT_INT, "prior means4", OPT_NONREQ, "0.0" },
    { "pprecision1", OPT_INT, "prior precisions1", OPT_NONREQ, "0.0" },
    { "pprecision2", OPT_INT, "prior precisions2", OPT_NONREQ, "0.0" },
    { "pprecision3", OPT_INT, "prior precisions3", OPT_NONREQ, "0.0" },
    { "pprecision4", OPT_INT, "prior precisions4", OPT_NONREQ, "0.0" },
    { "usepriorimg", OPT_BOOL, "use prior image or not",  OPT_NONREQ, "0.0"},
    { "bowsherthreshold", OPT_IMAGE, "bowsher threshold from idl", OPT_NONREQ, "0.0" },  
    { "bowsherlabel", OPT_IMAGE, "bowsher threshold from idl", OPT_NONREQ, "0.0" },  
    { "maskref", OPT_IMAGE, "the mask to indicate the voxels exclude from fitting", OPT_NONREQ, "0.0" },  
    // { "pmeanimg1", OPT_IMAGE, "prior mean image1", OPT_NONREQ, "0.0" },
    // { "pmeanimg2", OPT_IMAGE, "prior mean image2", OPT_NONREQ, "0.0" },
    // { "pmeanimg1", OPT_IMAGE, "prior mean image1", OPT_NONREQ, "0.0" },
    // { "pmeanimg2", OPT_IMAGE, "prior mean image2", OPT_NONREQ, "0.0" },
    { "" }
};

void PetFwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

// Options specified by the user are captured in the FabberRunData object which we use to set the variables in our model
//  class in the Initialize method. Initialize is called before the model will be used. Its purpose is to allow the model
//   to set up any internal variables based on the user-supplied options. Here we capture the time resolution option and
//    the number of exponentials - note that the latter has a default value:
void PetFwdModel::Initialize(FabberRunData &rundata)
{
    // int num_params = NumParams();
    usepriorimg = rundata.ReadBool("usepriorimg");
    // string designFile = rundata.GetString("pmeanimg1");
    // pmeanimg1 = rundata.GetVoxelData(designFile).AsRow();
    // designFile = rundata.GetString("pmeanimg2");
    // pmeanimg2 = rundata.GetVoxelData(designFile).AsRow();
    // pmeanimg1 = rundata.LoadVoxelData("pmeanimg1");
    // pmeanimg2 = rundata.LoadVoxelData("pmeanimg2");
    pmean1 = rundata.GetDoubleDefault("pmean1", 0.1);
    pmean2 = rundata.GetDoubleDefault("pmean2", 0.1);
    pmean3 = rundata.GetDoubleDefault("pmean3", 0.1);
    pmean4 = rundata.GetDoubleDefault("pmean4", 0.1);        
    pprecision1 = rundata.GetDoubleDefault("pprecision1", 1e-8); 
    pprecision2 = rundata.GetDoubleDefault("pprecision2", 1e-8); 
	pprecision3 = rundata.GetDoubleDefault("pprecision3", 1e-8); 
    pprecision4 = rundata.GetDoubleDefault("pprecision4", 1e-8); 
    nsample = rundata.GetDoubleDefault("nsample", 2);    //default is 2
    plasma_t = (double *)malloc(nsample*sizeof(double));
    plasma_c = (double *)malloc(nsample*sizeof(double));
        

    m_include_offset = rundata.GetBool("use-offset");
    m_usesrtm = rundata.GetBool("use-srtm");
    m_usertcm = rundata.GetBool("use-rtcm");
    string designFile = rundata.GetString("plasma_t");
    LOG << "LinearFwdModel::Reading design file: " << designFile << endl;
    FILE *myFile;
    myFile = fopen(designFile.c_str(), "r");
    LOG << designFile.c_str()<< endl;
     for (int i = 0; i < nsample; i++)
    {
        fscanf(myFile, "%lf", &plasma_t[i]);
        LOG << plasma_t[i] << endl;
    }
    // plasma_t = fabber::read_matrix_file(designFile);

    designFile = rundata.GetString("plasma_c");
    LOG << designFile.c_str()<< endl;
    LOG << "LinearFwdModel::Reading design file: " << designFile << endl;
    myFile = fopen(designFile.c_str(), "r");
     for (int i = 0; i < nsample; i++)
    {
        fscanf(myFile, "%lf", &plasma_c[i]);
        LOG << plasma_c[i] << endl;
    }
    // plasma_c = fabber::read_matrix_file(designFile);
     
    free(plasma_t);
	free(plasma_c);
}

int PetFwdModel::NumParams() const
{
    if (m_include_offset)
        return 5; //5;
    else
        return 4; //4;
}

void PetFwdModel::NameParams(vector<string> &names) const
{
    names.clear();
    names.push_back("K1");
    names.push_back("k2");
    names.push_back("k3");
    names.push_back("k4");
    if (m_include_offset)
        names.push_back("spillover");
}

// Next we need to specify what parameters our model includes.
// void ExpFwdModel::GetParameterDefaults(std::vector<Parameter> &params) const
// {
//     params.clear();

//     int p=0;
//     for (int i=0; i<m_num; i++) {
//         params.push_back(Parameter(p++, "amp" + stringify(i+1), DistParams(1, 100), DistParams(1, 100), PRIOR_NORMAL, TRANSFORM_LOG()));
//         params.push_back(Parameter(p++, "r" + stringify(i+1), DistParams(1, 100), DistParams(1, 100), PRIOR_NORMAL, TRANSFORM_LOG()));
//     }
// }
// GetParameterDefaults is quite important. It declares the parameters our model takes, and their prior and initial posterior distributions. 
// It is always called after Initialize so you can use whatever options you have set up to decide what parameters to include.

// Priors are central to Bayesian inference, and describe the extent of our belief about a parameter’s value before we have seen any data.
// For example if a parameter represents the T_1 value of grey matter in the brain there is a well known range of plausible values. 
// By declaring a suitable prior we ensure that probabilities are calculated correctly and unlikely values of the parameter are avoided 
// unless the data very strongly supports this.
// In our case we have no real prior information, so we are using an uninformative prior. This has a large variance so the model has a 
// lot of freedom in fitting the parameters and will try to get as close to matching the data as it can. This is reflected in the high 
// variance we are using (1e6). For the mean values, a and b are multiplicative so it makes sense to give them defaults of 1 wherease c
//  and d are additive so prior means of 0 seems more appropriate.
// The second DistParams instance represents the initial posterior. This is the starting point for the optimisation as it tries to find
//  the best values for each parameter. Usually this does not matter too much and can often be set to be identical to the prior.
// Sometimes, however, it may be helpful to give the initial posterior a more restrictive (lower) variance to avoid numerical instability.
// It is also possible to adjust the initial posterior on a per-voxel basis using the actual voxel data. We will not do that here, but it
//  can be useful when fitting, for example, a constant offset, where we can tell the optimisation to start with a value that is the mean
//   of the data. This may help avoid instability and local minima.
// In general it is against the spirit of the Bayesian approach to modify the priors on the basis of the data, and no means are provided to
//  do this. It is possible for the user to modify the priors on a global basis but this is not encouraged and in general a model should try to provide good priors that will not need modification.
void PetFwdModel::HardcodedInitialDists(MVNDist &prior, MVNDist &posterior) const
{
    int num_params = NumParams();
    // Check we have been given a distribution of the right number of parameters
    assert(prior.means.Nrows() == num_params);
    // if (usepriorimg == 1) {
    //     // GetParameterDefaults(std::vector<Parameter> &params)?
	// 	LOG<<"pan"<<voxel<<endl;
    //     prior.means(1) = pmeanimg1(voxel+1);              //0.6/2;
    //     prior.means(2) = pmeanimg2(voxel+1);
	// 		LOG<<"pan"<<endl;
    // } else {
        prior.means(1) = pmean1;              //0.6/2;
        prior.means(2) = pmean2;              //1.4;
        prior.means(3) = pmean3;
        prior.means(4) = pmean4; 
    // }

    SymmetricMatrix precisions = IdentityMatrix(num_params);   //* 1e1;     //1e-2;
    precisions(1,1) *= pprecision1; 
    precisions(2,2) *= pprecision2;     
	precisions(3,3) *= pprecision3; 
    precisions(4,4) *= pprecision4; 
    prior.SetPrecisions(precisions);   // 1e-12
    posterior = prior;
}

// We now go back to our model code where we finally reach the point where we write the code to calculate our model:
// We are given a list of parameter values (params) and need to produce a time series of predicted data values (result).
//  We do this by looping over the parameters and adding the result of each exponential to the output result.
// The additional argument key is not required in this case. It is used to allow a model to evaluate ‘alternative’ outputs 
// such as an interim residual or AIF curve
void PetFwdModel::Evaluate(const ColumnVector &params, ColumnVector &result) const
{

    int ret;
    double *results;
    // Check we have been given the right number of parameters
    assert(params.Nrows() == NumParams());
    result.ReSize(data.Nrows());

    results = (double *)malloc(nsample*sizeof(double));
    if (m_usesrtm) { ret=simSRTM(plasma_t, plasma_c, nsample, params(1),params(2),params(3), results); }
    else if (m_usertcm) { ret=simRTCM(plasma_t, plasma_c, nsample, params(1),params(2),params(3),params(4), results, NULL, NULL);  }
    else { ret=simC2(plasma_t, plasma_c, nsample, params(1),params(2),params(3),params(4), results, NULL, NULL); }
	
    std::vector<double> result1(results, results + nsample);
    for (int i = 1; i <= data.Nrows(); i++)
    {
        // result(i) = result1[i];
        result(i) = *(results+i-1);    // start from first image
    }
	
	free(results);


}

