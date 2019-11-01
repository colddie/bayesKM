/*  fwdmodel.cc - base class for generic forward models

 Adrian Groves and Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford  */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk

    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK


    LICENCE

    FMRIB Software Library, Release 6.0 (c) 2018, The University of
    Oxford (the "Software")

    The Software remains the property of the Oxford University Innovation
    ("the University").

    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.

    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.

    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.

    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Oxford
    University Innovation ("OUI"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    fsl@innovation.ox.ac.uk quoting Reference Project 9564, FSL.*/

#include "fwdmodel.h"

#include "easylog.h"
#include "priors.h"
#include "rundata.h"
#include "transforms.h"

#include <newmatio.h>

#include <memory>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

typedef int (*GetNumModelsFptr)(void);
typedef const char *(*GetModelNameFptr)(int);
typedef NewInstanceFptr (*GetNewInstanceFptrFptr)(const char *);

#ifdef _WIN32
// This stops Windows defining a load of macros which clash with FSL
#define WIN32_LEAN_AND_MEAN
#include "windows.h"
#define GETSYMBOL GetProcAddress
#define GETERROR GetLastErrorAsString

string GetLastErrorAsString()
{
    // Get the error message, if any.
    DWORD errorMessageID = ::GetLastError();
    if (errorMessageID == 0)
        return std::string(); // No error message has been recorded

    LPSTR messageBuffer = NULL;
    size_t size = FormatMessageA(
        FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
        NULL, errorMessageID, MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), (LPSTR)&messageBuffer, 0,
        NULL);

    std::string message(messageBuffer, size);

    // Free the buffer.
    LocalFree(messageBuffer);

    return message;
}
#else
// POSIX-style methods for shared libraries
#include <dlfcn.h>
#define GETSYMBOL dlsym
#define GETERROR dlerror
#endif

void FwdModel::LoadFromDynamicLibrary(const std::string &filename, EasyLog *log)
{
    FwdModelFactory *factory = FwdModelFactory::GetInstance();
    GetNumModelsFptr get_num_models;
    GetModelNameFptr get_model_name;
    GetNewInstanceFptrFptr get_new_instance_fptr;
    if (log)
        log->LogStream() << "Loading dynamic models from " << filename << endl;

#ifdef _WIN32
    HINSTANCE libptr = LoadLibrary(filename.c_str());
#else
    void *libptr = dlopen(filename.c_str(), RTLD_NOW);
#endif
    if (!libptr)
    {
        throw InvalidOptionValue(
            "loadmodels", filename, string("Failed to open library ") + GETERROR());
    }

    get_num_models = (GetNumModelsFptr)GETSYMBOL(libptr, "get_num_models");
    if (!get_num_models)
    {
        throw InvalidOptionValue("loadmodels", filename,
            string("Failed to resolve symbol 'get_num_models' ") + GETERROR());
    }

    get_model_name = (GetModelNameFptr)GETSYMBOL(libptr, "get_model_name");
    if (!get_model_name)
    {
        throw InvalidOptionValue("loadmodels", filename,
            string("Failed to resolve symbol 'get_model_name' ") + GETERROR());
    }

    get_new_instance_fptr = (GetNewInstanceFptrFptr)GETSYMBOL(libptr, "get_new_instance_func");
    if (!get_new_instance_fptr)
    {
        throw InvalidOptionValue("loadmodels", filename,
            string("Failed to resolve symbol 'get_new_instance_func' ") + GETERROR());
    }

    int num_models = get_num_models();
    if (log)
        log->LogStream() << "Loading " << num_models << " models" << endl;
    for (int i = 0; i < num_models; i++)
    {
        const char *model_name = get_model_name(i);
        if (!model_name)
        {
            throw InvalidOptionValue("loadmodels", filename,
                "Dynamic library failed to return model name for index " + stringify(i));
        }
        else
        {
            if (log)
                log->LogStream() << "Loading model " << model_name << endl;
            NewInstanceFptr new_instance_fptr = get_new_instance_fptr(model_name);
            if (!new_instance_fptr)
            {
                throw InvalidOptionValue("loadmodels", filename,
                    string("Dynamic library failed to return new instance function for model")
                        + model_name);
            }
            factory->Add(model_name, new_instance_fptr);
        }
    }
}

std::vector<std::string> FwdModel::GetKnown()
{
    FwdModelFactory *factory = FwdModelFactory::GetInstance();
    return factory->GetNames();
}

FwdModel *FwdModel::NewFromName(const string &name)
{
    FwdModelFactory *factory = FwdModelFactory::GetInstance();
    FwdModel *model = factory->Create(name);
    if (model == NULL)
    {
        throw InvalidOptionValue("model", name, "Unrecognized forward model");
    }
    return model;
}

void FwdModel::Initialize(FabberRunData &args)
{
    m_log = args.GetLogger();
}

void FwdModel::UsageFromName(const string &name, std::ostream &stream)
{
    std::auto_ptr<FwdModel> model(NewFromName(name));
    stream << name << ": " << model->ModelVersion() << endl << endl;
    stream << model->GetDescription() << endl << endl;
    stream << "Options: " << endl << endl;
    vector<OptionSpec> options;
    model->GetOptions(options);
    if (options.size() > 0)
    {
        for (vector<OptionSpec>::iterator iter = options.begin(); iter != options.end(); ++iter)
        {
            stream << *iter;
        }
    }
    else
    {
        model->Usage(stream);
    }
    vector<string> outputs;
    model->GetOutputs(outputs);
    if (outputs.size() > 0)
    {
        stream << endl << "Additional outputs: " << endl << endl;
        for (vector<string>::iterator iter = outputs.begin(); iter != outputs.end(); ++iter)
        {
            if (*iter != "")
                stream << "  " << *iter << endl;
        }
    }
}

string FwdModel::GetDescription() const
{
    return "No description available";
}
string FwdModel::ModelVersion() const
{
    return "No version info available.";
}
void FwdModel::Usage(std::ostream &stream) const
{
    stream << "No usage information available" << endl;
}

void FwdModel::PassData(unsigned int voxel_idx, const NEWMAT::ColumnVector &voxdata, const NEWMAT::ColumnVector &voxcoords,
    const NEWMAT::ColumnVector &voxsuppdata)
{
    voxel = voxel_idx;
    data = voxdata;
    suppdata = voxsuppdata;
    coords = voxcoords;
    coord_x = coords(1);
    coord_y = coords(2);
    coord_z = coords(3);
}

void FwdModel::GetParameters(FabberRunData &rundata, vector<Parameter> &params)
{

   GetParameterDefaults(params);
    m_params.clear();

    for (vector<Parameter>::iterator p = params.begin(); p < params.end(); ++p)
    {
        // Complexity below is due to there being two ways of specifying
        // priors. One is using the param-spatial-priors option which is
        // a sequence of chars in model parameter order, one for each
        // parameter. A + character means 'use the previous value for all
        // remaining parameters'. An 'I' means an image prior and
        // the filename is specified separately using an image-prior<n> option
        string types = Prior::ExpandPriorTypesString(
            rundata.GetStringDefault("param-spatial-priors", ""), params.size());
        assert(types.size() == params.size());
        if (types[p->idx] != PRIOR_DEFAULT)
        {
            p->prior_type = types[p->idx];
        }

        // Record the data key (filename) for an image prior. Note that the index is
        // conceptually different from the PSP_byname_image method use below - here
        // it is the parameter index in the model's list (starting at 1), below it depends on
        // the order in which the names are given in the options.
        p->options["image"] = "image-prior" + stringify(p->idx + 1);

        // Determine if we have any PSP_byname options for this parameter. These override the
        // options above
        int psp_idx = 1;
        while (true)
        {
            string name = rundata.GetStringDefault("PSP_byname" + stringify(psp_idx), "stop!");
            if (name == "stop!")
                break;
            else if (name == p->name)
            {
                string psp_idx_str = stringify(psp_idx);
                string transform_code
                    = rundata.GetStringDefault("PSP_byname" + psp_idx_str + "_transform", "");
                if (transform_code != "")
                    p->transform = GetTransform(transform_code);

                char prior_type = convertTo<char>(rundata.GetStringDefault(
                    "PSP_byname" + psp_idx_str + "_type", stringify(p->prior_type)));
                if (prior_type != PRIOR_DEFAULT)
                    p->prior_type = prior_type;

				
                double mean = rundata.GetDoubleDefault(
                    "PSP_byname" + psp_idx_str + "_mean", p->prior.mean());
                double prec = rundata.GetDoubleDefault(
                    "PSP_byname" + psp_idx_str + "_prec", p->prior.prec());
                p->prior = DistParams(mean, 1 / prec);
                p->options["image"] = "PSP_byname" + psp_idx_str + "_image";
                p->options["pimage"] = "PSP_byname" + psp_idx_str + "_pimage";            }
            psp_idx++;
        }

        if (p->prior.prec() > 1e12) {
            WARN_ONCE("Specified precision " + stringify(p->prior.prec()) + " is very high - this can trigger numerical instability. Using 1e12 instead");
            p->prior = DistParams(p->prior.mean(), 1e-12);
        }

        // FIXME do this here, or let the priors do it?
        //
        // Need to transform mean/precision as specified in the model into Fabber-space
        // Note that posterior is transformed in GetInitialPosterior
        p->prior = p->transform->ToFabber(p->prior);

        // Keep our own list of parameters
        m_params.push_back(*p);
    }
}

void FwdModel::GetInitialPosterior(MVNDist &posterior) const
{
    posterior.SetSize(m_params.size());

    // Set model defaults
    NEWMAT::SymmetricMatrix cov = posterior.GetCovariance();
    for (size_t p = 0; p < m_params.size(); p++)
    {
        posterior.means(p + 1) = m_params[p].post.mean();
        cov(p + 1, p + 1) = m_params[p].post.var();
    }
    posterior.SetCovariance(cov);

    // Do voxelwise initialization
    InitVoxelPosterior(posterior);

    // Finally, apply transforms
    ToFabber(posterior);
}

void FwdModel::ToFabber(MVNDist &mvn) const
{
    NEWMAT::SymmetricMatrix cov = mvn.GetCovariance();
    for (size_t p = 0; p < m_params.size(); p++)
    {
        mvn.means(p + 1) = m_params[p].transform->ToFabber(mvn.means(p + 1));
        cov(p + 1, p + 1) = m_params[p].transform->ToFabberVar(cov(p + 1, p + 1));
    }
    mvn.SetCovariance(cov);
}

void FwdModel::ToModel(MVNDist &mvn) const
{
    NEWMAT::SymmetricMatrix cov = mvn.GetCovariance();
    for (size_t p = 0; p < m_params.size(); p++)
    {
        DistParams dp(mvn.means(p + 1), cov(p + 1, p + 1));
        dp = m_params[p].transform->ToModel(dp);
        mvn.means(p + 1) = dp.mean();
        cov(p + 1, p + 1) = dp.var();
    }
    mvn.SetCovariance(cov);
}

void FwdModel::GetParameterDefaults(vector<Parameter> &params) const
{
    params.clear();
    vector<string> names;
    // Old method of naming parameters
    NameParams(names);

    // Old method of specifying default prior and posterior
    MVNDist priors(names.size()), posts(names.size());
	
    HardcodedInitialDists(priors, posts);

    for (unsigned int i = 0; i < names.size(); i++)
    {
        DistParams prior(priors.means(i + 1), priors.GetCovariance()(i + 1, i + 1));
        DistParams post(posts.means(i + 1), posts.GetCovariance()(i + 1, i + 1));
        Parameter p(i, names[i], prior, post, PRIOR_NORMAL, TRANSFORM_IDENTITY());

        //Old method of specifying ARD priors
        if (find(ardindices.begin(), ardindices.end(), i + 1) != ardindices.end())
        {
            p.prior_type = PRIOR_ARD;
        }
        params.push_back(p);
    }
}

void FwdModel::EvaluateFabber(
    const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result, const std::string &key) const
{
    assert((m_params.size() == 0) || (int(m_params.size()) == params.Nrows()));
    if (m_params.size() == 0)
    {
        EvaluateModel(params, result, key);
    }
    else
    {
        NEWMAT::ColumnVector tparams(params.Nrows());
        for (int i = 1; i <= params.Nrows(); i++)
        {
            tparams(i) = m_params[i - 1].transform->ToModel(params(i));
        }
        EvaluateModel(tparams, result, key);
    }
}

void FwdModel::DumpParameters(const NEWMAT::ColumnVector &params, const string &indent) const
{
    LOG << indent << "Parameters:" << endl;
    vector<string> names;
    NameParams(names);
    assert(int(names.size()) == params.Nrows());

    for (size_t i = 1; i <= names.size(); i++)
        LOG << indent << "  " << names[i - 1] << " = " << params(i) << endl;

    LOG << indent << "Total of " << names.size() << " parameters" << endl;
}
