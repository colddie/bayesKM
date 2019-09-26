#ifndef NO_NLLS

/* inference_nlls.cc - Non-Linear Least Squares class declarations

 Adrian Groves Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

 Copyright (C) 2007-2015 University of Oxford */

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

#include "inference_nlls.h"

#include "easylog.h"
#include "fwdmodel.h"
#include "priors.h"
#include "rundata.h"
#include "tools.h"
#include "version.h"

#include <newmat.h>

#include <string>
#include <vector>

using namespace std;
using namespace MISCMATHS;
using namespace NEWMAT;
using fabber::MaskRows;

static int NUM_OPTIONS = 1;
static OptionSpec OPTIONS[] = {
    { "vb-init", OPT_BOOL, "Whether NLLS is being run in isolation or as a pre-step for VB",
        OPT_NONREQ, "" },
    { "lm", OPT_BOOL, "Whether to use LM convergence (default is L)", OPT_NONREQ, "" },
};

void NLLSInferenceTechnique::GetOptions(vector<OptionSpec> &opts) const
{
    InferenceTechnique::GetOptions(opts);
    for (int i = 0; i < NUM_OPTIONS; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

std::string NLLSInferenceTechnique::GetDescription() const
{
    return "Non-linear least squares inference technique.";
}

string NLLSInferenceTechnique::GetVersion() const
{
    return fabber_version();
}
InferenceTechnique *NLLSInferenceTechnique::NewInstance()
{
    return new NLLSInferenceTechnique();
}
void NLLSInferenceTechnique::Initialize(FwdModel *fwd_model, FabberRunData &args)
{
    InferenceTechnique::Initialize(fwd_model, args);
    LOG << "NLLSInferenceTechnique::Initialising" << endl;

    // Determine whether NLLS is being run in isolation or as a pre-step for VB
    // This alters what we do if result is ill conditioned
    m_vbinit = args.GetBool("vb-init");

    // Initialize the model with MVN distributions for its parameters
    MVNDist *loadPosterior = new MVNDist(m_model->NumParams());
    MVNDist *junk = new MVNDist(m_model->NumParams());
    m_model->HardcodedInitialDists(*junk, *loadPosterior);

    // Option to load a 'posterior' which will allow the setting of intial parameter estimates for
    // NLLS
    string filePosterior = args.GetStringDefault("fwd-inital-posterior", "modeldefault");
    if (filePosterior != "modeldefault")
    {
        LOG << "NLLSInferenceTechnique::File posterior" << endl;
        loadPosterior->LoadFromMatrix(filePosterior);
    }

    initialFwdPosterior = loadPosterior;
    loadPosterior = NULL;

    // Determine whether we use L (default) or LM convergence
    m_lm = args.GetBool("lm");
    LOG << "NLLSInferenceTechnique::Done initialising" << endl;
}

void NLLSInferenceTechnique::DoCalculations(FabberRunData &allData)
{
    // Get basic voxel data
    const Matrix &data = allData.GetMainVoxelData();
    const Matrix &coords = allData.GetVoxelCoords();
    unsigned int Nvoxels = data.Ncols();

    // pass in some (dummy) data/coords here just in case the model relies upon it
    // use the first voxel values as our dummies
    if (Nvoxels > 0)
    {
        m_model->PassData(1, data.Column(1), coords.Column(1));
    }

    // Check how many samples in time series (ignoring any masked time points)
    int Nsamples = data.Nrows() - m_masked_tpoints.size();

    // Loop over voxels. The result for each voxel is
    // stored as a MVN distribution for its parameters
    // in resultMVNs.
    for (unsigned int voxel = 1; voxel <= Nvoxels; voxel++)
    {
        ColumnVector y = data.Column(voxel);
        ColumnVector vcoords = coords.Column(voxel);

        // Some models might want more information about the data
        m_model->PassData(voxel, y, vcoords);

        LinearizedFwdModel linear(m_model);

        // FIXME should be a single sensible way to get the
        // number of model parameters!
        int Nparams = initialFwdPosterior->GetSize();

        // FIXME how about a ctor for MVNDist which takes a size?
        MVNDist fwdPosterior;
        fwdPosterior.SetSize(Nparams);

        IdentityMatrix I(Nparams);

        // Create a cost function evaluator which will
        // measure the difference between the model
        // and the data
        NLLSCF costfn(y, m_model, m_masked_tpoints);

        // Set the convergence method
        // either Levenberg (L) or Levenberg-Marquardt (LM)
        NonlinParam nlinpar(Nparams, NL_LM);
        if (!m_lm)
        {
            nlinpar.SetGaussNewtonType(LM_L);
        }

        // set ics from 'posterior'
        ColumnVector nlinics = initialFwdPosterior->means;
        nlinpar.SetStartingEstimate(nlinics);
        nlinpar.LogPar(true);
        nlinpar.LogCF(true);

        try
        {
            // Run the nonlinear optimizer
            // output variable is unused - unsure if nonlin has any effect
            nonlin(nlinpar, costfn);

#if 0
			LOG << "NLLSInferenceTechnique::The solution is: " << nlinpar.Par() << endl;
			LOG << "NLLSInferenceTechnique::and this is the process " << endl;
			for (int i=0; i<nlinpar.CFHistory().size(); i++)
			{
				LOG << " cf: " << (nlinpar.CFHistory())[i] <<endl;
			}
			for (int i=0; i<nlinpar.ParHistory().size(); i++)
			{
				LOG << (nlinpar.ParHistory())[i] << ": :";
			}
#endif
            // Get the new parameters
            fwdPosterior.means = nlinpar.Par();

            // Recenter linearized model on new parameters
            linear.ReCentre(fwdPosterior.means);
            Matrix J = linear.Jacobian();
            MaskRows(J, m_masked_tpoints);

            // Calculate the NLLS precision
            // This is (J'*J)/mse
            // The covariance is the inverse
            SymmetricMatrix nllsprec;
            double sqerr = costfn.cf(fwdPosterior.means);
            double mse = sqerr / (Nsamples - Nparams);
            nllsprec << J.t() * J / mse;

            // Look for zero diagonal elements (implies parameter is not observable)
            // and set precision small, but non-zero - so that covariance can be calculated
            for (int i = 1; i <= nllsprec.Nrows(); i++)
            {
                if (nllsprec(i, i) < 1e-6)
                {
                    nllsprec(i, i) = 1e-6;
                }
            }
            fwdPosterior.SetPrecisions(nllsprec);
            fwdPosterior.GetCovariance();
        }
        catch (Exception &e)
        {
            LOG << "NLLSInferenceTechnique::NEWMAT Exception in this voxel:\n" << e.what() << endl;

            if (m_halt_bad_voxel)
                throw;

            LOG << "NLLSInferenceTechnique::Estimates in this voxel may be unreliable" << endl
                << "   (precision matrix will be set manually)" << endl
                << "   Going on to the next voxel" << endl;

            // output the results where we are
            fwdPosterior.means = nlinpar.Par();

            // recenter linearized model on new parameters
            linear.ReCentre(fwdPosterior.means);

            // precision matrix is probably singular so set manually
            fwdPosterior.SetPrecisions(I * 1e-12);
        }

        resultMVNs.push_back(new MVNDist(fwdPosterior));
        assert(resultMVNs.size() == voxel);
    }
}

NLLSCF::NLLSCF(
    const NEWMAT::ColumnVector &pdata, const FwdModel *pm, std::vector<int> masked_tpoints)
    : m_data(MaskRows(pdata, masked_tpoints))
    , m_model(pm)
    , m_linear(pm)
    , m_masked_tpoints(masked_tpoints)
{
}

double NLLSCF::cf(const ColumnVector &p) const
{
    // p = parameters
    // data_pred = data predicted by model
    ColumnVector data_pred;
    m_model->EvaluateFabber(p, data_pred, "");
    data_pred = MaskRows(data_pred, m_masked_tpoints);

    // m_data = actual data. Find sum of squares of differences
    // between this and the model data using a scalar product.
    double cfv = ((m_data - data_pred).t() * (m_data - data_pred)).AsScalar();
    return cfv;
}

ReturnMatrix NLLSCF::grad(const ColumnVector &p) const
{
    // Create an initial zero gradient vector
    ColumnVector gradv(p.Nrows());

    // Need to recenter the linearised model to the current parameter values
    m_linear.ReCentre(p);
    Matrix J = m_linear.Jacobian();
    J = MaskRows(J, m_masked_tpoints);

    // Evaluate the model given the parameters
    ColumnVector data_pred;
    m_model->EvaluateFabber(p, data_pred, "");
    data_pred = MaskRows(data_pred, m_masked_tpoints);

    gradv = -2 * J.t() * (m_data - data_pred);
    gradv.Release();
    return (gradv);
}

boost::shared_ptr<BFMatrix> NLLSCF::hess(
    const ColumnVector &p, boost::shared_ptr<BFMatrix> iptr) const
{
    boost::shared_ptr<BFMatrix> hessm;

    if (iptr && iptr->Nrows() == (unsigned)p.Nrows() && iptr->Ncols() == (unsigned)p.Nrows())
    {
        hessm = iptr;
    }
    else
    {
        hessm = boost::shared_ptr<BFMatrix>(new FullBFMatrix(p.Nrows(), p.Nrows()));
    }

    // need to recenter the linearised model to the current parameter values
    m_linear.ReCentre(p);
    Matrix J = m_linear.Jacobian();
    J = MaskRows(J, m_masked_tpoints);
    Matrix hesstemp = 2 * J.t() * J; // Make the G-N approximation to the hessian

    //(*hessm) = J.t()*J;

    for (int i = 1; i <= p.Nrows(); i++)
    {
        for (int j = 1; j <= p.Nrows(); j++)
            hessm->Set(i, j, hesstemp(i, j));
    }

    return (hessm);
}

#endif
