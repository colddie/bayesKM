/**
 * prior.cc
 *
 * Class for a parameter prior
 *
 * Copyright (C) 2007-2017 University of Oxford
 */

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

#include "priors.h"

#include "dist_mvn.h"
#include "rundata.h"
#include "tools.h"

#include <miscmaths/miscmaths.h>
#include <newmat.h>

#include <math.h>
#include <ostream>
#include <string>
#include <vector>

using namespace std;
using namespace NEWMAT;
using MISCMATHS::digamma;

std::ostream &operator<<(std::ostream &out, const Prior &prior)
{
    prior.DumpInfo(out);
    return out;
}

string Prior::ExpandPriorTypesString(string priors_str, unsigned int num_params)
{
    // Find out how many prior types are in the string, and what the + character
    // should be interpreted as
    unsigned int n_str_params = 0;
    char repeat_type = '-';
    bool plus_found = false;
    for (size_t i = 0; i < priors_str.size(); i++)
    {
        if (priors_str[i] != '+')
        {
            if (!plus_found)
                repeat_type = priors_str[i];
            n_str_params++;
        }
        else if (plus_found)
        {
            throw InvalidOptionValue(
                "param-spatial-priors", priors_str, "Only one + character allowed");
        }
        else
        {
            plus_found = true;
        }
    }

    if (n_str_params > num_params)
    {
        throw InvalidOptionValue("param-spatial-priors", priors_str, "Too many parameters");
    }
    else if (n_str_params < num_params)
    {
        // Expand '+' char, if present, to give correct number of parameters
        // If there is no +, append with '-', meaning 'model default'
        int deficit = num_params - n_str_params;
        size_t plus_pos = priors_str.find("+");
        if (plus_pos != std::string::npos)
        {
            priors_str.insert(plus_pos, deficit - 1, '+');
        }
        else
        {
            priors_str.insert(priors_str.end(), deficit, '-');
        }
    }
    else
    {
        // We already have enough types for all the parameters so erase any
        // pointless + char
        priors_str.erase(std::remove(priors_str.begin(), priors_str.end(), '+'), priors_str.end());
    }

    // Finally, replace all + chars with identified repeat type
    std::replace(priors_str.begin(), priors_str.end(), '+', repeat_type);
    assert(priors_str.size() == num_params);

    return priors_str;
}

DefaultPrior::DefaultPrior(const Parameter &p, FabberRunData &rundata)
    : m_param_name(p.name)
    , m_idx(p.idx)
    , m_type_code(p.prior_type)
    , m_params(p.prior)
{
}

void DefaultPrior::DumpInfo(std::ostream &out) const
{
    out << "DefaultPrior: Parameter " << m_idx << " '" << m_param_name << "'"
        << " mean: " << m_params.mean() << " precision: " << m_params.prec();
}

//////////////////////////////////
double DefaultPrior::SetImgPrior(MVNDist *prior, MVNDist *posterior, const RunContext &ctx, FabberRunData &rundata,Parameter &param)
{
    return 0;
}

double DefaultPrior::ApplyToMVN_(MVNDist *prior, const RunContext &ctx, FabberRunData &rundata,Parameter &param)
{
    return 0;
}
//////////////////////////////////



double DefaultPrior::ApplyToMVN(MVNDist *prior, const RunContext &ctx, FabberRunData &rundata,Parameter &param)
{
    prior->means(m_idx + 1) = m_params.mean();

    SymmetricMatrix prec = prior->GetPrecisions();
    prec(m_idx + 1, m_idx + 1) = m_params.prec();
    prior->SetPrecisions(prec);

    return 0;
}

ImagePrior::ImagePrior(const Parameter &p, FabberRunData &rundata)
    : DefaultPrior(p,rundata)
{
    m_log = rundata.GetLogger();
    m_filename = p.options.find("image")->second;
    m_image = rundata.GetVoxelData(m_filename).AsRow();
}

void ImagePrior::DumpInfo(std::ostream &out) const
{
    out << "ImagePrior: Parameter " << m_idx << " '" << m_param_name << "'"
        << " filename: " << m_filename << " precision: " << m_params.prec();
}

double ImagePrior::ApplyToMVN(MVNDist *prior, const RunContext &ctx, FabberRunData &rundata,Parameter &param)
{
    prior->means(m_idx + 1) = m_image(ctx.v);

    SymmetricMatrix prec = prior->GetPrecisions();
    prec(m_idx + 1, m_idx + 1) = m_params.prec();
    prior->SetPrecisions(prec);

    return 0;
}

void ARDPrior::DumpInfo(std::ostream &out) const
{
    out << "ARDPrior: Parameter " << m_idx << " '" << m_param_name << "'"
        << " initial mean: " << m_params.mean() << " initial precision: " << m_params.prec();
}

double ARDPrior::ApplyToMVN(MVNDist *prior, const RunContext &ctx, FabberRunData &rundata,Parameter &param)
{
    SymmetricMatrix cov = prior->GetCovariance();
    double post_mean = ctx.fwd_post[ctx.v - 1].means(m_idx + 1);
    double post_cov = ctx.fwd_post[ctx.v - 1].GetCovariance()(m_idx + 1, m_idx + 1);
    // (Chappel et al 2009 Eq D4)
    double new_cov = post_mean * post_mean + post_cov;

    if (ctx.it == 0)
    {
        // Special case for first iter
        // Initially set prior to model default. The other alternative is to set
        // it to be initially non-informative, however this can still be achieved
        // by specifying the model prior mean/precision by options.
        cov(m_idx + 1, m_idx + 1) = m_params.var();
        prior->means(m_idx + 1) = m_params.mean();
        // LOG << "first iter ARD: " << m_params.var() << ", " << m_params.mean() << endl;
    }
    else
    {
        // LOG << "post: " << post_cov << ", " << post_mean << endl;
        // Update covariance on subsequent iterations
        cov(m_idx + 1, m_idx + 1) = new_cov;
        // LOG << "subs iter ARD: " << new_cov << ", " << prior->means(m_idx+1) << endl;
    }
    prior->SetCovariance(cov);

    // Calculate the free energy contribution from ARD term
    // (Chappel et al 2009, end of Appendix D)
    double b = 2 / new_cov;
    return -1.5 * (log(b) + digamma(0.5)) - 0.5 - gammaln(0.5) - 0.5 * log(b);
}

SpatialPrior::SpatialPrior(const Parameter &p, FabberRunData &rundata)
    : DefaultPrior(p,rundata)
    , m_akmean(1e-8)
    , m_spatial_dims(3)
    , m_spatial_speed(-1)
{
    m_log = rundata.GetLogger();
    m_spatial_dims = rundata.GetIntDefault("spatial-dims", 3);
    if (m_spatial_dims < 0 || m_spatial_dims > 3)
    {
        throw InvalidOptionValue("spatial-dims", stringify(m_spatial_dims), "Must be 0, 1, 2 or 3");
    }
    else if (m_spatial_dims == 1)
    {
        WARN_ONCE("spatial-dims=1 is very weird... hope you're just testing!");
    }
    else if (m_spatial_dims == 2)
    {
        WARN_ONCE("spatial-dims=2 doesn't decompose into slices");
    }

    // FIXME check valid range
    m_spatial_speed = rundata.GetDoubleDefault("spatial-speed", -1);

    // FIXME still needed?
    m_update_first_iter = rundata.GetBool("update-spatial-prior-on-first-iteration");

    // extract data (and the coords) from rundata 
    // Rows are volumes
    // Columns are (time) series
    // num Rows is size of (time) series
    // num Cols is size of volumes
    m_origdata = &rundata.GetMainVoxelData();
    m_nvoxels = m_origdata->Ncols();

}


//////////////////////////////////
double SpatialPrior::SetImgPrior(MVNDist *prior, MVNDist *posterior, const RunContext &ctx, FabberRunData &rundata,Parameter &param)
{
    // sptial dependent prior
    // LOG <<"test243"<<param.options.find("pmeanimg")->second;
    string m_filename = param.options.find("image")->second;
    NEWMAT::RowVector m_image = rundata.GetVoxelData(m_filename).AsRow();
    prior->means(m_idx + 1) = m_image(ctx.v);
    posterior->means(m_idx+1) = m_image(ctx.v); 

}//////////////////////////////////



double SpatialPrior::CalculateAkmean(const RunContext &ctx)
{
    // The following calculates Tr[Sigmak*S'*S]
    double tmp1 = 0.0;
    double tmp2 = 0.0;
    for (int v = 1; v <= ctx.nvoxels; v++)
    {
        // Ignore voxels where numerical issues have occurred
        if (std::find(ctx.ignore_voxels.begin(), ctx.ignore_voxels.end(), v)
            != ctx.ignore_voxels.end()) continue;

        double sigmak = ctx.fwd_post.at(v - 1).GetCovariance()(m_idx + 1, m_idx + 1);
        int nn = ctx.neighbours.at(v - 1).size();
        if (m_type_code == PRIOR_SPATIAL_m) // useMRF)
            tmp1 += sigmak * m_spatial_dims * 2;
        else if (m_type_code == PRIOR_SPATIAL_M) // useMRF2)
            tmp1 += sigmak * (nn + 1e-8);
        else if (m_type_code == PRIOR_SPATIAL_p)
            tmp1 += sigmak * (4 * m_spatial_dims * m_spatial_dims + nn);
        else // P
            tmp1 += sigmak * (nn * nn + nn);

        double wk = ctx.fwd_post.at(v - 1).means(m_idx + 1);
        double Swk = 0.0;
        for (vector<int>::const_iterator v2It = ctx.neighbours[v - 1].begin();
             v2It != ctx.neighbours.at(v - 1).end(); ++v2It)
        {
            Swk += wk - ctx.fwd_post.at(*v2It - 1).means(m_idx + 1);
        }
        if (m_type_code == PRIOR_SPATIAL_p || m_type_code == PRIOR_SPATIAL_m)
            Swk += wk * (m_spatial_dims * 2 - ctx.neighbours.at(v - 1).size());

        if (m_type_code == PRIOR_SPATIAL_m || m_type_code == PRIOR_SPATIAL_M)
            tmp2 += Swk * wk;
        else
            tmp2 += Swk * Swk;
    }

    LOG << "SpatialPrior::UpdateAkmean " << m_idx << ": tmp1=" << tmp1 << ", tmp2=" << tmp2 << endl;

    double gk = 1 / (0.5 * tmp1 + 0.5 * tmp2 + 0.1); // prior q1 == 10 (1/q1 == 0.1)
    double akmean = gk * (ctx.nvoxels * 0.5 + 1.0);  // prior q2 == 1.0
    double akmeanMax = akmean * m_spatial_speed;
    if (akmean < 1e-50)
    {
        LOG << "SpatialPrior::UpdateAkmean " << m_idx << ": was " << akmean << endl;
        WARN_ONCE("SpatialPrior::UpdateAkmean akmean value was tiny!");
        akmean = 1e-50; // prevent crashes
    }

    if (akmeanMax < 0.5)
    {
        akmeanMax = 0.5; // totally the wrong scaling.. oh well
    }

    if (m_spatial_speed > 0 && akmean > akmeanMax)
    {
        LOG << "SpatialPrior::UpdateAkmean " << m_idx
            << ": Rate-limiting the increase on akmean: was " << akmean << ", now " << akmeanMax
            << endl;
        akmean = akmeanMax;
    }

    LOG << "SpatialPrior::UpdateAkmean " << m_idx << ": New akmean: " << akmean << endl;
    return akmean;
}

void SpatialPrior::DumpInfo(std::ostream &out) const
{
    out << "SpatialPrior: Parameter " << m_idx << " '" << m_param_name << "'"
        << " type " << m_type_code << " mean: " << m_params.mean()
        << " precision: " << m_params.prec();
}


double SpatialPrior::ApplyToMVN(MVNDist *prior, const RunContext &ctx, FabberRunData &rundata, Parameter &param)
{


//////////////////////////////////
    // use bowsher
    NEWMAT::RowVector m_bowsherthresh;
    NEWMAT::RowVector m_bowsherlabel;
    string bowsher_fname = rundata.GetStringDefault("bowsherlabel", "");
    bool m_have_bowsher = (bowsher_fname != "");
    if (m_have_bowsher)
    {
        LOG << "used bowsher prior..." << endl;
        // string m_filename0 = param.options.find("threshold")->second;
        bowsher_fname = rundata.GetStringDefault("bowsherthreshold", "");
        m_bowsherthresh = rundata.GetVoxelData(bowsher_fname).AsRow();
        	
        // m_filename0 = param.options.find("label")->second;
        bowsher_fname = rundata.GetStringDefault("bowsherlabel", "");
        m_bowsherlabel = rundata.GetVoxelData(bowsher_fname).AsRow();
    }
/////////////////////////////////



if (m_type_code == PRIOR_SPATIAL_m || m_type_code == PRIOR_SPATIAL_M || 
    m_type_code == PRIOR_SPATIAL_p || m_type_code == PRIOR_SPATIAL_P)   // self-define priors
{

    if (ctx.v == 1 && (ctx.it > 0 || m_update_first_iter))
    {
        m_akmean = CalculateAkmean(ctx);
    }
	
//////////////////////////////////
    double anat1 = 0.0;
    double anat2 = 0.0;
    double thresh = 0.0;
    double labelsWeight = 1.0;
    unsigned int nn_bowsher = 0;
    if (m_have_bowsher)
    {
        anat1 = m_bowsherlabel[ctx.v-1];
        thresh = m_bowsherthresh[ctx.v-1];
    }
//////////////////////////////////

    double weight8 = 0; // weighted +8
    double contrib8 = 0.0;
    for (vector<int>::const_iterator nidIt = ctx.neighbours[ctx.v - 1].begin();
         nidIt != ctx.neighbours[ctx.v - 1].end(); ++nidIt)
    {
        int nid = *nidIt;
        const MVNDist &neighbourPost = ctx.fwd_post[nid - 1];

//////////////////////////////////
        if (m_have_bowsher)
        {
            anat2 = m_bowsherlabel[nid-1];
            labelsWeight = (fabs(anat1 - anat2) <= thresh);
            LOG << anat1 << anat2 << thresh << " " << labelsWeight << " "<< endl;
        }
//////////////////////////////////
        contrib8 += 8 * neighbourPost.means(m_idx + 1) *labelsWeight     ;  
        if (labelsWeight) { weight8 += 8; nn_bowsher += 1; }                            

        LOG << "incrementing neighbours " << nid << " " << contrib8 << " " <<labelsWeight<< " " << nn_bowsher << endl;
    }

    double weight12 = 0; // weighted -1, may be duplicated
    double contrib12 = 0.0;
    for (vector<int>::const_iterator nidIt = ctx.neighbours2[ctx.v - 1].begin();
         nidIt != ctx.neighbours2[ctx.v - 1].end(); ++nidIt)
    {
        int nid = *nidIt;
        const MVNDist &neighbourPost = ctx.fwd_post[nid - 1];
        contrib12 += -neighbourPost.means(m_idx + 1);
        weight12 += -1;
    }

    int nn = ctx.neighbours[ctx.v - 1].size();
    //////////////////
    if (m_have_bowsher) { nn = nn_bowsher;}
    //////////////////

    if (m_type_code == PRIOR_SPATIAL_p)
    {
        assert(nn <= m_spatial_dims * 2);
        weight8 = 8 * 2 * m_spatial_dims;
        weight12 = -1 * (4 * m_spatial_dims * m_spatial_dims - nn);
    }

    double spatial_prec = 0;

    if (m_type_code == PRIOR_SPATIAL_P)
        spatial_prec = m_akmean * (nn * nn + nn);
    else if (m_type_code == PRIOR_SPATIAL_m)
        spatial_prec = m_akmean * m_spatial_dims * 2;
    else if (m_type_code == PRIOR_SPATIAL_M)
        spatial_prec = m_akmean * (nn + 1e-8);
    else if (m_type_code == PRIOR_SPATIAL_p)
        spatial_prec = m_akmean * (4 * m_spatial_dims * m_spatial_dims + nn);
    else
        assert(false);

    // Set the prior precision for this parameter
    SymmetricMatrix precs = prior->GetPrecisions();
    if (m_type_code == PRIOR_SPATIAL_p || m_type_code == PRIOR_SPATIAL_m)
    {
        //	Penny-style DirichletBC priors -- ignoring initialFwdPrior completely!
        precs(m_idx + 1, m_idx + 1) = spatial_prec;
    }
    else
    {
        precs(m_idx + 1, m_idx + 1) = m_params.prec() + spatial_prec;
    }
    prior->SetPrecisions(precs);

    // Set the prior mean for this parameter
    // Note that we multiply by reciprocals rather than dividing. This is
    // to maximise numerical compatibility with NEWMAT which presumably
    // does it as an optimization when dividing a whole matrix by a constant
    double mTmp;
    if (m_type_code == PRIOR_SPATIAL_m)
    {
        // Dirichlet BCs on MRF
        double rec = 1 / (8 * m_spatial_dims * 2);
        mTmp = contrib8 * rec;
    }
    else if (m_type_code == PRIOR_SPATIAL_M)
    {
        double rec = 1 / (8 * (double(nn) + 1e-8));
        mTmp = contrib8 * rec;
    }
    else if (weight8 != 0)
    {
        double rec = 1 / (weight8 + weight12);
        mTmp = (contrib8 + contrib12) * rec;
    }
    else
        mTmp = 0;

    LOG << "SpatialPrior:: at voxel" << ctx.v  << ", " << 
    ", " << prior->GetCovariance()(m_idx+1, m_idx+1) << ", " << spatial_prec
    << ", " << contrib8 << ", " << mTmp << " : "  << endl;

    if (m_type_code == PRIOR_SPATIAL_m || m_type_code == PRIOR_SPATIAL_M)
        prior->means(m_idx + 1) = prior->GetCovariance()(m_idx + 1, m_idx + 1) * spatial_prec *mTmp  ; // = mTmp for p or m
    else
    {
        // equivalent, when non-spatial priors are very weak: m_fwd_prior[v-1].means = mTmp;
        prior->means(m_idx + 1) = prior->GetCovariance()(m_idx + 1, m_idx + 1)
            * (spatial_prec * mTmp + m_params.prec() * m_params.mean());
    }
			
}
else   // self-define priors
{	
		
if (ctx.v == 1 && (ctx.it > 0 || m_update_first_iter))
    {
        m_akmean = CalculateAkmean(ctx) ;     ////////////////////
    }

	// get the akmean
	double spatial_prec = 0;

    if (m_type_code == PRIOR_SPATIAL_n || m_type_code == PRIOR_SPATIAL_k)                 //;;;;;;;;;;;;;;;
        spatial_prec = m_akmean * (m_spatial_dims + 1e-8) ;      
    else 
        spatial_prec = m_akmean * m_spatial_dims * 2;


    // Set the prior precision for this parameter
    SymmetricMatrix precs = prior->GetPrecisions();
    if (m_type_code == PRIOR_SPATIAL_n || m_type_code == PRIOR_SPATIAL_k)     //;;;;;;;;;;;;;;;
    {
        //	Penny-style DirichletBC priors -- ignoring initialFwdPrior completely!
        precs(m_idx + 1, m_idx + 1) = m_params.prec() + spatial_prec;
    }
    else
    {
        precs(m_idx + 1, m_idx + 1) = spatial_prec;
    }
    prior->SetPrecisions(precs);
	

    // calculate neighbor weigths
    double tmp = 0.0;
    double mTmp = 0.0;
    double vnm = 0.0;
    double vm = 0.0;    //ctx.fwd_post[ctx.v - 1].means(m_idx + 1);         ///
    double weight = 0.0;
    double Sweight = 0.0;

    const int v = ctx.v;
//////////////////////////  
    double SD = rundata.GetDoubleDefault("sd", 1.0);     // 1.0/10; 
    // int size = ctx.neighboursnn[ctx.v - 1].size();   
    std::vector<double> vm1;
    std::vector<double> vnm1;
    // int search_w = 3;
    // int similarity_w = 3;
    // double vm1 = 0.0;
    // double vnm1 = 0.0;
//////////////////////////
    if (m_type_code == PRIOR_SPATIAL_n)
    { 
    LOG << "voxel " << ctx.v << " has neighboursn " <<  ctx.neighboursn[ctx.v - 1] << endl;
    LOG << "voxel " << ctx.v << " has neighboursnn " <<  ctx.neighboursnn[ctx.v - 1] << endl;
    LOG << "voxel " << ctx.v << " has neighbours1 " <<  ctx.neighbours[ctx.v - 1] << endl;

        for (vector<int>::const_iterator nidIt = ctx.neighboursnn[ctx.v - 1].begin();     // similarity window at current voxel
            nidIt != ctx.neighboursnn.at(ctx.v - 1).end(); ++nidIt) 
        {    
            int nid = *nidIt;  
            vm1.push_back(ctx.fwd_post[nid-1].means(m_idx + 1));
            // vm += ctx.fwd_post[nid- 1].means(m_idx + 1);      
        }

        for (vector<int>::const_iterator nidIt1 = ctx.neighboursn[ctx.v - 1].begin();    // search window at current voxel
            nidIt1 != ctx.neighboursn.at(ctx.v - 1).end(); ++nidIt1)
        {     
            int nid1 = *nidIt1;    
            for (vector<int>::const_iterator nidIt2 = ctx.neighboursnn[nid1 - 1].begin();    // similarity window at current voxel   
            nidIt2 != ctx.neighboursnn.at(nid1-1).end(); ++nidIt2) 
            {   
                int nid2 = *nidIt2; 
                vnm1.push_back(ctx.fwd_post[nid2-1].means(m_idx + 1));
                // const MVNDist &neighbourPost = ctx.fwd_post[nid2 - 1];
                // vnm = neighbourPost.means(m_idx + 1);
                // vnm += ctx.fwd_post[nid2 - 1].means(m_idx + 1);  
            }

            for (size_t ind=0; ind<vm1.size(); ind++) 
            {   
               tmp = vm1[ind] - vnm1[ind];
               weight += exp(-tmp * tmp /SD/SD);
            }
            Sweight += weight;    
            weight = 0.0;
            // tmp = vm - vnm;
            // Sweight += exp(-tmp * tmp /SD/SD);
            // LOG << "check here " << nid1 << " " << vnm << " ";
        }

        for (vector<int>::const_iterator nidIt1 = ctx.neighboursn[ctx.v - 1].begin();      // search window at current voxel
            nidIt1 != ctx.neighboursn.at(ctx.v - 1).end(); ++nidIt1)
        {     
            int nid1 = *nidIt1;     
            // for (vector<int>::const_iterator nidIt2 = ctx.neighboursnn[nid1 - 1].begin();    // similarity window at current voxel
            // nidIt2 != ctx.neighboursnn.at(nid1 - 1).end(); ++nidIt2) 
            // {   
            //     int nid2 = *nidIt2; 
            //     // vnm1[nid2] = ctx.fwd_post[nid2-1].means(m_idx + 1);
            //     // const MVNDist &neighbourPost = ctx.fwd_post[nid2 - 1];
            //     // vnm = neighbourPost.means(m_idx + 1);
            //     // vnm += ctx.fwd_post[nid2 - 1].means(m_idx + 1);  
            // }

            for (size_t ind=0; ind<vm1.size(); ind++) 
            {   
               tmp = vm1[ind] - vnm1[ind];
               weight += exp(-tmp * tmp /SD/SD);
            }
            mTmp +=  weight / Sweight * ctx.fwd_post[nid1 - 1].means(m_idx + 1) ;               /////////
            // LOG << weight << " " << endl;
            weight = 0.0;

            // tmp = vm - vnm;
            // mTmp += ctx.fwd_post[nid1 - 1].means(m_idx + 1) * exp(-tmp *tmp /SD/SD) / Sweight;
        }

    LOG << "SpatialPrior:: " << prior->GetCovariance()(m_idx+1, m_idx+1) << ", " << spatial_prec
     << ", " << vm << ", " << vnm << ", " << Sweight << ", " << mTmp << " : " << SD << endl;
    } 
    else if (m_type_code == PRIOR_SPATIAL_k)
    {
        ColumnVector vdata = m_origdata->Column(ctx.v);
        ColumnVector vndata;
        for (vector<int>::const_iterator nidIt = ctx.neighbours[v - 1].begin();
            nidIt != ctx.neighbours.at(v - 1).end(); ++nidIt)
        {
            vndata = m_origdata->Column(*nidIt);
            tmp = vdata.Sum() - vndata.Sum();         //last summation
            Sweight += exp(-tmp * tmp);
        }
        for (vector<int>::const_iterator nidIt = ctx.neighbours[v - 1].begin();
                nidIt != ctx.neighbours.at(v - 1).end(); ++nidIt)
        {
            vndata = m_origdata->Column(*nidIt);
            tmp = vdata.Sum() - vndata.Sum();
            int nid = *nidIt;     
            const MVNDist &neighbourPost = ctx.fwd_post[nid - 1];
            vnm = neighbourPost.means(m_idx + 1);
            mTmp += vnm * exp(-tmp *tmp /SD/SD) / Sweight;
        }
    LOG << "SpatialPrior:: " << prior->GetCovariance()(m_idx+1, m_idx+1) << ", " << spatial_prec
     << ", " << Sweight << ", " << mTmp << " : " << endl;
    }
        

    if (m_type_code == PRIOR_SPATIAL_n || m_type_code == PRIOR_SPATIAL_k) // non-local means prior
        prior->means(m_idx + 1) = prior->GetCovariance()(m_idx + 1, m_idx + 1) * spatial_prec * mTmp; 
    else  // delete later?
    {
        prior->means(m_idx + 1) = prior->GetCovariance()(m_idx + 1, m_idx + 1)
        * (spatial_prec * mTmp + m_params.prec() * m_params.mean());
    }
	
	
}   // end prior choice
		
    return 0;
}


double SpatialPrior::ApplyToMVN_(MVNDist *prior, const RunContext &ctx, FabberRunData &rundata,Parameter &param)
{
    prior->means(m_idx + 1) = m_params.mean();

    SymmetricMatrix prec = prior->GetPrecisions();
    prec(m_idx + 1, m_idx + 1) = m_params.prec();
    prior->SetPrecisions(prec);

    return 0;
}



std::vector<Prior *> PriorFactory::CreatePriors(const std::vector<Parameter> &params)
{
    vector<Prior *> priors;
    for (size_t i = 0; i < params.size(); i++)
    {
        Prior *prior = PriorFactory::CreatePrior(params[i]);
        LOG << "PriorFactory::CreatePriors " << *prior << endl;
        priors.push_back(prior);
    }

    return priors;
}

PriorFactory::PriorFactory(FabberRunData &rundata)
    : Loggable(rundata.GetLogger())
    , m_rundata(rundata)
{
}

Prior *PriorFactory::CreatePrior(Parameter p)
{
    LOG << "test11" << stringify(p.prior_type) << endl;
    switch (p.prior_type)
    {
    case PRIOR_NORMAL:
    case PRIOR_DEFAULT:
        return new DefaultPrior(p, m_rundata);
    case PRIOR_IMAGE:
        return new ImagePrior(p, m_rundata);
    case PRIOR_SPATIAL_M:
    case PRIOR_SPATIAL_m:
    case PRIOR_SPATIAL_P:
    case PRIOR_SPATIAL_p:
    case PRIOR_SPATIAL_n:
    case PRIOR_SPATIAL_k:
        return new SpatialPrior(p, m_rundata);
    case PRIOR_ARD:
        return new ARDPrior(p, m_rundata);
    default:
        throw InvalidOptionValue("Prior type", stringify(p.prior_type), "Supported types: NMmPpAI");
    }
}
