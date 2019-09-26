/*  inference_spatialvb.h - implementation of VB with spatial priors

 Adrian Groves & Michael Chappell, FMRIB Image Analysis Group & IBME QuBIc Group

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

#include "convergence.h"
#include "inference.h"
#include "run_context.h"
#include "newimage/newimage.h"
#include <newmat.h>

// #include "nlm.cxx"

#include <string>
#include <vector>

class Vb : public InferenceTechnique
{
public:
    static InferenceTechnique *NewInstance();

    Vb()
        : m_nvoxels(0)
        , m_noise_params(0)
        , m_needF(false)
        , m_printF(false)
        , m_saveF(false)
        , m_origdata(NULL)
        , m_coords(NULL)
        , m_suppdata(NULL)
        // , m_pimg (NULL)
        , m_num_mcsteps(0)
        , m_spatial_dims(-1)
        , m_locked_linear(false)
    {
    }

    virtual void GetOptions(vector<OptionSpec> &opts) const;
    virtual std::string GetDescription() const;
    virtual string GetVersion() const;

    virtual void Initialize(FwdModel *fwd_model, FabberRunData &args);
    virtual void DoCalculations(FabberRunData &data);

    virtual void SaveResults(FabberRunData &rundata) const;

protected:
    /**
     * Initialize noise prior or posterior distribution from a file stored in the
     * rundata under the given parameter key
     */
    void InitializeNoiseFromParam(FabberRunData &args, NoiseParams *dist, string param_key);

    /**
     * Pass the model the data, coords and suppdata for a voxel.
     *
     * FIXME this is not very nice and should not be necessary. Need to
     * audit what models are using this info and find alternatives, e.g.
     * reading suppdata in Initialize instead
     */
    void PassModelData(int voxel);

    /**
     * Determine whether we need spatial VB mode
     *
     * It is required either because it has been asked for (--method=spatialvb) or
     * if any spatial priors have been specified (types mMpP)
     */
    bool IsSpatial(FabberRunData &rundata) const;

    /**
     * Do calculations loop in voxelwise mode (i.e. all iterations for
     * one voxel, then all iterations for the next voxel, etc)
     */
    virtual void DoCalculationsVoxelwise(FabberRunData &data);

    /**
     * Do calculations loop in spatial mode (i.e. one iteration of all
     * voxels, then next iteration of all voxels, etc)
     */
    virtual void DoCalculationsSpatial(FabberRunData &data);

    /**
     * Calculate free energy if required, and display if required
     */
    double CalculateF(int v, std::string label, double Fprior);

    /**
     * Output detailed debugging information for a voxel
     */
    void DebugVoxel(int v, const string &where);

    /**
     * Setup per-voxel data for Spatial VB
     *
     * Spatial VB needs each voxel's prior/posterior and other
     * data stored as it affects neighbouring voxels. This sets
     * up the vectors which store these things which are just
     * created on the fly for normal VB and throw away after each
     * voxel is done.
     */
    void SetupPerVoxelDists(FabberRunData &allData);

    /**
    * Check voxels are listed in order
    *
    * Order must be increasing in z value, or if same
    * increasing in y value, and if y and z are same
    * increasing in x value.
    *
    * This is basically column-major (Fortran) ordering - used as default by NEWIMAGE.
    */
    void CheckCoordMatrixCorrectlyOrdered(const NEWMAT::Matrix &voxelCoords);

    /**
     * Calculate first and second nearest neighbours of each voxel
    */
    void CalcNeighbours(const NEWMAT::Matrix &voxelCoords);

    /**
     * Ignore this voxel in future updates.
     *
     * No calculation of priors or posteriors will occur for this voxel
     * and it will be removed from the lists of neighbours for other voxels.
     * The effect should be as if it were masked
     */
    void IgnoreVoxel(int v);

    /** Number of voxels in data */
    int m_nvoxels;

    /**
     * Noise model in use. This is created by the inference
     * method deleted on destruction
     */
    std::auto_ptr<NoiseModel> m_noise;
    /**
     * Number of noise parameters.
     */
    int m_noise_params;

    /** True if convergence detector requires the free energy */
    bool m_needF;

    /** True if we need to print the free energy at each iteration */
    bool m_printF;

    /** True if we need to to save the final free energy */
    bool m_saveF;

    /** Free energy for each voxel */
    std::vector<double> resultFs;

    /** Voxelwise input data */
    const NEWMAT::Matrix *m_origdata;

    /** Voxelwise co-ordinates */
    const NEWMAT::Matrix *m_coords;

    /** Voxelwise supplementary data */
    const NEWMAT::Matrix *m_suppdata;
	
	// //////////////////////////////////////////
    // /** Voxelwise prior data */	
	// const NEWMAT::Matrix *m_pimg;
	// /////////////////////////////////////////

    /** Number of motion correction steps to run */
    int m_num_mcsteps;

    /** Stores current run state (parameters, MVNs, linearization centres etc */
    RunContext *m_ctx;

    /** Linearized wrapper around the forward model */
    std::vector<LinearizedFwdModel> m_lin_model;

    /** Convergence detector for each voxel */
    std::vector<ConvergenceDetector *> m_conv;

    /**
     * Number of spatial dimensions
     *
     * 0 = no spatial smoothing
     * 1 = Probably not sensible!
     * 2 = Smoothing in slices only
     * 3 = Smoothing by volume
     */
    int m_spatial_dims;

    /**
     * Fix the linearization centres of the linearized forward model.
     *
     * This reduces the inference to a purely linear problem. The fixed
     * centres are generally loaded from an MVN file
     */
    bool m_locked_linear;

    bool m_have_mask;

    NEWMAT::RowVector m_mask;
};
