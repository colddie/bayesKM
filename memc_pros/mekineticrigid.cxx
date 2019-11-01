// SimpleITK includes
#include "SimpleITK.h"

// ITK includes
#include "itkImage.h"
#include "itkCurvatureFlowImageFilter.h"

#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <armadillo>
#include "optim.hpp"
// #include "tgo.h"

#include 

using namespace std;
namespace sitk = itk::simple;



// global
//
struct ll_data
{
	int parNr = 6;
	int nframe;
	float *rigmotion;
	float tstart;
	float tstop;
	float *plasma_tt;
	float *plasma_t;
	float *plasma_c;
    sitk::Image images;
};



int meKineticRigid(int argc, float* argv[])
{
	const char *imgfilename;
	float *parms0;
	int model;
	const char *debugfile  = "debug.txt"; 
	int verbose = 1; 
    // float *rigotion;


    /* debug */
    if (argc > 9 || argc <= 1) {
      FILE *pfile = fopen(debugfile, "a+");
	  fprintf(pfile, "meKineticRigid: 18arguments required, %d supplied\n", argc);
	  fclose(pfile);
      return -1;
    }


	// Read the image
	//	
	ll_data opt_data;
    optdata.nframe       =  *(size_t *)  argv[0];
	imgfilename          =  string((*(idls *) argv[1]).s);
	parms0               =  (float *)  argv[2];		
	opt_data.tstart      =  *(float *) argv[3];
	opt_data.tstop       =  *(float *) argv[4];
	opt_data.plasma_tt   =  (float *)  argv[5];
	opt_data.plasma_t    =  (float *)  argv[6];
	opt_data.plasma_c    =  *(unsigned int*) argv[7];
	opt_data.rigmotion   =  (float *)  argv[8];

	sitk::ImageFileReader reader;
	reader.SetFileName( imgfilename );
	opt_data.images = reader.Execute();

    if (verbose == 1 ) {
		FILE *pfile = fopen(debugfile, "a+");
		fprintf(pfile, "test, %d supplied\n", argc);
		fprintf(pfile, "size: %d, %d %d \n" , optdata.nframe , opt_data.tstart, opt_data.tstop);
		fprintf(pfile, "spacing: %f, %f %f \n",  opt_data.plasma_tt[0], opt_data.plasma_t[0], opt_data.plasma_c[0]);
		// fprintf(pfile, "offset: %f, %f %f \n",  offsetx, offsety, offsetz);
		fclose(pfile);
	}
	

    arma::vec x = arma::ones(opt_data.parNr,1); // initial values: (2,2)
	for (int i=0;i<opt_data.parNr;i++) {
		x[i] = opt_data.parms0[i]; }

    settings.verbose_print_level = verbose;
    int success = optim::de(x,Func,opt_data,settings);
 
    if (success) {
        std::cout << "de: Booth test completed successfully." << std::endl;
    } else {
        std::cout << "de: Booth test completed unsuccessfully." << std::endl;
    }
    arma::cout << "de: solution to Booth test:\n" << x << arma::endl;


    return 1;
}





float Func (const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)
{

	ll_data_t* objfn_data = reinterpret_cast<ll_data_t*>(opt_data);
	float *rigmotion = objfn_data.rigmotion;
	sitk::Image imgs = objfn_data.imgs;
	std::vector<unsigned int> dims = imgs.GetDimension();
	sitk::Image imgs1 = sitk::Image( dims , sitk::sitkFloat32);

	//
	std::vector<unsigned int> rotation_center(3);
	rotation_center[0] = dims[0]/2;
	rotation_center[1] = dims[1]/2;
	rotation_center[2] = dims[2]/2;
	float theta_x = rigmotion[0];
	float theta_y = rigmotion[1];
	float theta_z = rigmotion[2];
	std::vector<unsigned int> translation(3);
	translation[0] = rigmotion[3];
	translation[1] = rigmotion[4];
	translation[2] = rigmotion[5];

	sitk::Euler3DTransform euler_transform;
	//  = sitk::Euler3DTransform(rotation_center,
	//  theta_x, theta_y, theta_z, translation)
	euler_transform.SetCenter(rotation_center);
	euler_transform.SetRotation(theta_x, theta_y, theta_z);
	euler_transform.SetOffet(translation);
	// resampled_image = resample(image, euler_transform)
	//   parameterMap = sitk.GetDefaultParameterMap('translation');
	//   transformixImageFilter = sitk.TransformixImageFilter();
	//   transformixImageFilter.SetTransformParameterMap(transformParameterMap);


	for(int iframe=objfn_data.nframe/2;iframe<objfn_data.nframe;iframe++ ) {   ///?????
		// tmpimg = imgs[:,:,:iframe];
		std::vector<unsigned int> extractSize(4);
		std::vector<unsigned int> extractIndex(4);
		extractSize = imgs.Getsize();
		extractSize[3] = 0;
		extractIndex[0] = 0;
		extractIndex[1] = 0;
		extractIndex[2] = 0;
		extractIndex[3] = iframe;

		sitk::Image out = sitk::Extract(imgs, extractSize,extractIndex);
		sitk::Image resampled_out = sitk::resample(out, euler_transform);

		imgs1 = sitk::Paste(imgs, resampled_out, resampled_out.GetSize(), destinationIndex=[0,0,0,iframe]);
	}
	sitk::writeImage(imgs,'images_before_transform.nii');
	sitk::writeImage(imgs1,'images_afer_transform.nii');

	// Patlak
	float var = 0.0;
	float *tac = (float *)malloc(nframe*sizeof(float));    // cast float to double?
	float *output = (float *)malloc(5*sizeof(float));
	for (int jcol=0;jcol<dims[0];jcol++) {
	for (int jrow=0;jrow<dims[1];jrow++) {
	for (int jplane=0;jplane<dims[2];jplane++)
	{
		extractSize[0] = 1;
		extractSize[1] = 1;
		extractSize[2] = 1;
		extractSize[3] = nframe;
		extractIndex[0] = jcol;
		extractIndex[1] = jrow;
		extractIndex[2] = jplane;
		extractIndex[3] = iframe;
		sitk::Image out = sitk::Extract(imgs1,extractSize,extractIndex);

		for (int iframe=0;iframe<nframe;iframe++) { tac[iframe]=out[iframe]; }

		// success = patlak_c(nframe, objfn_data.plasma_t, objfn_data.plasma_tt, tac, 
		// 				objfn_data.plasma_c,objfn_data.tstart,objfn_data.tstop,output,0,0,0);    //debug,llsq_model,isweight
		var += output[4];
	}

	free(tac);
	free(output);

	return var;

}
}
}


// ImportFilterType readinFilter_voxel(fiximg,defimg,sizex,sizey,sizez,spacingx,spacingy
//                  spacingz,offsetx,offsety,offsetz,outputdir,paramfile)
// {

// 	/*  Here start the ITK part... */
//     // ------------------------------------------------------
// 	// get parameters ...
//     ImportFilterType::Pointer importerfix = ImportFilterType::New();

// 	ImageType::SizeType size;
// 	size[0] = sizex; //406;        //int(dim[0]);
// 	size[1] = sizey; //408;        //int(dim[1]);
// 	size[2] = sizez;

// 	ImageType::IndexType start;
// 	start[0] = 0;
// 	start[1] = 0;
// 	start[2] = 0;

// 	ImageType::RegionType region;
// 	region.SetSize(size);
// 	region.SetIndex(start);
// 	importerfix->SetRegion(region);

// 	double spacing[3];
// 	spacing[0] = spacingx;
// 	spacing[1] = spacingy;
// 	spacing[2] = spacingz;
// 	importerfix->SetSpacing(spacing);

// 	double origin[3];
// 	origin[0] = offsetx;
// 	origin[1] = offsety;
// 	origin[2] = offsetz;
// 	importerfix->SetOrigin(origin);

// 	const bool importFilterWillDeleteTheInputBuffer = false;
	
// 	PixelType * pixelDatafix = static_cast< PixelType * > (fiximg); 
// 	const unsigned int totalNumberOfPixels = size[0] * size[1] * size[2];
	
// 	importerfix->SetImportPointer(pixelDatafix, 
// 	                           totalNumberOfPixels,
// 							   importFilterWillDeleteTheInputBuffer);
// 	importerfix->Update();

//     return importerfix;

// }