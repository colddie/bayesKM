// SimpleITK includes
#include <SimpleITK.h>

// ITK includes
#include "itkImage.h"
#include "itkCurvatureFlowImageFilter.h"

#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include "optim.hpp"
#include "spline.h"
// #include "tgo.h"
#include "libtpcmodel.h"

#include "meKineticRigid.h"



using namespace std;
namespace sitk = itk::simple;



// global
//
struct ll_data
{
	int verbose;
	int parNr = 6;
	int nframe;
	int fitmodel;
	int fitmethod;
	int *index;
	float *rigmotion;
	float tstart;
	float tstop;
	double *plasma_tt;
	double *plasma_t;
	double *plasma_c;
    sitk::Image imgs;
};
// Definition of IDL string
typedef struct{
	short slen;
	short stype;
	char* s;
} idls;



extern "C" double Func0(int parNr, double *vals, void *opt_data)
{
// printf("evaluate func once...");
	const char *debugfile  = "debug.txt";
	// FILE *pfile = fopen(debugfile, "a+");
	ll_data* objfn_data = reinterpret_cast<ll_data*>(opt_data);
	int verbose = objfn_data->verbose;
	int nframe = objfn_data->nframe;
	int *index = objfn_data->index; 
	float *rigmotion = objfn_data->rigmotion;
	sitk::Image imgs = objfn_data->imgs;
	std::vector<unsigned int> dims = imgs.GetSize();
	sitk::Image imgs1 = sitk::Image( dims , sitk::sitkFloat32);


//    // b-spline interpoolate rigmotion
// 	rigmotions1 = (float *) malloc(nframe*sizeof(float));
// 	rigmotions2 = (float *) malloc(nframe*sizeof(float));
// 	rigmotions3 = (float *) malloc(nframe*sizeof(float));
// 	rigmotions4 = (float *) malloc(nframe*sizeof(float));
// 	rigmotions5 = (float *) malloc(nframe*sizeof(float));
// 	rigmotions6 = (float *) malloc(nframe*sizeof(float));  
// 	for(int iframe=0;iframe<objfn_data->nframe;iframe++ ) { 
// 		int tmpIndex = index[iframe] - 1;
// 		rigmotions1[iframe] = rigmotion[tmpIndex*6];
// 		rigmotions2[iframe] = rigmotion[tmpIndex*6+1];
// 		rigmotions3[iframe] = rigmotion[tmpIndex*6+2];
// 		rigmotions4[iframe] = rigmotion[tmpIndex*6+3];
// 		rigmotions5[iframe] = rigmotion[tmpIndex*6+4];
// 		rigmotions6[iframe] = rigmotion[tmpIndex*6+5];
// 	}

//    std::vector<double> X(5), Y(5);
//    X[0]=0.1; X[1]=0.4; X[2]=1.2; X[3]=1.8; X[4]=2.0;
//    Y[0]=0.1; Y[1]=0.7; Y[2]=0.6; Y[3]=1.1; Y[4]=0.9;

//    tk::spline s;
//    s.set_points(X,Y);    // currently it is required that X is already sorted

//    printf("spline at %f is %f\n", x, s(x));

//    return EXIT_SUCCESS;

    // sitk::Image resampled_vectorout(dims,sitk::sitkVectorFloat32,objfn_data->nframe);
	std::vector<sitk::Image> resampled_vectorout;
	std::vector<unsigned int> extractSize(4);
	std::vector<int> extractIndex(4);
	// std::vector<int> destinationIndex(4);
	for(int iframe=0;iframe<objfn_data->nframe;iframe++ ) {  

		if (index[iframe] == 0) {
			extractSize = imgs.GetSize();
			extractSize[3] = 0;
			extractIndex[0] = 0;
			extractIndex[1] = 0;
			extractIndex[2] = 0;
			extractIndex[3] = iframe;
			sitk::Image out0 = sitk::Extract(imgs, extractSize,extractIndex);
			resampled_vectorout.push_back(out0);
	     	continue;
		}

		int tmpIndex = index[iframe] - 1;
		std::vector<double> rotation_center(3);
		rotation_center[0] = dims[0]/double(2);
		rotation_center[1] = dims[1]/double(2);
		rotation_center[2] = dims[2]/double(2);
		double theta_x = vals[tmpIndex*6];
		double theta_y = vals[tmpIndex*6+1];
		double theta_z = vals[tmpIndex*6+2];
		std::vector<double> translation(3);
		translation[0] = vals[tmpIndex*6+3];
		translation[1] = vals[tmpIndex*6+4];
		translation[2] = vals[tmpIndex*6+5];


		sitk::Euler3DTransform euler_transform;
		// euler_transform.SetCenter(rotation_center);     // not neccessary as auto center to image center
		euler_transform.SetRotation(theta_x, theta_y, theta_z);
		euler_transform.SetTranslation(translation);
		// resampled_image = resample(image, euler_transform)
		//   parameterMap = sitk.GetDefaultParameterMap('translation');
		//   transformixImageFilter = sitk.TransformixImageFilter();
		//   transformixImageFilter.SetTransformParameterMap(transformParameterMap);

    if (verbose == 1 ) {
		FILE *pfile = fopen(debugfile, "a+");
		fprintf(pfile, "index: %d %d, \n", iframe, tmpIndex);
		fprintf(pfile, "offset: %f, %f %f \n",  rotation_center[0], rotation_center[1], rotation_center[2]);
		fprintf(pfile, "theta: %f, %f %f \n",  theta_x, theta_y, theta_z);
		fprintf(pfile, "translation: %f, %f %f \n",  translation[0], translation[1], translation[2]);
		fclose(pfile);
	}

		extractSize = imgs.GetSize();
		extractSize[3] = 0;
		extractIndex[0] = 0;
		extractIndex[1] = 0;
		extractIndex[2] = 0;
		extractIndex[3] = iframe;

		sitk::Image out = sitk::Extract(imgs, extractSize,extractIndex);
		sitk::Image resampled_out = sitk::Resample(out, euler_transform);
        resampled_vectorout.push_back(resampled_out);
		// sitk::ImageFileWriter writer2;
		// writer2.SetFileName( std::string("slice_after_transform.nii") );
		// writer2.Execute(resampled_out);
		// sitk::ImageFileWriter writer3;
		// writer3.SetFileName( std::string("slice_before_transform.nii") );
		// writer3.Execute(out);
		// return 100;
		// sitk::Image tmpimgs = sitk::JoinSeries(resampled_out);
		// imgs1 = sitk::Paste(imgs, tmpimgs, tmpimgs.GetSize(), extractIndex);

	}   // endfor frame

    // join all 3d images to a 4d volume
	imgs1 = sitk::JoinSeries(resampled_vectorout);

    if (verbose == 1 ) {
		sitk::ImageFileWriter writer;
		writer.SetFileName( std::string("images_before_transform.nii") );
		writer.Execute(imgs);
		sitk::ImageFileWriter writer1;
		writer1.SetFileName( std::string("images_after_transform.nii") );
		writer1.Execute(imgs1);
	}

	// Patlak
    std::vector<unsigned int> dim = { dims[0],dims[1],dims[2] };
	sitk::Image Kiimg = sitk::Image( dim, sitk::sitkFloat32); 
	double var = 0.0;
	double output[5];
	double *weights = (double *)malloc(nframe*sizeof(double));
	// double *plasma_t = reinterpret_cast<double*>(objfn_data->plasma_t);
	// double *plasma_tt = reinterpret_cast<double*>(objfn_data->plasma_tt);
	// double *plasma_c = reinterpret_cast<double*>(objfn_data->plasma_c);
	for (int iframe=0;iframe<nframe;iframe++) { weights[iframe] = 1.0; }
	for (int jplane=1;jplane<dims[2];jplane++)	{    //evrey other two
	for (int jrow=0;jrow<dims[1]-6;jrow+=6) {
	for (int jcol=0;jcol<dims[0]-6;jcol+=6) 
	{
		// if ( (jrow != 120) && (jcol != 120) ) { continue; }
		double *tac = (double *)malloc(nframe*sizeof(double));    // cast float to double?
		// sitk::Image out = sitk::Extract(imgs1,extractSize,extractIndex);
		std::vector<unsigned int> tacIndex(4);
		double tacsum = 0.0;
		for (int iframe=0;iframe<nframe;iframe++) { 
			tacIndex[0] = jcol;
			tacIndex[1] = jrow;
			tacIndex[2] = jplane;
			tacIndex[3] = iframe;
			tac[iframe]=(double)imgs1.GetPixelAsFloat(tacIndex); 
			tacsum +=tac[iframe];
		}
        if (tacsum < 0.1) { continue ;}      // skip the background

		int success = patlak_c(nframe, objfn_data->plasma_tt, objfn_data->plasma_t, tac, 
						objfn_data->plasma_c,(double)objfn_data->tstart,(double)objfn_data->tstop,output,0,0,0,weights);    //debug,llsq_model,isweight
		var += abs(output[2])+abs(output[3]);   // absolute

		if (verbose==1) { 
			// if (output[0] > 0.0) { printf("positive detected! \n"); }
			              std::vector<unsigned int> vIndex(3);
		                  vIndex[0]=jcol; vIndex[1]=jrow; vIndex[2]=jplane; 
						  Kiimg.SetPixelAsFloat(vIndex,output[1]); }

		if (verbose==1) { 
			FILE *pfile = fopen(debugfile, "a+");
			fprintf(pfile, "pixel index: %d %d %d \n", jcol, jrow, jplane);
			fprintf(pfile, "offset: %f %f %f %f %f %f \n", output[0], output[1], output[2],
			        output[3], output[4], output[5]);
			fprintf(pfile, "offset:  %f %f %f %f %f \n", tac[0], tac[1], tac[2], tac[nframe-2],tac[nframe-1] );
			fclose(pfile);
		}
	    free(tac);		
	}
	}
	}

	if (verbose == 1) {		
		sitk::ImageFileWriter writerK;
		writerK.SetFileName( std::string("Kiimg.nii") );
		writerK.Execute(Kiimg);
	}
    free(weights);	
	// FILE *pfile = fopen(debugfile, "a+");
	printf("total variance: %f current estimate %f %f %f %f %f %f \n", var,
	                vals[0],vals[1],vals[2],vals[3],vals[4],vals[5]);
	return var;
}



extern "C" double Func (const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)
{
	int parNr = 6;
	double vals[6];
	for (int ipar=0;ipar<parNr;ipar++)  { vals[ipar] = vals_inp[ipar]; }
	double var = Func0(parNr, vals, opt_data);

	return var;
}





extern "C"  int meKineticRigid(int argc, float* argv[])     //meKineticRigid
{

    std::string imgfilename;
	float *parms0;
	int fitmodel;
	int   *index;
	const char *debugfile  = "debug.txt"; 
	// int verbose = 1; 
    // float *rigotion;

    /* debug */
    if (argc > 13 || argc <= 1) {
      FILE *pfile = fopen(debugfile, "a+");
	  fprintf(pfile, "meKineticRigid: 18 arguments required, %d supplied\n", argc);
	  fclose(pfile);
      return -1;
    }

	// Read the image
	//	
	ll_data opt_data;
    opt_data.nframe      =  *(size_t *)  argv[0];
	imgfilename          =  std::string((*(idls *) argv[1]).s);
	parms0               =  (float *)  argv[2];		
	opt_data.tstart      =  *(float *) argv[3];
	opt_data.tstop       =  *(float *) argv[4];
	opt_data.plasma_tt   =  (double *)  argv[5];
	opt_data.plasma_t    =  (double *)  argv[6];
	opt_data.plasma_c    =  (double*)   argv[7];
	opt_data.fitmodel    =  *(int *)   argv[8];
	opt_data.fitmethod   =  *(int *)   argv[9];
	opt_data.rigmotion   =  (float *)  argv[10];
	opt_data.index       =  (int *)   argv[11];
	opt_data.verbose     =  *(int *)   argv[12];

	sitk::ImageFileReader reader;
	reader.SetFileName( imgfilename );
	opt_data.imgs = reader.Execute();

    if (opt_data.verbose == 1 ) {
		FILE *pfile = fopen(debugfile, "a+");
		fprintf(pfile, "test, %d supplied\n", argc);
		fprintf(pfile, "size: %d, %f %f \n" , opt_data.nframe , opt_data.tstart, opt_data.tstop);
		fprintf(pfile, "spacing: %f, %f %f \n",  opt_data.plasma_tt[0], opt_data.plasma_t[0], opt_data.plasma_c[opt_data.nframe-1]);
		// fprintf(pfile, "offset: %f, %f %f \n",  offsetx, offsety, offsetz);
		fclose(pfile);
	}
	


if (opt_data.fitmethod == 1) {
    arma::vec x = arma::ones(opt_data.parNr,1); // initial values: (2,2)
	for (int i=0;i<opt_data.parNr;i++) {
		x[i] = parms0[i]; }

    optim::algo_settings_t settings;
	arma::vec upper_bounds(opt_data.parNr);
	arma::vec lower_bounds(opt_data.parNr);
	lower_bounds[0] = 0.0;	upper_bounds[0] =  0.0;
	lower_bounds[1] = 0.0;	upper_bounds[1] =  0.0;
	lower_bounds[2] = -0.3;	upper_bounds[2] =  0.3;
	lower_bounds[3] = 0.0;	upper_bounds[3] =  0.0;	
	lower_bounds[4] = -8.0;	upper_bounds[4] =  8.0;	
	lower_bounds[5] = -8.0;	upper_bounds[5] =  8.0;		

	settings.vals_bound = true;
    settings.lower_bounds = lower_bounds;
	settings.upper_bounds = upper_bounds;
    settings.verbose_print_level = 1;   //opt_data.verbose;
	settings.de_max_fn_eval = 100;
	settings.de_check_freq = 100;
    int success = optim::de(x,Func,&opt_data,settings);

    if (success) {
        std::cout << "de: Booth test completed successfully." << std::endl;
    } else {
        std::cout << "de: Booth test completed unsuccessfully." << std::endl;
    }
    arma::cout << "de: solution to Booth test:\n" << x << arma::endl;
}

if (opt_data.fitmethod==2) {

    int tgoNr,iterNr,neighNr,parNr;
    double pmin[6], pmax[6];
    bool TGO_LOCAL_INSIDE=0;
    bool TGO_SQUARED_TRANSF=1;
    tgoNr=300; iterNr=0; neighNr=5;
    parNr=6;
    iterNr=0;
	pmin[0] = 0.0;	pmax[0] =  0.0;
	pmin[1] = 0.0;	pmax[1] =  0.0;
	pmin[2] = -0.3;	pmax[2] =  0.3;
	pmin[3] = 0.0;	pmax[3] =  0.0;	
	pmin[4] = -8.0;	pmax[4] =  8.0;
	pmin[5] = -8.0; pmax[5] =  8.0;

    double wss=0;
    double *output = (double *)malloc(parNr*sizeof(double));   // not needed!
    // double output[2];

    // bool success_0 = 1;
    bool success_0 = tgo(
      pmin, pmax, Func0, &opt_data, parNr, neighNr,
      &wss, output, tgoNr, iterNr, 10);  // opt_data.verbose);

    if (!success_0) {
        std::cout << "powell: Booth test completed successfully." << std::endl;
    } else {
        std::cout << "powell: Booth test completed unsuccessfully." << std::endl;
    }
    
    std::cout << "powell: solution to Booth test: \n"  << output[0] <<" " << output[1] << " "<<
	    output[2] << " "<< output[3] << " " << output[4] << " "<<output[5]<< "\n" << std::endl;
	free(output);
}


    return 1;
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