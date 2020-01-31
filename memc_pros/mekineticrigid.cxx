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
	int parNr;    // = 6;
	int nframe;
	int slowmotion;
	int fitmodel;
	int fitmethod;
	int *index;
	int *Findex;
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
	ll_data* objfn_data = reinterpret_cast<ll_data*>(opt_data);
	int verbose = objfn_data->verbose;
	int slowmotion = objfn_data->slowmotion;
	int nframe = objfn_data->nframe;
	// int parNr  = objfn_data->parNr;
	int pparNr = (objfn_data->parNr)/6;
	int *index = objfn_data->index; 
	int *Findex = objfn_data->Findex;
	float *rigmotion = objfn_data->rigmotion;
	sitk::Image imgs = objfn_data->imgs;
	std::vector<unsigned int> dims = imgs.GetSize();
	sitk::Image imgs1 = sitk::Image( dims , sitk::sitkFloat32);

	std::vector<float> rigmotions(nframe*6);    //nframe*6
	std::vector<double> X(pparNr),Y1(pparNr),Y2(pparNr),Y3(pparNr),Y4(pparNr),Y5(pparNr),Y6(pparNr);
    if (slowmotion) {
		// b-spline interpoolate rigmotion
		std::vector<float> rigmotions1(nframe); 
		std::vector<float> rigmotions2(nframe); 
		std::vector<float> rigmotions3(nframe); 
		std::vector<float> rigmotions4(nframe); 
		std::vector<float> rigmotions5(nframe); 
		std::vector<float> rigmotions6(nframe); 

		for(int iframe=0;iframe<nframe;iframe++ ) { 
			int tmpIndex = index[iframe];     // - 1;
			rigmotions1[iframe] = vals[tmpIndex*6];
			rigmotions2[iframe] = vals[tmpIndex*6+1];
			rigmotions3[iframe] = vals[tmpIndex*6+2];
			rigmotions4[iframe] = vals[tmpIndex*6+3];
			rigmotions5[iframe] = vals[tmpIndex*6+4];
			rigmotions6[iframe] = vals[tmpIndex*6+5]; printf("%d \n", iframe);
		}

		// interpoalte for each individual frame
		for (int ipar;ipar<pparNr;ipar++  ) {
			int tmpIndex2 = Findex[ipar];      // - 1;
			X[ipar] = objfn_data->plasma_t[tmpIndex2]; 
			Y1[ipar] = rigmotions1[tmpIndex2]; 
			Y2[ipar] = rigmotions2[tmpIndex2]; 
			Y3[ipar] = rigmotions3[tmpIndex2]; 
			Y4[ipar] = rigmotions4[tmpIndex2]; 
			Y5[ipar] = rigmotions5[tmpIndex2]; 
			Y6[ipar] = rigmotions6[tmpIndex2]; 
		}

		tk::spline s1; s1.set_points(X,Y1); tk::spline s2; s2.set_points(X,Y2);
		tk::spline s3; s3.set_points(X,Y3); tk::spline s4; s4.set_points(X,Y4);
		tk::spline s5; s5.set_points(X,Y5); tk::spline s6; s6.set_points(X,Y6);

		for (int iframe=0;iframe<nframe;iframe++) {
			rigmotions[iframe*6]   = s1(objfn_data->plasma_t[iframe]);
			rigmotions[iframe*6+1] = s2(objfn_data->plasma_t[iframe]);
			rigmotions[iframe*6+2] = s3(objfn_data->plasma_t[iframe]);
			rigmotions[iframe*6+3] = s4(objfn_data->plasma_t[iframe]);
			rigmotions[iframe*6+4] = s5(objfn_data->plasma_t[iframe]);
			rigmotions[iframe*6+5] = s6(objfn_data->plasma_t[iframe]);printf("%d \n", iframe);
		}
	}

    // sitk::Image resampled_vectorout(dims,sitk::sitkVectorFloat32,objfn_data->nframe);
	std::vector<sitk::Image> resampled_vectorout;
	std::vector<unsigned int> extractSize(4);
	std::vector<int> extractIndex(4);
	// std::vector<int> destinationIndex(4);
	for(int iframe=0;iframe<objfn_data->nframe;iframe++ ) {  
// printf("%d \n", iframe);
		if (index[iframe] == 0) {
			extractSize = imgs.GetSize();
			extractSize[3]  = 0;
			extractIndex[0] = 0;
			extractIndex[1] = 0;
			extractIndex[2] = 0;
			extractIndex[3] = iframe;
			sitk::Image out0 = sitk::Extract(imgs, extractSize,extractIndex);
			resampled_vectorout.push_back(out0);
	     	continue;
		}

        // get motion for each frame
		double theta_x;
		double theta_y;
		double theta_z;
		std::vector<double> translation(3);
		std::vector<double> rotation_center(3);
		rotation_center[0] = dims[0]/double(2);
		rotation_center[1] = dims[1]/double(2);
		rotation_center[2] = dims[2]/double(2);
		if (slowmotion) {
			theta_x        = rigmotions[iframe*6];
			theta_y        = rigmotions[iframe*6+1];
			theta_z        = rigmotions[iframe*6+2];
			translation[0] = rigmotions[iframe*6+3];
			translation[1] = rigmotions[iframe*6+4];
			translation[2] = rigmotions[iframe*6+5];
			// print rigmotions[]
		} else {
			int tmpIndex   = index[iframe];   /// - 1;
			theta_x        = vals[tmpIndex*6];
			theta_y        = vals[tmpIndex*6+1];
			theta_z        = vals[tmpIndex*6+2];
			translation[0] = vals[tmpIndex*6+3];
			translation[1] = vals[tmpIndex*6+4];
			translation[2] = vals[tmpIndex*6+5];
		}

		// sitk::AffineTransform euler_transform(3);
		// euler_transform.SetCenter(rotation_center);  
		// std::vector<double> matrix0 = euler_transform.GetMatrix();
		// std::vector<double> matrix = { cos(theta_z),-sin(theta_z),0,sin(theta_z),cos(theta_z),0,0,0,1 }; 
        // euler_transform.SetMatrix(matrix);

		sitk::Euler3DTransform euler_transform;
		// euler_transform.SetCenter(rotation_center);     // not neccessary as auto center to image center
		euler_transform.SetRotation(theta_x, theta_y, theta_z);		
		euler_transform.SetTranslation(translation);



    if (verbose == 1 ) {
		FILE *pfile = fopen(debugfile, "a+");
		// fprintf(pfile, "motions: %f %f %f %f %f %f %f %f \n", rigmotions[0],rigmotions[(nframe-1)*6+5], 
		//                        X[0], X[pparNr-1], Y1[0], Y1[pparNr-1], rigmotions1[0],rigmotions1[nframe-1]);
		fprintf(pfile, "nframe: %d , \n", nframe);		
		fprintf(pfile, "index: %d , \n", iframe);
		fprintf(pfile, "offset: %f, %f %f \n",  rotation_center[0], rotation_center[1], rotation_center[2]);
		fprintf(pfile, "theta: %f, %f %f \n",  theta_x, theta_y, theta_z);
		fprintf(pfile, "translation: %f, %f %f \n",  translation[0], translation[1], translation[2]);
		fclose(pfile);
	}

		extractSize     = imgs.GetSize();
		extractSize[3]  = 0;
		extractIndex[0] = 0;
		extractIndex[1] = 0;
		extractIndex[2] = 0;
		extractIndex[3] = iframe;

		sitk::Image out = sitk::Extract(imgs, extractSize,extractIndex);
		euler_transform.SetCenter(out.TransformContinuousIndexToPhysicalPoint(rotation_center));  //perform transformation and resampling in physical space and not index space
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
	double output[5];  output[0]=0; output[1]=0; output[2]=0; output[3]=0; output[4]=0;    
	double *weights = (double *)malloc(nframe*sizeof(double));
	// double *plasma_t = reinterpret_cast<double*>(objfn_data->plasma_t);
	// double *plasma_tt = reinterpret_cast<double*>(objfn_data->plasma_tt);
	// double *plasma_c = reinterpret_cast<double*>(objfn_data->plasma_c);
	for (int iframe=0;iframe<nframe;iframe++) { weights[iframe] = 1.0; }
	for (int jplane=0;jplane<dims[2];jplane++)	{    //evrey other two
	for (int jrow=0;jrow<dims[1];jrow+=1) {
	for (int jcol=0;jcol<dims[0];jcol+=1) 
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
        if (tacsum < 0.1) { 	    free(tac);	continue ;}      // skip the background, do not forget free memory

		int success = patlak_c(nframe, objfn_data->plasma_tt, objfn_data->plasma_t, tac, 
						objfn_data->plasma_c,(double)objfn_data->tstart,(double)objfn_data->tstop,output,0,0,0,weights);    //debug,llsq_model,isweight
		var += abs(output[2]);   //+abs(output[3]);   // absolute
        
		if (success != 1) { printf("patlak fitting error! \n"); }
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

	if (verbose == 1) {			
		// FILE *pfile = fopen(debugfile, "a+");
		// printf("total variance: %f current estimate %f %f %f %f %f %f \n", var,
		//                 vals[0],vals[1],vals[2],vals[3],vals[4],vals[5]);
		printf("total variance: %f  ",var);
		printf("current estimate ");
		for(int i=0;i<objfn_data->parNr;i++) { printf("%f ", vals[i]); }
		printf("\n");
	}

	free(weights);	
	return var;
}



extern "C" double Func (const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data)
{

	// int parNr = 6;
	ll_data* objfn_data = reinterpret_cast<ll_data*>(opt_data);
	double vals[objfn_data->parNr];
	for (int ipar=0;ipar<objfn_data->parNr;ipar++)  { vals[ipar] = vals_inp[ipar]; }
	double var = Func0(objfn_data->parNr, vals, objfn_data);

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
    if (argc > 16 || argc <= 1) {
      FILE *pfile = fopen(debugfile, "a+");
	  fprintf(pfile, "meKineticRigid: 18 arguments required, %d supplied\n", argc);
	  fclose(pfile);
      return -1;
    }

	// Read the image
	//	
	ll_data opt_data;
	opt_data.parNr       =  *(size_t *)  argv[0];
    opt_data.nframe      =  *(size_t *)  argv[1];
	imgfilename          =  std::string((*(idls *) argv[2]).s);
	parms0               =  (float *)  argv[3];		
	opt_data.tstart      =  *(float *) argv[4];
	opt_data.tstop       =  *(float *) argv[5];
	opt_data.plasma_tt   =  (double *)  argv[6];
	opt_data.plasma_t    =  (double *)  argv[7];
	opt_data.plasma_c    =  (double *)   argv[8];
	opt_data.fitmodel    =  *(int *)   argv[9];
	opt_data.fitmethod   =  *(int *)   argv[10];
	opt_data.rigmotion   =  (float *)  argv[11];
	opt_data.index       =  (int *)   argv[12];
	opt_data.Findex      =  (int *)   argv[13];
	opt_data.verbose     =  *(int *)   argv[14];
	opt_data.slowmotion  =  *(int *)   argv[15];




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
	lower_bounds[3] = -8.0;	upper_bounds[3] =  8.0;	
	lower_bounds[4] = -8.0;	upper_bounds[4] =  8.0;	
	lower_bounds[5] = 0.0;	upper_bounds[5] =  0.0;		
	// lower_bounds[6]

	settings.vals_bound     = true;
    settings.lower_bounds   = lower_bounds;
	settings.upper_bounds   = upper_bounds;
    settings.verbose_print_level = 1;   //opt_data.verbose;
	settings.de_max_fn_eval = 100;
	settings.de_check_freq  = 100;
    int success = optim::de(x,Func,&opt_data,settings);

    if (success) {
        std::cout << "de: Booth test completed successfully." << std::endl;
    } else {
        std::cout << "de: Booth test completed unsuccessfully." << std::endl;
    }
    arma::cout << "de: solution to Booth test:\n" << x << arma::endl;
}

if (opt_data.fitmethod==2) {

    int tgoNr,iterNr,neighNr;   //,parNr;
    double pmin[opt_data.parNr], pmax[opt_data.parNr];
    bool TGO_LOCAL_INSIDE=0;
    bool TGO_SQUARED_TRANSF=1;
    tgoNr=30; iterNr=10; neighNr=5;
	for (int i=0;i<opt_data.parNr;i++) { 
		if (i%6==0 || i%6==1 || i%6==5) { pmin[i] = 0.0;	pmax[i] =  0.0; } 
		if (i%6==3 || i%6==4)           { pmin[i] = -10.0;	pmax[i] =  10.0; } 
		if (i%6==2)                     { pmin[i] = -0.3;	pmax[i] =  0.3; } 
		// printf("pmin and pmax at voxel %f %f",pmin[i],pmax[i]);			
	}

    double wss=0;
    // double *output = (double *)malloc(opt_data.parNr*sizeof(double));   // not needed!
    double output[opt_data.parNr];

    bool success_0 = tgo(
      pmin, pmax, Func0, &opt_data, opt_data.parNr, neighNr,
      &wss, output, tgoNr, iterNr, 10);  // opt_data.verbose);

    if (!success_0) {
        std::cout << "powell: Booth test completed successfully." << std::endl;
    } else {
        std::cout << "powell: Booth test completed unsuccessfully." << std::endl;
    }
    
    // std::cout << "powell: solution to Booth test: \n"  << output[0] <<" " << output[1] << " "<<
	//     output[2] << " "<< output[3] << " " << output[4] << " "<<output[5]<< "\n" << std::endl;
	std::cout << "powell: solution to Booth test: \n"<< std::endl;
	for (int i=0;i<opt_data.parNr;i++) { std::cout <<output[i]<<" "<< std::endl; } 
	// free(output);
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