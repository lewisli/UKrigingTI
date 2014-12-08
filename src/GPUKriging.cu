/**
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */

#include "GPUKriging.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include "cutil_math.h"
#include "vector_functions.h"

// A wrapper for making GPU calls that lists the error and location
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, char *file, int line, bool abort=true)
{
	if (code != cudaSuccess)
	{
		fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

bool lessThanX(const datConfig &elem1, const datConfig &elem2)
{
	return elem1.coord.x < elem2.coord.x;
}
bool lessThanY(const datConfig &elem1, const datConfig &elem2)
{
	return elem1.coord.y < elem2.coord.y;
}

__device__ int sub2ind(short2 coord, short2 dim)
{
	return coord.y*dim.x + coord.x;
}

__global__ void computeSOP(float **dTIValues, short2 *dh_alpha, short2 *TISize,
		short2 *erodedSize, short2 *erodedStart, float scaleFactor,
		float *output, int numImages)
{
	// Shared memory (closest power of 2)
	extern __shared__ float sData[];
	sData[threadIdx.x] = 0;

	short2 index0,index1,index2;

	// Iterate over each training image
	for (unsigned int n = 0; n < numImages; n++)
	{
		// Read current training image size and start coordinate
		short2 currentImageSz = erodedSize[n];
		short2 currentImageStart = erodedStart[n];

		// Current training image sum
		float sum = 0;

		// Get pointer to actual training image
		float *currentImage = dTIValues[n];

		if (threadIdx.x < currentImageSz.x)
		{
			// We use each thread to compute each row of all training images
			for (unsigned int i = 0; i < currentImageSz.y; ++i)
			{
				// Index0 has coordinates (threadIdx.x, i)
				index0.x = threadIdx.x + currentImageStart.x;
				index0.y = i + currentImageStart.y;

				index1 = index0 + dh_alpha[blockIdx.x];
				index2 = index0 + dh_alpha[blockIdx.y];


				sum += currentImage[sub2ind(index1,TISize[n])]*
						currentImage[sub2ind(index2,TISize[n])];
			}

			sData[threadIdx.x] += sum*scaleFactor;
		}
		else
			sData[threadIdx.x] += 0;
	}

	__syncthreads();
	float globalSum = 0;


	if (threadIdx.x == 0)
	{
		for (unsigned int i = 0; i < blockDim.x; i++)
			globalSum += sData[i];

		int outputIndex = blockIdx.x*gridDim.y + blockIdx.y;
		output[outputIndex] = globalSum;
	}
}

__global__ void computeSOP0AVG(float **dTIValues, short2 *dh_alpha,
		short2 *TISize, short2 *erodedSize, short2 *erodedStart,
		float scaleFactor, float *SOP0, float *AVG, int numImages)
{
	// Shared memory
	extern __shared__ float sData[];
	sData[threadIdx.x] = 0;
	sData[threadIdx.x+blockDim.x] = 0;

	short2 index0,index1;

	for (unsigned int n = 0; n < numImages; n++)
	{
		// Read current training image size and start coordinate
		short2 currentImageSz = erodedSize[n];
		short2 currentImageStart = erodedStart[n];

		// Each thread computes a row
		float sumSOP0 = 0;
		float sumAVG  = 0;

		// Get pointer to actual training image
		float *currentImage = dTIValues[n];

		if (threadIdx.x < currentImageSz.x)
		{
			// We use each thread to compute each row of all training images
			for (unsigned int i = 0; i < currentImageSz.y; ++i)
			{
				// Index0 has coordinates (threadIdx.x, i)
				index0.x = threadIdx.x + currentImageStart.x;
				index0.y = i + currentImageStart.y;

				index1 = index0 + dh_alpha[blockIdx.x];

				sumSOP0 += currentImage[sub2ind(index0,TISize[n])]*
						currentImage[sub2ind(index1,TISize[n])];

				sumAVG  += currentImage[sub2ind(index1,TISize[n])];
			}

			sData[threadIdx.x] += sumSOP0*scaleFactor;
			sData[threadIdx.x+blockDim.x] += sumAVG*scaleFactor;
		}
		else
		{
			sData[threadIdx.x] += 0;
			sData[threadIdx.x+blockDim.x] += 0;
		}
	}

	__syncthreads();

	// Sum over the shared memory registers
	float globalSumSOP0 = 0;
	float globalSumAVG  = 0;
	if (threadIdx.x == 0)
	{
		for (unsigned int i = 0; i < blockDim.x; i++)
		{
			globalSumSOP0 += sData[i];
			globalSumAVG  += sData[i+blockDim.x];
		}

		SOP0[blockIdx.x] = globalSumSOP0;
		AVG[blockIdx.x]  = globalSumAVG;
	}
}

__global__ void computeAVG(float **dTIValues, short2 *dh_alpha,
		short2 *TISize, short2 *erodedSize, short2 *erodedStart,
		float scaleFactor,float *AVG, float *AVGSq, int numImages)

{
	// Shared memory
	extern __shared__ float sData[];
	sData[threadIdx.x] = 0;
	sData[threadIdx.x+blockDim.x] = 0;

	short2 index0;

	for (unsigned int n = 0; n < numImages; n++)
	{
		// Read current training image size and start coordinate
		short2 currentImageSz = erodedSize[n];
		short2 currentImageStart = erodedStart[n];

		// Each thread computes a row
		float sum = 0;
		float sumSq = 0;

		// Get pointer to actual training image
		float *currentImage = dTIValues[n];

		if (threadIdx.x < currentImageSz.x)
		{
			// We use each thread to compute each row of all training images
			for (unsigned int i = 0; i < currentImageSz.y; ++i)
			{
				// Index0 has coordinates (threadIdx.x, i)
				index0.x = threadIdx.x + currentImageStart.x;
				index0.y = i + currentImageStart.y;

				sum  += currentImage[sub2ind(index0,TISize[n])];
				sumSq += currentImage[sub2ind(index0,TISize[n])]*
						currentImage[sub2ind(index0,TISize[n])];
			}

			sData[threadIdx.x] += sum*scaleFactor;
			sData[threadIdx.x+blockDim.x] += sumSq*scaleFactor;
		}
		else
		{
			sData[threadIdx.x] += 0;
			sData[threadIdx.x+blockDim.x] += 0;
		}
	}

	__syncthreads();

	//	// Do a reduction here instead!
	float globalSum = 0;
	float globalSumSq = 0;
	if (threadIdx.x == 0)
	{
		for (unsigned int i = 0; i < blockDim.x; i++)
		{
			globalSum += sData[i];
			globalSumSq += sData[i+blockDim.x];
		}

		AVG[blockIdx.x] = globalSum;
		AVGSq[blockIdx.x] = globalSumSq;
	}
}

GPUKriging::GPUKriging(const TIKrigingParam &parameters, int GPURank)
{
	// Step 0: Copy parameters
	globalParameters = parameters;

	// Step 1: Load Training Image(s)
	numTrainingImages = parameters.TINames.size();

	// Allocate memory for enough pointers to each training image
	hTIVals = new float*[numTrainingImages];
	dTIVals = new float*[numTrainingImages];
	TISizes.clear();
	TISizes.reserve(numTrainingImages);

	// Read training image
	for (unsigned int i = 0; i < numTrainingImages; i++)
	{
		// First copy to temp buffer
		std::string TIFilename = parameters.TIRoot + parameters.TINames[i];
		TIValues.clear();
		loadGrid(TIFilename, TISize, TIValues);

		// Figure out size of training image
		TISizes.push_back(TISize);
		int numPixels = TISize.x * TISize.y;

		// Allocate memory for Training Images
		hTIVals[i] = new float[numPixels];
		std::copy(TIValues.begin(), TIValues.end(),hTIVals[i]);
	}

	// Read hard data
	loadHardData(parameters.hardDataPath,hardDataVals);

	// Initialize the GPU
	initializeGPU(GPURank);


}

GPUKriging::~GPUKriging()
{
	// Free memory
	for (unsigned int i = 0; i < TISizes.size(); i++)
	{
		delete [] hTIVals[i];
		cudaFree(dTIVals[i]);
	}
	delete [] hTIVals;
	delete [] dTIVals;

}

bool GPUKriging::initializeGPU(int GPURank)
{
	// Start GPU
	cudaGetLastError();
	cudaSetDevice(GPURank);
	cudaDeviceProp deviceProp;
	gpuErrchk(cudaGetDeviceProperties(&deviceProp, GPURank));

	std::cout << "Copying Training Image(s) To GPU " << GPURank << " : "
			<< deviceProp.name << std::endl;

	// Allocate an array on GPU (dTIVals to store numTrainingImages worth of
	// pointers

	// float **dTIVals
	gpuErrchk( cudaMalloc(&dTIVals, sizeof(float*)*numTrainingImages));

	// Need a dummy array of pointers that we allocate from the host side then
	// pass to the GPU
	float *deviceTIs[5];

	// Copy training image to GPU
	for (unsigned int i = 0; i < numTrainingImages; i++)
	{
		int numBytes = TISizes[i].x*TISizes[i].y*sizeof(float);

		// First allocate memory on the GPU
		gpuErrchk(cudaMalloc(&deviceTIs[i],numBytes));

		// Copy memory to the GPU
		gpuErrchk(cudaMemcpy(deviceTIs[i],hTIVals[i],numBytes,
				cudaMemcpyHostToDevice));
	}

	// Copy the pointers to the GPU variable
	gpuErrchk(cudaMemcpy(dTIVals, deviceTIs, sizeof(float*)*numTrainingImages, cudaMemcpyHostToDevice));


	// Copy TI Sizes to GPU
	dTISizes = TISizes;

	return true;
}

int GPUKriging::computeKrigingMatrices(int x, int y,
		int numPts2Use)
{
	assert(numPts2Use <= dcf.size());

	// Copy the TIValues to device_vector
	short2 loc = make_short2(x,y);

	// Compute eroded image dimensions into erodedImageStartLoc and
	// erodedImageSize
	erodeImages(loc, numPts2Use);


	float scaleFactor = 1.0f/((erodedSize.x-1)*(erodedSize.y-1));

	float invScaleFactor = 0;
	for (unsigned int k = 0; k < numTrainingImages; k++)
	{
		invScaleFactor += erodedImageSize[k].x*erodedImageSize[k].y;
	}
	invScaleFactor = 1.0f/invScaleFactor;


	// Compute h_data
	h_alpha.clear();
	for (unsigned int k = 0; k < numPts2Use; k++)
		h_alpha.push_back(dcf[k].coord-loc);

	try
	{
		// Copy the DCF and erodedImages info to GPU
		dh_alpha = h_alpha;
		dErodedImageStartLoc = erodedImageStartLoc;
		dErodedImageSize = erodedImageSize;
	}
	catch(thrust::system_error &e)
	{
		// output an error message and exit
		std::cerr << "Error: " << e.what() << std::endl;
		exit(-1);
	}

	short2* dStartLocPtr = thrust::raw_pointer_cast(dErodedImageStartLoc.data());
	short2* dSizePtr 	 = thrust::raw_pointer_cast(dErodedImageSize.data());
	short2* dTISizesPtr	 = thrust::raw_pointer_cast(dTISizes.data());
	// Resize GPU memory if we are using finite Kriging correction
	resizeGPUMatrix(numPts2Use);

	// Define kernel structure
	dim3 gridSizeSOP(numPts2Use,numPts2Use);
	dim3 gridSizeSOP0(numPts2Use,1);
	dim3 gridSizeAvg(1,1);

	// Needs to converted to largest ErodedSize

	dim3 threadSize(erodedSize.x,1,1);
	size_t shardMemSize = erodedSize.x*sizeof(float);

	// Run kernels
	computeSOP<<<gridSizeSOP, threadSize, shardMemSize >>>
			(dTIVals,dh_alphaPtr,dTISizesPtr, dSizePtr, dStartLocPtr,
					invScaleFactor,dSOPMatrixPtr, numTrainingImages);
	cudaDeviceSynchronize();
	gpuErrchk( cudaPeekAtLastError() );
//
//	thrust::host_vector<float> h1 = dSOPMatrix;
//	std::cout << h1[0] << std::endl;
//	std::cout << "Crashed here" << std::endl;
//	exit(0);

	computeSOP0AVG<<<gridSizeSOP0, threadSize, shardMemSize*2>>>
			(dTIVals, dh_alphaPtr, dTISizesPtr, dSizePtr, dStartLocPtr,
					invScaleFactor, dSOP0VectorPtr,dHVectorPtr,numTrainingImages);
	gpuErrchk( cudaPeekAtLastError() );



	computeAVG<<<gridSizeSOP0, threadSize, shardMemSize*2>>>
			(dTIVals, dh_alphaPtr, dTISizesPtr, dSizePtr, dStartLocPtr,
					invScaleFactor,dAvgPtr,dAvgSqPtr,numTrainingImages);
	gpuErrchk( cudaPeekAtLastError() );

	cudaDeviceSynchronize();



	thrust::host_vector<float> C1 = dAvgSq;


	computePenaltyTerms(C1[0]);

	// Free the eroded images information since the next pixel will probably
	// have a different number of points and layout.

	return numPts2Use;
}

void GPUKriging::loadGrid(const std::string &path, short2 &gridSize,
		thrust::host_vector<float> &values)
{
	std::string line;
	std::ifstream myfile;
	myfile.open(path.c_str());
	int numRead = 0;
	int numProperties = 0;

	if (myfile.is_open())
	{
		// First line is name of grid and dimension
		getline (myfile,line);
		unsigned found = line.find_last_of("(");
		line = line.substr(found+1);
		line.erase(line.end()-1, line.end());
		char str[10];
		int null;
		int gridX, gridY;
		std::sscanf (line.c_str(),"%d %[x] %d %[x] %d",
				&gridX, str,&gridY,str, &null);
		gridSize.x = gridX;
		gridSize.y = gridY;

		// Second line is number of properties are in the file
		getline(myfile,line);
		numProperties = std::atoi(line.c_str());

		// Skip lines containing property names
		for (unsigned int i = 0; i < numProperties; i++)
			getline(myfile,line);

		while (myfile.good())
		{
			getline (myfile,line);
			float value = std::atof(line.c_str());
			values.push_back(value);
			numRead ++;
		}
		values.pop_back();
		myfile.close();

		// Make sure we have read proper number of points
		if (values.size() != TISize.x*TISize.y)
			std::cerr << "Misread Grid!" << std::endl;

	}
	else
	{
		std::cerr << "Unable to open file" << std::endl;
	}
}

void GPUKriging::loadHardData(const std::string &path,
		thrust::host_vector<float4> &values)
{
	std::string line;
	std::ifstream myfile;
	myfile.open(path.c_str());
	int numRead = 0;
	int numProperties = 0;
	std::string temp;
	if (myfile.is_open())
	{
		// First line is name of grid and dimension
		getline (myfile,line);

		// Second line is number of properties are in the file
		getline(myfile,line);
		numProperties = std::atoi(line.c_str());

		// Skip lines containing property names
		for (unsigned int i = 0; i < numProperties; i++)
			getline(myfile,line);

		while (!myfile.eof())
		{
			float4 value;

			// Read in files and decrement index by 1 since MATLAB uses
			// 1-indexing and we are using 0-indexing
			myfile >> temp;
			value.x = std::atof(temp.c_str());
			myfile >> temp;
			value.y = std::atof(temp.c_str());
			myfile >> temp;
			value.z = std::atof(temp.c_str());
			myfile >> temp;
			value.w = std::atof(temp.c_str());
			values.push_back(value);
			numRead++;
		}
		values.pop_back();
	}
	else
	{
		std::cerr << "Couldn't open hard data file" << std::endl;
	}

}

void GPUKriging::searchEllipsoid(const int &x, const int &y)
{
	dcf.clear();
	datConfig dataPoint;

	std::vector<datConfig> allPoints;

	for (unsigned int i = 0; i < hardDataVals.size(); i++)
	{

		float rhs = (x - hardDataVals[i].x)*(x - hardDataVals[i].x)/
				(globalParameters.searchEllipseDim.x*
						globalParameters.searchEllipseDim.x) +
						(y - hardDataVals[i].y)*(y - hardDataVals[i].y)/
						(globalParameters.searchEllipseDim.y*
								globalParameters.searchEllipseDim.y);

		// Write coordinate of point
		dataPoint.coord.x = hardDataVals[i].x;
		dataPoint.coord.y = hardDataVals[i].y;
		dataPoint.dist = (x - hardDataVals[i].x)*
				(x - hardDataVals[i].x) + (y - hardDataVals[i].y)*
				(y - hardDataVals[i].y);
		dataPoint.dist = sqrt(dataPoint.dist);
		dataPoint.val = hardDataVals[i].w;
		dataPoint.ellipRad = rhs;

		// Push onto vector
		allPoints.push_back(dataPoint);
	}

	// Sort allPoints vector based off distance
	std::sort(allPoints.begin(),allPoints.end());

	// Decide which points we want to keep in the DCF
	int minPoints = globalParameters.minConditionNum;

	// First copy all minPoints to DCF
	for (unsigned int i = 0; i < minPoints; i++)
		dcf.push_back(allPoints[i]);

	// Copy regions that are less than 100;
	for (unsigned int i = minPoints; i < allPoints.size(); i++)
	{
		if (allPoints[i].ellipRad < 1.0f)
			dcf.push_back(allPoints[i]);

		if (dcf.size() >= globalParameters.maxConditionNum)
			break;
	}


	// Compute data configuration for current pixel
	numDCFPts = dcf.size();
}

void GPUKriging::erodeImages(const short2 &ptLoc, int numClosestPts)
{
	std::pair<short2, short2> erodedCoords;
	datConfig temp;
	short2 maxCoordinates, minCoordinates;
	int x1,x2,y1,y2;

	// Sort the DCF according to distance from current point
	temp = *std::max_element(dcf.begin(),dcf.begin()+numClosestPts,lessThanX);
	maxCoordinates.x = temp.coord.x;
	temp = *std::max_element(dcf.begin(),dcf.begin()+numClosestPts,lessThanY);
	maxCoordinates.y = temp.coord.y;
	temp = *std::min_element(dcf.begin(),dcf.begin()+numClosestPts,lessThanX);
	minCoordinates.x = temp.coord.x;
	temp = *std::min_element(dcf.begin(),dcf.begin()+numClosestPts,lessThanY);
	minCoordinates.y = temp.coord.y;

	// Compute the size of the kernel
	x1 = abs(min(ptLoc.x, minCoordinates.x) - ptLoc.x);
	x2 = abs(max(ptLoc.x, maxCoordinates.x) - ptLoc.x);
	y1 = abs(min(ptLoc.y, minCoordinates.y) - ptLoc.y);
	y2 = abs(max(ptLoc.y, maxCoordinates.y) - ptLoc.y);


	// Clear vectors
	erodedImageSize.clear();
	erodedImageStartLoc.clear();

	for (unsigned int i = 0; i < numTrainingImages; i++)
	{
		// Compute the top left and bottom right coordinates of erodedImage
		short2 xCoords, yCoords;
		xCoords.x = x1;
		xCoords.y = TISizes[i].x - x2 - 1;
		yCoords.x = y1 ;
		yCoords.y = TISizes[i].y - y2 - 1;

		erodedCoords = std::make_pair(xCoords,yCoords);

		erodedSize.x = erodedCoords.first.y - erodedCoords.first.x;
		erodedSize.y = erodedCoords.second.y - erodedCoords.second.x;
		erodedStart.x = erodedCoords.first.x;
		erodedStart.y = erodedCoords.second.x;
		erodedSize.x = max(erodedSize.x, 2);
		erodedSize.y = max(erodedSize.y, 2);

		erodedImageSize.push_back(erodedSize);
		erodedImageStartLoc.push_back(erodedStart);
	}
}

void GPUKriging::resizeGPUMatrix(int size)
{
	// Copy h_data to GPU
	dh_alphaPtr = thrust::raw_pointer_cast(dh_alpha.data());


	// Resize output matrices
	dSOPMatrix.resize(size*size);
	dSOP0Vector.resize(size);
	dHVector.resize(size);
	dAvg.resize(1);
	dAvgSq.resize(1);

	dSOPMatrixPtr  = thrust::raw_pointer_cast(dSOPMatrix.data());
	dSOP0VectorPtr = thrust::raw_pointer_cast(dSOP0Vector.data());
	dHVectorPtr    = thrust::raw_pointer_cast(dHVector.data());
	dAvgPtr 	   = thrust::raw_pointer_cast(dAvg.data());
	dAvgSqPtr	   = thrust::raw_pointer_cast(dAvgSq.data());

	dErodedImageSzPtr = thrust::raw_pointer_cast(dErodedImageSize.data());
	dErodedImageStartPtr = thrust::raw_pointer_cast(dErodedImageStartLoc.data());
}

void GPUKriging::computePenaltyTerms(const float &C0)
{
	// Compute penalty terms
	penaltyValues.clear();
	for (unsigned int i = 0; i < numDCFPts; ++i)
	{
		float hNorm = h_alpha[i].x*h_alpha[i].x + h_alpha[i].y*h_alpha[i].y;
		hNorm = sqrt(hNorm);

		if (hNorm > globalParameters.r)
		{
			float penaltyTerm = (hNorm - globalParameters.r)/
					(globalParameters.R - hNorm);
			penaltyTerm *= penaltyTerm;
			penaltyTerm *= C0;

			penaltyValues.push_back(4*penaltyTerm);
		}
		else
			penaltyValues.push_back(0);
	}
}
