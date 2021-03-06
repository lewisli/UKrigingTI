/*
 * MPSKriging.h
 *
 *  Created on: Jul 15, 2013
 *      Author: lewisli
 */

#ifndef MPSKRIGING_H_
#define MPSKRIGING_H_

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/reduce.h>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <assert.h>


typedef struct
{
	std::string TIRoot;
	std::vector<std::string> TINames;
	std::string hardDataPath;
	int numHardPts;
	float3 searchEllipseDim;
	int maxConditionNum;
	int minConditionNum;
	float R;
	float r;
	bool usePenalty;
	bool useFinite;
	bool useMultipleTI;
    std::string resultsPath;
}TIKrigingParam;

struct datConfig
{
	// The coordinates of the hard data location
	short2 coord;

	// The actual value of the hard data point
	float val;

	// The distance to the point we are trying to estimate
	float dist;

	// The RHS to the ellipse equation
	float ellipRad;

	// Operator used for sorting vector by ascending
	bool operator < (const datConfig &a) const
	{
		return (dist < a.dist);
	}
};

class GPUKriging
{
public:
	GPUKriging(const TIKrigingParam &parameters, int GPURank);
	~GPUKriging();

	// NVCC has issues compiling Eigen, hencce we need this friend class
	friend class KrigingSolver;

protected:
	int computeKrigingMatrices(int x, int y, int numPts2Use, int TINum);
	int computeKrigingMatrices(int x, int y, int numPts2Use);
	void getGPUInformation();

private:
	void loadGrid(const std::string &path, short2 &gridSize,
			thrust::host_vector<float> &values);

	void loadHardData(const std::string &path,
			thrust::host_vector<float4> &values);

	void searchEllipsoid(const int &x, const int &y);

	void erodeImages(const short2 &ptLoc, int numClosestPts);

	void resizeGPUMatrix(int size);

	void computePenaltyTerms(const float &C0);

	bool initializeGPU(int GPURank);

	// This needs to be converted into a set of pointers to pointers
	float **hTIVals, **dTIVals;


	thrust::host_vector<short2>   erodedImageSize, 		erodedImageStartLoc;
	thrust::device_vector<short2> dErodedImageStartLoc, dErodedImageSize;
	thrust::device_vector<short2> dTISizes;
	thrust::host_vector<short2> TISizes;

	thrust::host_vector<float> TIValues;
	thrust::device_vector<short2> dh_alpha;
	thrust::host_vector<float4> hardDataVals;
	thrust::host_vector<short2> h_alpha;
	std::vector<datConfig> dcf;

	std::vector<float> penaltyValues;

	thrust::device_vector<float> dSOPMatrix;
	thrust::device_vector<float> dSOP0Vector;
	thrust::device_vector<float> dHVector;
	thrust::device_vector<float> dAvg;
	thrust::device_vector<float> dAvgSq;

	float *dSOPMatrixPtr, *dSOP0VectorPtr, *dHVectorPtr, *dAvgPtr, *dAvgSqPtr;
	short2 *dErodedImageSzPtr, *dErodedImageStartPtr;

	short2 TISize;
	TIKrigingParam globalParameters;

	float * dTIValuesPtr;
	short2 *dh_alphaPtr;

	int numDCFPts;
	int numTrainingImages;

	short2 erodedSize,erodedStart;

};

#endif /* MPSKRIGING_H_ */
