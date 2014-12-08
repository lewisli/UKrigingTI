/*
 * KrigingSolver.h
 *
 *  Created on: Jul 17, 2013
 *      Author: lewisli
 */

#ifndef KRIGINGSOLVER_H_
#define KRIGINGSOLVER_H_

class GPUKriging;

#include "Eigen/Dense"
#include <sys/time.h>
#include <iomanip>
#include "GPUKriging.h"

using namespace Eigen;

typedef unsigned long long timestamp_t;

static timestamp_t
get_timestamp ()
{
	struct timeval now;
	gettimeofday (&now, NULL);
	return  now.tv_usec + (timestamp_t)now.tv_sec * 1000000;
}


class KrigingSolver {
public:
	KrigingSolver(char* parameterPath);
	virtual ~KrigingSolver();

private:
	void solve();
	int getSystemInformation();
	float computeEstimate(const VectorXf &w, GPUKriging *gpuKrig);
	float computeVariance(const VectorXf &w, const VectorXf &d, GPUKriging *gpuKrig);
	void computeWeights(int numPts, VectorXf &w, MatrixXf &A,
			VectorXf &b, bool usePenalty,GPUKriging *gpuKrig);
	void solveSystem(int i, int j, GPUKriging *gpuKrig);
	void processParameters(char *parameterPath);
	std::string generateOutputFilename(const TIKrigingParam &parameters);
	void writeToFile();
	void estimateTimeRemaining(float completePercentage, double timeElapsed);

	bool createFolder(std::string path);

	timestamp_t t0;
	std::string filename;
	int gridSizeX, gridSizeY;
	float *finalEstimate, *finalVariance;
	int *finalConditionNumber;
	TIKrigingParam parameters;
};

#endif /* KRIGINGSOLVER_H_ */
