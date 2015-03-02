/*
 * KrigingSolver.cpp
 *
 *  Created on: Jul 17, 2013
 *      Author: lewisli
 */

#include "KrigingSolver.h"
#include "util/rapidxml.hpp"
#include "util/rapidxml_utils.hpp"
#include "util/helper_cuda.h"

#include <assert.h>
#include <iostream>
#include <string>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef USEOPENMP
#include <omp.h>
#endif
struct stat sb;

using namespace rapidxml;
using namespace std;

std::string stoi(int a)
{
    std::stringstream ss;
    ss << a;
    return ss.str();
}
KrigingSolver::KrigingSolver(char* parameterPath)
{
    // Read parameterXML file
    processParameters(parameterPath);

    // Generate results path
    filename = generateOutputFilename(parameters);

    // Query system for GPU information
    if (getSystemInformation() == 0)
    {
        std::cerr << "Can't find any CUDA-Enabled GPUs on this machine."
                     " Perhaps try using MATLAB version" << std::endl;
        throw ("No CUDA-enabled GPUs");
    }

    // Allocate Memory for output file
    finalEstimate = new float[gridSizeX*gridSizeY];
    finalVariance = new float[gridSizeX*gridSizeY];
    finalConditionNumber = new int[gridSizeX*gridSizeY];

    // Start timer
    t0 = get_timestamp();
    solve();
}

void KrigingSolver::processParameters(char *parameterPath)
{
    // Open XML file for reading format
    rapidxml::file<> xmlFile(parameterPath);
    rapidxml::xml_document<> doc;
    doc.parse<0>(xmlFile.data());
    xml_node<> *node, *childrenNode;


    node = doc.first_node()->first_node("TIRoot");
    parameters.TIRoot = std::string(node->value());
    parameters.TINames.clear();


    // Assume that there is at least one Training image specified
    childrenNode = node->first_node("TIName0");

    if (childrenNode == 0)
    {
        std::cerr << "Missing Training Image Name" << std::endl;
        exit(-1);
    }

    while(childrenNode != 0)
    {
        parameters.TINames.push_back(std::string(childrenNode->value()));
        childrenNode = childrenNode->next_sibling();
    }

    node = doc.first_node()->first_node("hardDataPath");
    parameters.hardDataPath = std::string(node->value());
    node = doc.first_node()->first_node("searchEllipseX");
    parameters.searchEllipseDim.x = std::atoi(node->value());
    node = doc.first_node()->first_node("searchEllipseY");
    parameters.searchEllipseDim.y = std::atoi(node->value());
    node = doc.first_node()->first_node("searchEllipseZ");
    parameters.searchEllipseDim.z = std::atoi(node->value());
    node = doc.first_node()->first_node("maxConditionNum");
    parameters.maxConditionNum = std::atoi(node->value());
    node = doc.first_node()->first_node("minConditionNum");
    parameters.minConditionNum = std::atoi(node->value());
    node = doc.first_node()->first_node("gridSizeX");
    gridSizeX = std::atoi(node->value());
    node = doc.first_node()->first_node("gridSizeY");
    gridSizeY = std::atoi(node->value());
    node = doc.first_node()->first_node("r");
    parameters.r = std::atoi(node->value());
    node = doc.first_node()->first_node("R");
    parameters.R =  std::atoi(node->value());
    node = doc.first_node()->first_node("usePenalty");
    parameters.usePenalty = bool(std::atoi(node->value()));
    node = doc.first_node()->first_node("UseFiniteKrig");
    parameters.useFinite = bool(std::atoi(node->value()));
    node = doc.first_node()->first_node("hardDataPoints");
    parameters.numHardPts = std::atoi(node->value());
    node = doc.first_node()->first_node("UseMultipleTI");
    parameters.useMultipleTI = bool(std::atoi(node->value()));
    node = doc.first_node()->first_node("ResultsPath");
    parameters.resultsPath = std::string(node->value());

}

bool KrigingSolver::createFolder(std::string path)
{
    // Check if path exists, if not create it
    if (stat(path.c_str(), &sb) == -1)
    {
        size_t pre=0,pos;
        std::string dir;
        int mdret;

        if(path[path.size()-1]!='/'){
            // force trailing / so we can handle everything in loop
            path+='/';
        }

        while((pos=path.find_first_of('/',pre))!=std::string::npos){
            dir=path.substr(0,pos++);
            pre=pos;
            if(dir.size()==0) continue; // if leading / first time is 0 length
            if((mdret=mkdir(dir.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH))
                    && errno!=EEXIST){
                return true;
            }
        }
        return true;
    }

    // if folder already exists
    return true;
}
std::string KrigingSolver::generateOutputFilename(const TIKrigingParam
                                                  &parameters)
{
    unsigned found = parameters.resultsPath.find_last_of("/\\");
    std::string resultsDir = parameters.resultsPath.substr(0,found);

    if (!createFolder(resultsDir))
          exit(-1);

    return parameters.resultsPath;
}

KrigingSolver::~KrigingSolver()
{
    // TODO Auto-generated destructor stub
    delete [] finalEstimate;
    delete [] finalVariance;
    delete [] finalConditionNumber;
}

void KrigingSolver::computeWeights(int numPts, VectorXf &w,
                                   MatrixXf &A, VectorXf &b, bool usePenalty,
                                   GPUKriging *gpuKrig)
{
    assert(numPts <= gpuKrig->numDCFPts);

    int systemSize = numPts + 1;
    int SOPSize = numPts;

    A.resize(systemSize,systemSize);
    b.resize(systemSize);
    w.resize(b.rows());

    thrust::host_vector<float> h1 = gpuKrig->dSOPMatrix;
    thrust::host_vector<float> h2 = gpuKrig->dSOP0Vector;
    thrust::host_vector<float> h3 = gpuKrig->dHVector;
    thrust::host_vector<float> h4 = gpuKrig->dAvg;

    for (unsigned int ii = 0; ii < SOPSize; ii++)
        for (unsigned int jj = 0; jj < SOPSize; jj++)
            A(ii,jj) = h1[ii*SOPSize+jj];

    for (unsigned int ii = 0; ii < SOPSize; ii++)
    {
        b(ii) = h2[ii];
        A(SOPSize,ii) = h3[ii];
        A(ii,SOPSize) = h3[ii];

        // Add penalty term
        if (usePenalty) A(ii,ii) += gpuKrig->penaltyValues[ii];
    }

    A(SOPSize,SOPSize) = 0;
    b(SOPSize) = h4[0];

    // Decomposition to smooth out covariacnes
    //	JacobiSVD<MatrixXf> svd(A, ComputeFullU | ComputeFullV);
    //	MatrixXf r = svd.singularValues();
    //
    //	for (int rr = r.rows()-1; rr < r.rows(); rr++)
    //		r(rr,0) = 0;
    //	A = svd.matrixU()*r.asDiagonal()*  svd.matrixV().transpose();
    //	A(SOPSize,SOPSize) = 0;

    w = A.fullPivLu().solve(b);
}

int KrigingSolver::getSystemInformation()
{
    // determine the number of CUDA capable GPUs
    int num_gpus; cudaGetDeviceCount(&num_gpus);
    if(num_gpus < 1)
    {
        printf("no CUDA capable devices were detected\n");
        return 0;
    }

    /////////////////////////////////////////////////////////////////
    // display CPU and GPU configuration
    //
#ifdef USEOPENMP
    printf("Number of host CPUs:\t%d\n", omp_get_num_procs());
#endif
    printf("Number of CUDA devices:\t%d\n", num_gpus);

    for(int i = 0; i < num_gpus; i++)
    {
        cudaDeviceProp dprop;
        cudaGetDeviceProperties(&dprop, i);
        printf("  GPU# %d: %s\n", i, dprop.name);
    }

    printf("---------------------------\n");

    return num_gpus;
}
void KrigingSolver::solveSystem(int i, int j, GPUKriging *gpuKrig)
{
    // Declare Eigen objects
    VectorXf b,lambda, tempLambda, fullLambda;
    MatrixXf A;

    // Perform search ellipse
    gpuKrig->searchEllipsoid(i,j);

    int numCondData = gpuKrig->numDCFPts;
    int index = i*gridSizeY + j;

    bool useAverage = gpuKrig->globalParameters.useFinite;
    bool usePenalty = gpuKrig->globalParameters.usePenalty;

    // Make sure we have conditioning points in the region
    if (numCondData > 0)
    {
        //unsigned int cpu_thread_id = omp_get_thread_num();

        // Without finite Kriging correction
        if (!useAverage)
        {
            // Compute complete Kriging system
            gpuKrig->computeKrigingMatrices(i,j,numCondData);
            computeWeights(gpuKrig->numDCFPts, lambda, A, b,
                           usePenalty,gpuKrig);
        }
        else
        {
            lambda.resize(numCondData+1);
            fullLambda.resize(numCondData+1);
            lambda.setZero();

            int kStart = 0;
            if (numCondData >=2) kStart = 1;

            for (unsigned int k = kStart; k < numCondData; k++)
            {
                fullLambda.setZero();
                gpuKrig->computeKrigingMatrices(i,j,k+1);
                computeWeights(k+1, tempLambda, A,b, usePenalty,gpuKrig);

                for (unsigned int kk = 0; kk < k+1; kk++)
                    fullLambda(kk) = tempLambda(kk);

                // Write Lagragian multiplier
                fullLambda(numCondData) = tempLambda(k+1);

                lambda = lambda + fullLambda;
            }
            // Average the weights
            lambda = lambda/numCondData;
        }

        float result = computeEstimate(lambda,gpuKrig);
        float variance = computeVariance(lambda,b,gpuKrig);

        finalEstimate[index] = result;
        finalVariance[index] = variance;
        finalConditionNumber[index] = numCondData;
    }
    else
    {
        finalEstimate[index] = 0;
        finalVariance[index] = 0;
        finalConditionNumber[index] = 0;
    }
}

void KrigingSolver::writeToFile()
{
    std::ofstream myfile;
    myfile.open (filename.c_str());

    if (!myfile.is_open())
    {
        std::cerr << "Error Creating Output File" << std::endl;
        return;
    }

    for (unsigned int i = 0; i < gridSizeX; i++)
        for (unsigned int j = 0; j < gridSizeY; j++)
        {
            int index = i*gridSizeY + j;
            myfile << i << " " << j << " " << finalEstimate[index] << " "
                   << finalVariance[index] << " "
                   << finalConditionNumber[index] << std::endl;
        }

    myfile.close();

}
void KrigingSolver::solve()
{
    // Number of GPU cores on machine
    int numGPUs; cudaGetDeviceCount(&numGPUs);

    // For now, use 1 CPU thread for each GPU
    std::cout << "Creating " << numGPUs << " Threads" << std::endl;

#ifdef USEOPENMP
    // Number of CPU cores on machine
    int numHostThreads = omp_get_num_procs();
    omp_set_num_threads(numGPUs);
#endif




    // Set up parallel threads
    int numComplete = 0;
#pragma omp parallel
    {
#ifdef USEOPENMP
        unsigned int cpu_thread_id = omp_get_thread_num();
        unsigned int num_cpu_threads = omp_get_num_threads();
#else
        unsigned int cpu_thread_id = 0;
        unsigned int num_cpu_threads =1;
#endif


        // set and check the CUDA device for this CPU thread
        int gpu_id = -1;
        cudaSetDevice(cpu_thread_id % numGPUs);
        cudaGetDevice(&gpu_id);

        printf("CPU thread %d (of %d) uses CUDA device %d\n",
               cpu_thread_id+1, num_cpu_threads, gpu_id);

        GPUKriging *gpuKrig = new GPUKriging(parameters,gpu_id);

#pragma omp for
        for (int index = 0; index < gridSizeX*gridSizeY; index++)
        {
            int i = index/gridSizeY;
            int j = index % gridSizeY;
            solveSystem(i,j,gpuKrig);

#pragma omp atomic
            numComplete++;

            if (numComplete % (gridSizeX*3) == 0)
            {
                double secs;
                timestamp_t t1 = get_timestamp();
                secs = (t1 - t0) / 1000000.0L;
                estimateTimeRemaining(float(numComplete*100)/
                                      (gridSizeX*gridSizeY),secs);
            }
        }
    }

    std::cout << "Writing Output To: " << filename << std::endl;
    writeToFile();
    std::cout << "Done" << std::endl;
}

void KrigingSolver::estimateTimeRemaining(float completePercentage,
                                          double timeElapsed)
{
    double timeEstimate = timeElapsed * (1.0/(completePercentage/100)) -
            timeElapsed;

    int time,hour,min,sec;
    time = int(timeElapsed);
    hour=time/3600;
    time=time%3600;
    min=time/60;
    time=time%60;
    sec=time;

    std::cout << std::setw(8) << completePercentage << " %. Running for "
              << std::setw(8) << hour << " h, " << min << " m, " << sec << " s.";

    time = int(timeEstimate);
    hour=time/3600;
    time=time%3600;
    min=time/60;
    time=time%60;
    sec=time;
    std::cout << " Estimated Remaining: " <<
                 hour << " h, " << min << " m, " << sec << " s." << std::endl;
}

float KrigingSolver::computeEstimate(const VectorXf &w,
                                     GPUKriging *gpuKrig)
{
    float estimate = 0;

    for (unsigned int i = 0; i < w.rows() - 1; i++)
        estimate += w(i)*gpuKrig->dcf[i].val;

    return estimate;
}

float KrigingSolver::computeVariance(const VectorXf &w, const VectorXf &d,
                                     GPUKriging *gpuKrig)
{
    thrust::host_vector<float> SOP00 = gpuKrig->dAvgSq;
    return SOP00[0] - w.dot(d);
}
