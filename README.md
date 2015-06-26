UKrigingTI
==========
This is the source code for **Li, L., Romary, T., & Caers, J. (2015). *Universal kriging with training images. Spatial Statistics.*** It provides a method to perform kriging without variograms when a representative training image is available.

Dependencies
============
1. CUDA V5.5+
2. Python 2.7.x 
3. CMake 2.8+

Build
=====
Should be able to generate Makefile with CMake. Tested on Ubuntu 14.04 with CUDA V5.5, NVIDIA GTX670 Graphics card, and on Stanford ICME GPU Cluster (NVIDIA C2070)

Usage
=====
Parameters for the program are passed using an XML format.

    ./UKTI parameters.xml 

Refer to [here](../scripts/python/TrialTest.xml) for an example of a parameter xml file. Furthermore, this [python script](../scripts/python/UKTIParameterFileGeneator.py) can be used to generate parameter files.
