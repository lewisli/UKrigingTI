/* *
 * Copyright 1993-2012 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 */
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "KrigingSolver.h"


int main(int argc, char *argv[])
{
	// Ensure that filename is provided
	if (argc != 2)
	{
		std::cerr << "Invalid parameters" << std::endl;
		return -1;
	}

	KrigingSolver *a = new KrigingSolver(argv[1]);

	return 0;
}

