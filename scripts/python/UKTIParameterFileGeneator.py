#!/usr/bin/env python
import sys
import getopt
import string
import math
import os
from lxml import etree

import sys
import subprocess

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg


class UKTIParamterFileGenerator:
	def __init__(self, par):
		self.par = par

	def writefile(self, outputPath):
		root = etree.Element("parameters")
		root.clear()

		# Get UKTI RootDir
		UKTIDirChild = etree.SubElement(root, "UKTIRootDir")
		UKTIDirChild.set("name", "UKTIRootDir")
		UKTIDirChild.text = self.par['UKTIRootDir']

		DataDir = self.par['UKTIRootDir'] + '/' + self.par['DataDir']
		ResultsDir = self.par['UKTIRootDir'] + '/' + self.par['ResultsDir']
		ResultsPath = ResultsDir + self.par['DataSetName'] + '/' + str(self.par['NumHardData']) + '/' + self.par["HDataName"] + '_result'

		# Write DataSetName
		DataSetChild = etree.SubElement(root, "DataSetName")
		DataSetChild.set("name","DataSetName")
		DataSetChild.text = self.par['DataSetName']

		# Output path
		ResultsPathChild = etree.SubElement(root,"ResultsPath")
		ResultsPathChild.set("name","ResultsPath")
		ResultsPathChild.text = ResultsPath

		# Write path to where the TIs are stored
		TIRootChild = etree.SubElement(root, "TIRoot")
		TIRootChild.set("name", "TIRoot")
		TIRootChild.text = DataDir + self.par['TIDir'] + self.par['DataSetName'] + '/'

		# Write TI Name
		TIName = etree.SubElement(TIRootChild, 'TIName0')
		TIName.text = 'TI'+str(self.par['TINum'])

		# Write path where hard data is stored
		hardDataPathChild = etree.SubElement(root, "hardDataPath")
		hardDataPathChild.set("name","hardDataPath")
		hardDataPathChild.text = DataDir + self.par['HDDir'] + self.par['DataSetName'] + '/' + self.par['HDataName']
		hardDataNumChild = etree.SubElement(root, "hardDataPoints")
		hardDataNumChild.set("name","hardDataPoints")
		hardDataNumChild.text = str(self.par['NumHardData'])

		# Write search ellipse dimension
		searchEllipseXChild = etree.SubElement(root, "searchEllipseX")
		searchEllipseXChild.set("name","searchEllipseX")
		searchEllipseXChild.text =  str(self.par['SearchEllipseDim'])
		searchEllipseYChild = etree.SubElement(root, "searchEllipseY")
		searchEllipseYChild.set("name","searchEllipseY")
		searchEllipseYChild.text =  str(self.par['SearchEllipseDim'])
		searchEllipseZChild = etree.SubElement(root, "searchEllipseZ")
		searchEllipseZChild.set("name","searchEllipseZ")
		searchEllipseZChild.text =  '1'

		# Write Min/Max conditioning number
		minCondNumChild = etree.SubElement(root, "minConditionNum")
		minCondNumChild.set("name","minCondNum")
		minCondNumChild.text = str(self.par['MinCondNum'])
		maxCondNumChild = etree.SubElement(root, "maxConditionNum")
		maxCondNumChild.set("name","maxCondNum")
		maxCondNumChild.text = str(self.par['MaxCondNum'])

        # Write grid size
		gridSizeXChild = etree.SubElement(root, "gridSizeX")
		gridSizeXChild.set("name","gridSizeX")
		gridSizeXChild.text = str(self.par['GridSize'][0])
		gridSizeYChild = etree.SubElement(root, "gridSizeY")
		gridSizeYChild.set("name","gridSizeY")
		gridSizeYChild.text = str(self.par['GridSize'][1])

		# Flag to turn on Finite Kriging Correction
		finiteFlagChild = etree.SubElement(root, "UseFiniteKrig")
		finiteFlagChild.set("name", "UseFiniteKrig")
		finiteFlagChild.text = str(int(self.par['UseFiniteKrig']))

		# Write multiple TI flag
		useMultipleTIChild = etree.SubElement(root, "UseMultipleTI")
		useMultipleTIChild.set("name", "UseMultipleTI")
		useMultipleTIChild.text = str(int(self.par['UseMultipleTI']))

		# Write penalty terms
		penaltyFlagChild = etree.SubElement(root, "usePenalty")
		penaltyFlagChild.set("name", "usePenalty")
		penaltyFlagChild.text = str(int(self.par['UsePenaltyCorrection']))

        # Write R and r terms
		rFlag = etree.SubElement(root, "r")
		rFlag.set("name", "r")
		rFlag.text = str(int(self.par['rPenalty']))
		RFlag = etree.SubElement(root, "R")
		RFlag.set("name", "R")
		RFlag.text = str(int(self.par['RPenalty']))

		with open(outputPath, "w") as text_file:
			text_file.write(etree.tostring(root, pretty_print=True))

	def createSlurmfile(self,name,email,nodes,ppn,time,content,nodetype='defq'):
	    #This function creates the batch job for the Slurm scheduler

		pbs_dirt = 'pbs'

		# Allow for sub hour timing
		hours = int(math.floor(time))
		minPercent = time-hours
		minutes = int(round(minPercent*60))

		lines = []
		lines.append('#!/bin/bash')
		lines.append('#SBATCH -J %s' % (name))
		lines.append('#SBATCH --nodes=%d' % (nodes))
		lines.append('#SBATCH --cpus-per-task=%d '% (ppn))
		lines.append('#SBATCH --time=%d:%02d:00' % (hours,minutes))
		lines.append('#SBATCH --job-name=%s' % name)
		lines.append('#SBATCH -o %s/%s.log' % (pbs_dirt,name))
		lines.append('#SBATCH --ntasks-per-node=1') # 1 task per node, then it gives all cpus to task (OMP) 
		if nodetype:
		  lines.append('#SBATCH --partition=%s'%nodetype)
		if email:
		    lines.append('#SBATCH --mail-type=FAIL')
		    lines.append('#SBATCH --mail-user=%s' % email)
		lines.append('#-----------')
		lines.append('cd %s' % os.environ.get('UKTIPATH'))
		lines.append('srun -n %d bin/UKTI %s' %(ppn, content))

		file = open('%s/%s' % (pbs_dirt,name),'w')
		text = string.join(lines,'\n')
		file.write(text)
		file.close()
         
def main(arg=None):

    par = dict( DataDir = 'data/',
    			ResultsDir = 'results/',
    			TIDir = 'trainingImages/',
    			HDDir = 'hardData/',
    			UKTIRootDir = os.environ.get('UKTIPATH'),
    			HDataName = 'hData0',
    			DataSetName = 'WalkerLake',
    			GridSize = [260, 300],
    			TIName = 'TI',
    			SearchEllipseDim = 25,
    			MinCondNum = 20,
    			MaxCondNum = 50,
    			NumHardData = 400,
    			UseFiniteKrig = True,
    			UseMultipleTI = False,
    			UsePenaltyCorrection = False,
    			HardDataName = 'hData',
    			rPenalty = 5,
    			RPenalty = 20,
    			TINum  = 1,
        )

    Gen = UKTIParamterFileGenerator(par)
    Gen.writefile('TrialTest.xml')
    Gen.createSlurmfile('TestJobName','lewisli@stanford.edu',1,1,2,'pbs/TrialText.xml','defq')

if __name__ == "__main__":
    main(sys.argv[1:])