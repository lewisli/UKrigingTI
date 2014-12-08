#!/usr/bin/python
import lxml.etree
import lxml.builder
import subprocess
import sys


E = lxml.builder.ElementMaker()
ROOT = E.parameters
FIELD1 = E.TIPath
FIELD2 = E.hardDataPath
FIELD2A = E.hardDataPoints
FIELD3 = E.searchEllipseDimX
FIELD4 = E.searchEllipseDimY
FIELD5 = E.searchEllipseDimZ
FIELD6 = E.maxConditionNum
FIELD7 = E.gridSizeX
FIELD8 = E.gridSizeY
FIELD9 = E.r
FIELD10 = E.usePenalty
FIELD13 = E.useFinite
FIELD11 = E.R
FIELD12 = E.minConditionNum

TIRoot = ['../data/DS-NonStationary/']
trialNumber = [3,4]

HardData = [TIRoot[0] + 'hardData']
HardDatPts = ['25','50','100','200']

searchEllipseDimX = ['50']
maxConditionNumber = ['15']
minConditionNumber = ['1','12']
usePenaltyVal = ['0','1'];
useFinite = ['0','1']
gridSize = ['100','200']

TIPath = []
for a in range(0,len(trialNumber)):
    TIPath.append(TIRoot[0] + 'TI' + str(trialNumber[a]))


for a in range(0,len(TIPath)):
    for i in range(0,len(minConditionNumber)):
        for j in range(0,len(HardDatPts)):
            for k in range(0,len(searchEllipseDimX)):
                for l in range(0,len(usePenaltyVal)):
                    for m in range(0,len(useFinite)):
                        the_doc = ROOT(
                               
                            FIELD1(TIPath[a], name=''),
                            FIELD2( (HardData[0]+HardDatPts[j]), name=''),
                            FIELD2A(HardDatPts[j], name = ''),
                            FIELD3(searchEllipseDimX[k], name=''),
                            FIELD4(searchEllipseDimX[k], name=''),
                            FIELD5('1', name=''),
                            FIELD6(maxConditionNumber[0], name=''),
                            FIELD12(minConditionNumber[i], name=''),
                            FIELD13(useFinite[m],name = ''),
                            FIELD7(gridSize[0], name=''),
                            FIELD8(gridSize[1], name=''),
                            FIELD9('5', name= ''),
                            FIELD10(usePenaltyVal[l],name = ''),
                            FIELD11('20',name = '')
                            
                            )   
                        print lxml.etree.tostring(the_doc, pretty_print=True)
                        
                        outputFileName = ['trial' + str(trialNumber[a]) + '.xml']
                        with open(outputFileName[0], "w") as text_file:
                            text_file.write(lxml.etree.tostring(the_doc, pretty_print=True))
                            
                        subprocess.call(["../Release/MPSKriging", outputFileName[0]])

