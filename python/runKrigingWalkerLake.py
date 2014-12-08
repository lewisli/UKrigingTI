from lxml import etree
import sys
import subprocess

TIRoot = '../data/WLake/'
trialNumber = [1,2,3]
HardData = [TIRoot + 'hardData']
HardDatPts = ['25','50','100','200']
searchEllipseDimX = ['25','50','100']
maxConditionNumber = ['50']
minConditionNumber = ['1','12', '20']
usePenaltyVal = ['0'];
useFinite = ['0','1']
gridSize = ['260', '300']
useMultipleTI = 0;
R = 20
r = 5

root = etree.Element("parameters")
    
# Generate TI Paths
TIPath = []
if useMultipleTI == 1:
    # Use multiple training images per run
    trialNumberLength = 1;
    imagePerRun = len(trialNumber)
else:
    # Using multiple training images
    trialNumberLength = len(trialNumber)
    imagePerRun = 1

for a in range(0,trialNumberLength):
    for b in range(0, len(HardDatPts)):
        for c in range(0, len(searchEllipseDimX)):
            for d in range(0, len(minConditionNumber)):
                for e in range(0, len(useFinite)):
                    for f in range(0, len(usePenaltyVal)):
                        # Create new document
                        root.clear()
                    
                        # Write TI Root 
                        TIRootChild = etree.SubElement(root, "TIRoot")
                        TIRootChild.set("name", "TIRoot")
                        TIRootChild.text = TIRoot
                        
                        for z in range(0, imagePerRun):
                            child = etree.SubElement(TIRootChild, 'TIName'+str(z))
                            child.text = 'TI'+str(trialNumber[a+z])
                            
                        # Write Hard Data Path
                        hardDataPathChild = etree.SubElement(root, "hardDataPath")
                        hardDataPathChild.set("name","hardDataPath")
                        hardDataPathChild.text = HardData[0] + HardDatPts[b]
                        hardDataNumChild = etree.SubElement(root, "hardDataPoints")
                        hardDataNumChild.set("name","hardDataPoints")
                        hardDataNumChild.text = HardDatPts[b]
                        
                        # Write search ellipse
                        searchEllipseXChild = etree.SubElement(root, "searchEllipseX")
                        searchEllipseXChild.set("name","searchEllipseX")
                        searchEllipseXChild.text =  searchEllipseDimX[c]
                        searchEllipseYChild = etree.SubElement(root, "searchEllipseY")
                        searchEllipseYChild.set("name","searchEllipseY")
                        searchEllipseYChild.text =  searchEllipseDimX[c]       
                        searchEllipseZChild = etree.SubElement(root, "searchEllipseZ")
                        searchEllipseZChild.set("name","searchEllipseZ")
                        searchEllipseZChild.text =  '1'
                        
                        # Write Min/Max conditioning number
                        minCondNumChild = etree.SubElement(root, "minConditionNum")
                        minCondNumChild.set("name","minCondNum")
                        minCondNumChild.text = minConditionNumber[d]
                        maxCondNumChild = etree.SubElement(root, "maxConditionNum")
                        maxCondNumChild.set("name","maxCondNum")
                        maxCondNumChild.text = maxConditionNumber[0]
                        
                        # Write grid size
                        gridSizeXChild = etree.SubElement(root, "gridSizeX")
                        gridSizeXChild.set("name","gridSizeX")
                        gridSizeXChild.text = gridSize[0]
                        gridSizeYChild = etree.SubElement(root, "gridSizeY")
                        gridSizeYChild.set("name","gridSizeY")
                        gridSizeYChild.text = gridSize[1]
                        
                        # Write finite flag
                        finiteFlagChild = etree.SubElement(root, "UseFiniteKrig")
                        finiteFlagChild.set("name", "UseFiniteKrig")
                        finiteFlagChild.text = useFinite[e]
                        
                        # Write multiple TI flag
                        useMultipleTIChild = etree.SubElement(root, "UseMultipleTI")
                        useMultipleTIChild.set("name", "UseMultipleTI")
                        useMultipleTIChild.text = str(useMultipleTI)
                        
                        # Write penalty terms
                        penaltyFlagChild = etree.SubElement(root, "usePenalty")
                        penaltyFlagChild.set("name", "usePenalty")
                        penaltyFlagChild.text = str(usePenaltyVal[f])
                        
                        # Write R and r terms
                        rFlag = etree.SubElement(root, "r")
                        rFlag.set("name", "r")
                        rFlag.text = str(r)
                        RFlag = etree.SubElement(root, "R")
                        RFlag.set("name", "R")
                        RFlag.text = str(R)
                        
                        
                        
                        # Output To Console
                        print etree.tostring(root, pretty_print=True)
                        
                        # Output To File
                        outputFilename = ['trial' + str(trialNumber[a]) + '.xml']
                        
                        # Run file
                        with open(outputFilename[0], "w") as text_file:
                            text_file.write(etree.tostring(root, pretty_print=True))
                        subprocess.call(["../Release/MPSKriging", outputFilename[0]])
                        
