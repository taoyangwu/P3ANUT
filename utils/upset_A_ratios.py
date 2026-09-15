import os
import sys

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
SRC_DIR = os.path.join(PROJECT_ROOT, "src")

if SRC_DIR not in sys.path:
    sys.path.insert(0, SRC_DIR)

from volcanoPlot import supportingLogic as vPlotLogic

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def upsetplot(fileNames, output,key):
        
        #Create a binary number where each bit represents if a set is in or not
        refrenceNumber = np.arange(0, len(fileNames))
        refrenceNumber = np.power(2, refrenceNumber)
        
        #Flip the bits so that the array represents the sets
        refrenceNumber = np.flip(refrenceNumber)
        
        #Match the key to the sets
        refrenceNumber = np.bitwise_and(key, refrenceNumber)
        refrenceNumber = np.sign(refrenceNumber)
        
        files = []
        dfs = []
        for i, fileName in enumerate(fileNames):
            pdData = pd.read_csv(fileName)
            files.append(pdData.loc[:, "sequence"].values)
            if(refrenceNumber[i] == 1):
                pdData.set_index("sequence", inplace=True)
                dfs.append(pdData)
            
            
        t = fileInsersection(key, files)
        return len(t)

def fileInsersection(key, fileSets):
        #Create a binary number where each bit represents if a set is in or not
        refrenceNumber = np.arange(0, len(fileSets))
        refrenceNumber = np.power(2, refrenceNumber)
        
        #Flip the bits so that the array represents the sets
        refrenceNumber = np.flip(refrenceNumber)
        
        #Match the key to the sets
        refrenceNumber = np.bitwise_and(key, refrenceNumber)
        refrenceNumber = np.sign(refrenceNumber)
        
        #Seperate the files to be unioned and the files to be differenced
        unionOrDifferencedToggle = np.argsort(refrenceNumber * -1, kind="stable")

        intersectionSet = pd.Index(fileSets[unionOrDifferencedToggle[0]])
        difference = pd.Index([])
        
   
        if(max(refrenceNumber) == 0):
            return []

        for i in range(1, len(unionOrDifferencedToggle)):

            comparision = pd.Index(fileSets[unionOrDifferencedToggle[i]])
            
            t = refrenceNumber[unionOrDifferencedToggle[i]]
 

            if(refrenceNumber[unionOrDifferencedToggle[i]] == 1):
                
                intersectionSet = intersectionSet.intersection(comparision)
            else:
                # s = np.setdiff1d(s, c)
                difference = difference.union(comparision)
                
        # print(len(unions), len(difference))
        return intersectionSet.difference(difference)

def volcanoPlot(fileA, fileB, output, ratio = 0.5, pvalue = 0.05, quadrant = 1):
    data = vPlotLogic.csvComparision(fileA, fileB)

    
        
    quadDf = vPlotLogic.returnQuadrant(data, ratio, pvalue, quadrant)
  
        
    exportDF = vPlotLogic.exportDF(fileA, fileB, quadDf)

    exportDF.to_csv(output)

    return output

def main():
    fileNames = ["data/P3ANUT_R26/UnifiedFiles/m7_man_CONA.csv", "data/P3ANUT_R26/UnifiedFiles/m7_man_BSA.csv", "data/P3ANUT_R26/UnifiedFiles/m7_gal_CONA.csv"]
    
    ratio_range = np.arange(3,9, 0.5)
    counts = np.zeros(len(ratio_range))
    
    for i ,ratio in enumerate(ratio_range):
        man_gal = volcanoPlot(fileNames[0], fileNames[2], "MAN_GAL.csv", ratio=ratio, pvalue=0.05, quadrant=1)
        man_bsa = volcanoPlot(fileNames[0], fileNames[1], "MAN_BSA.csv", ratio=ratio, pvalue=0.05, quadrant=1)
        gal_bsa = volcanoPlot(fileNames[2], fileNames[1], "GAL_BSA.csv", ratio=ratio, pvalue=0.05, quadrant=1)

        counts[i] = upsetplot([man_gal, man_bsa, gal_bsa], "upsetPlot.csv", 7)

    fig, ax = plt.subplots()
    ax.plot(ratio_range, counts, marker="o")
    ax.set_xlabel("Ratio Threshold")
    ax.set_ylabel("Count of Sequences in Intersection")
    ax.set_title("Count of Sequences in Intersection vs Ratio Threshold")
    plt.savefig("ratio_vs_count.png")
    plt.show()

if __name__ == "__main__":
    main()

