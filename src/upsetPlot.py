"""
Upset Plot set logic.

The Tkinter interface this module used to carry has been replaced by the PyQt
node graph in src/p3anut_ui. What remains is the set comparison logic, kept on
the upsetPlot class so existing callers such as CLI_upsetplot.py continue to
work unchanged.
"""

import os

import numpy as np
import pandas as pd


class upsetPlot:
    """Set intersection helpers shared by the Upset Plot block and the CLI."""

    @staticmethod
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
    
    @staticmethod
    def fileSetComparision(fileNames, intersectionDict = {}):
        maxNumber = 2 ** len(fileNames)
        insersectionConuts = []
        fileLengths = []
        
        totalFileUnion = pd.Index([])
        files = []
        for i in fileNames:
            pdData = pd.read_csv(i)
            fileLengths.append(len(pdData))
            files.append(sequenceValues := pdData.loc[:, "sequence"].values)
            totalFileUnion = totalFileUnion.union(sequenceValues)
            
        totalUnionSize = len(totalFileUnion)
        
        for i in range(maxNumber):
            
            insersectionConuts.append(t := len(upsetPlot.fileInsersection(i, files)))
            intersectionDict[i] = t
            
        return insersectionConuts, fileLengths, totalUnionSize
    
    @staticmethod
    def fileOutput(fileNames, output,key):
        
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
        selectedNames = []
        for i, fileName in enumerate(fileNames):
            pdData = pd.read_csv(fileName)
            files.append(pdData.loc[:, "sequence"].values)
            if(refrenceNumber[i] == 1):
                pdData.set_index("sequence", inplace=True)
                dfs.append(pdData)
                #Only the selected files contribute a column, so their names are
                #tracked alongside rather than indexed back out of fileNames
                selectedNames.append(fileName)
            
            
        t = upsetPlot.fileInsersection(key, files)
        print(f"Number of sequences in output: {len(t)}")
        
        #get the data from the files
        df_counts = []
        for df in dfs:
            df_counts.append(np.array(df.loc[t]["m_index"].values))
            
        df_counts_np = np.array(df_counts)
        
        m_index_counts = np.mean(df_counts_np, axis=0)
        s_index_counts = np.std(df_counts_np, axis=0)
        
        nedDF = pd.DataFrame({"sequence" : t, "m_index" : m_index_counts, "s_index" : s_index_counts})

        #Two selected files can share a name while living in different folders,
        #so a repeated column gets a numeric suffix rather than colliding.
        usedNames = set(nedDF.columns)
        for i, df_count in enumerate(df_counts):
            shortName = os.path.basename(selectedNames[i]).split(".")[0]
            columnName = f"{shortName}_m_index"

            suffix = 2
            while columnName in usedNames:
                columnName = f"{shortName}_m_index ({suffix})"
                suffix += 1

            usedNames.add(columnName)
            nedDF.insert(i + 3, columnName, df_count)
        nedDF.set_index("sequence", inplace=True)
        nedDF.sort_values(by=["m_index"], inplace=True, ascending=False)
        nedDF.to_csv(output)
        print(f"Output saved to {output}")
