"""
Run Unifier data logic.

The Tkinter interface this module used to carry has been replaced by the PyQt
node graph in src/p3anut_ui. What remains is the merge routine, which the Run
Unifier block, CLI_runUnifier.py and any other caller share.
"""

import os

import numpy as np
import pandas as pd


def merge(filePaths, outputPath="merged.csv"):
    
    
    if(len(filePaths) < 1):
        return None
    
    DFs = []
    columnNames = ["m_index", "s_index"]
    for path in filePaths:
        DFs.append(pd.read_csv(path))

        #Two inputs can share a file name while living in different folders, so
        #a repeat gets a numeric suffix rather than colliding on the join.
        name = os.path.basename(path)
        if name in columnNames:
            suffix = 2
            while f"{name} ({suffix})" in columnNames:
                suffix += 1
            name = f"{name} ({suffix})"

        columnNames.append(name)

    
        
    currentDF = DFs[0]
    currentDF.set_index('sequence', inplace=True)
    currentDF.rename(columns={'m_index': columnNames[2]}, inplace=True)
    currentDF.drop(['s_index'], axis=1, inplace=True)   
    
    for i in range(1, len(DFs)):
        DF = DFs[i]
        DF.set_index('sequence', inplace=True)
        DF.drop(['s_index'], axis=1, inplace=True)
        #Rename the mean col to count
        DF.rename(columns={'m_index': columnNames[i + 2]}, inplace=True)
        currentDF = currentDF.join(DF, how='outer')
        
      
    currentDF.fillna(0, inplace=True)
    print(currentDF)
    
    colSums = currentDF.sum(axis=0)
    currentDF = currentDF.div(colSums, axis=1)
    
    
    means = np.mean(currentDF.values, axis=1)
    stds = np.std(currentDF.values, axis=1)
    
    #Add the means and stds to the currentDF
    
    currentDF['m_index'] = means
    currentDF['s_index'] = stds
    
    total_row = [ 1/ (len(filePaths) * x) for x in colSums.values]
    t =  np.mean(total_row)
    total_row = ["NORMALIZED_ONE_COUNT", t, 0.0] + total_row
    
    dfColumns = columnNames
    dfColumns.insert(0, "sequence")
    
    header_DF = pd.DataFrame([total_row], columns=dfColumns)
    header_DF.set_index('sequence', inplace=True)
    
    
    currentDF.sort_values(by=['m_index'], inplace=True, ascending=False)
    
    totalDF = pd.concat([header_DF, currentDF])
    totalDF.to_csv(outputPath, index=True, columns = columnNames[1:])
