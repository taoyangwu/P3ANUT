"""
Volcano Plot data logic.

The Tkinter interface this module used to carry has been replaced by the PyQt
node graph in src/p3anut_ui. What remains is the comparison logic, which the
Volcano Plot block, CLI_VolcanoPlot.py and any other caller share.
"""

import os

import numpy as np
import pandas as pd
from scipy import stats as sp


class supportingLogic:


    def  csvComparision(fileA, fileB):
        
        
        dataFrameA = pd.read_csv(fileA, index_col=0)
        dataFrameA.drop([x for x in dataFrameA.columns.values if x not in ["sequence", "m_index", "s_index"]], axis=1, inplace=True)
        
        #Check in the NORMALIZED_ONE_COUNT row is within the file
        if("NORMALIZED_ONE_COUNT" in dataFrameA.index):
            t = dataFrameA.loc["NORMALIZED_ONE_COUNT"]
            dfA_oneCount = dataFrameA.loc["NORMALIZED_ONE_COUNT"][0]
            dataFrameA.drop("NORMALIZED_ONE_COUNT", inplace=True)
        else:
            dfA_oneCount = dataFrameA["m_index"].min()
        
        dataFrameB = pd.read_csv(fileB, index_col=0)
        dataFrameB.drop([x for x in dataFrameB.columns.values if x not in ["sequence", "m_index", "s_index"]], axis=1, inplace=True)
        
        if("NORMALIZED_ONE_COUNT" in dataFrameB.index):
            dfB_oneCount = dataFrameB.loc["NORMALIZED_ONE_COUNT"][0]
            dataFrameB.drop("NORMALIZED_ONE_COUNT", inplace=True)
        else:
            dfB_oneCount = dataFrameB["m_index"].min()
        
        dfa_ColumnCount = 3
        dfb_ColumnCount = 3
        
        
        joined = dataFrameA.join(dataFrameB, how='outer', lsuffix='_a', rsuffix='_b')
        
        joined.fillna(0, inplace=True)
        
        sumA, sumB = joined.iloc[:, 0].sum(), joined.iloc[:, 2].sum()
        
        [statistic, pValue] = sp.ttest_ind_from_stats(joined['m_index_a'], joined['s_index_a'], dfa_ColumnCount, joined['m_index_b'], joined['s_index_b'], dfb_ColumnCount, equal_var=False, alternative='greater')
        
        
        joined['m_index_a'] = joined['m_index_a'].replace(0, dfA_oneCount)
        joined['m_index_b'] = joined['m_index_b'].replace(0, dfB_oneCount)
        
        joined['AvB_Ratio'] = (joined['m_index_a'] / sumA) / (joined['m_index_b'] / sumB)
        joined['-log10(P-Value)'] = -(np.log10(pValue + 1e-10))
        
        #Drop columns
        joined.drop(['s_index_a', 's_index_b', 'm_index_a', 'm_index_b'], axis=1, inplace=True)
        
        
        return joined

    def quarantCounts(df, ratio, pvalue):
        # df1T = df1.
        q1 = df[(df['AvB_Ratio'] >= ratio) & (df['-log10(P-Value)'] >= pvalue)]
        q2 = df[(df['AvB_Ratio'] < ratio) & (df['-log10(P-Value)'] >= pvalue)]
        q3 = df[(df['AvB_Ratio'] < ratio) & (df['-log10(P-Value)'] < pvalue)]
        q4 = df[(df['AvB_Ratio'] >= ratio) & (df['-log10(P-Value)'] < pvalue)]
        
        return [len(q1), len(q2), len(q3), len(q4)]

    def returnQuadrant(df, ratio, pvalue, quadrant = 1):
        if(quadrant == 1):
            return df[(df['AvB_Ratio'] >= ratio) & (df['-log10(P-Value)'] >= pvalue)]
        elif(quadrant == 2):
            return df[(df['AvB_Ratio'] < ratio) & (df['-log10(P-Value)'] >= pvalue)]
        elif(quadrant == 3):
            return df[(df['AvB_Ratio'] < ratio) & (df['-log10(P-Value)'] < pvalue)]
        elif(quadrant == 4):
            return df[(df['AvB_Ratio'] >= ratio) & (df['-log10(P-Value)'] < pvalue)]
        else:
            raise ValueError('Invalid Quadrant Number. Must be 1, 2, 3, or 4.')
        
    def trimmedDF(csv, quarantDF):
        dataFrameA = pd.read_csv(csv)
        dataFrameA.set_index('sequence', inplace=True)
        return dataFrameA.loc[quarantDF.index]

    def createDisto(list, maxValue = 1, minValue = 0, bins = 10):
        disto = np.zeros(bins)
        counts = np.zeros(bins)
        
        for i in list:
            t = int((i - minValue) / (maxValue - minValue) * bins)
            t = min(t, bins-1)
            disto[t] += 1
            counts[0:t] += 1
            
        maxDisto = np.max(disto)
        disto = disto / maxDisto
        counts = counts / counts[0] if counts[0] != 0 else counts
        
        axes = (np.arange(bins) * maxValue)/bins + 1/(2*bins)  
            
        return disto, counts, axes

    def exportDF(fileA, fileB,  quadrantDF):
        
        dataFrameA = pd.read_csv(fileA)
        dataFrameB = pd.read_csv(fileB)
        
        # print(dataFrameA)
        
        dataFrameA.set_index('sequence', inplace=True)
        dataFrameB.set_index('sequence', inplace=True)
        
        #Drop all columns except for m_index and s_index
        dataFrameA.drop([x for x in dataFrameA.columns.values if x not in ["m_index", "s_index"]], axis=1, inplace=True)
        dataFrameB.drop([x for x in dataFrameB.columns.values if x not in ["m_index", "s_index"]], axis=1, inplace=True)
        
        joined = dataFrameA.join(dataFrameB, how='outer', lsuffix='_a', rsuffix='_b')
        
        joined.fillna(0, inplace=True)
        
        joined['m_index'] = joined['m_index_a'] + joined['m_index_b']
        
        
        joined.drop(['s_index_a', 's_index_b', 'm_index_a', 'm_index_b'], axis=1, inplace=True)
        joined.insert(1, 's_index', 0)
        
        trimmed = joined.loc[quadrantDF.index]
        trimmed.sort_values(by='m_index', inplace=True, ascending=False)
        print(trimmed)
        
        
        return trimmed
    
    def createDistrabution(fileA, fileB, fileName):
        joined = supportingLogic.csvComparision(fileA, fileB)
        
        ratioCount = {}
        
        start, end, step = 0, 10, 0.1
        stepDecimal = len(str(step).split(".")[-1])
        rolling_count, rolling_above, rolling_below = 0, 0, 0
        for i in np.linspace(start, end, int((end - start)/step) + 1):
            
            count_withinRange_belowPvalue = len(joined[(joined['AvB_Ratio'] >= i) & (joined['AvB_Ratio'] < i + step) & (joined['-log10(P-Value)'] >= -np.log10(0.05))])
            count_withinRange_abovePvalue = len(joined[(joined['AvB_Ratio'] >= i) & (joined['AvB_Ratio'] < i + step) & (joined['-log10(P-Value)'] < -np.log10(0.05))])
            
            count_withinRange = count_withinRange_belowPvalue + count_withinRange_abovePvalue
            ratioCount[float(i)] = {"belowPvalue": count_withinRange_belowPvalue, "abovePvalue": count_withinRange_abovePvalue, "total": count_withinRange}
            rolling_count += count_withinRange
            rolling_below += count_withinRange_belowPvalue
            rolling_above += count_withinRange_abovePvalue
        
        count_withinRange_belowPvalue = len(joined[(joined['AvB_Ratio'] > i + step) & (joined['-log10(P-Value)'] >= -np.log10(0.05))])
        count_withinRange_abovePvalue = len(joined[(joined['AvB_Ratio'] > i + step) & (joined['-log10(P-Value)'] < -np.log10(0.05))])
        
        count_withinRange = count_withinRange_belowPvalue + count_withinRange_abovePvalue
        
        ratioCount[float(end)] = {"belowPvalue": count_withinRange_belowPvalue, "abovePvalue": count_withinRange_abovePvalue, "total": count_withinRange}
        rolling_count += count_withinRange
        rolling_below += count_withinRange_belowPvalue
        rolling_above += count_withinRange_abovePvalue
            
            
        with open(fileName, 'w') as f:
            f.write("Ratio_Range, Total_Count, Rolling_Count, Count_Below_Pvalue, Count_Above_Pvalue\n")
             
            for key, value in ratioCount.items():
                roundedKey = round(key, stepDecimal)
                if(key == end):
                    f.write(f"{roundedKey}-{roundedKey}+, {value['total']}, {rolling_count}, {value['belowPvalue']}, {value['abovePvalue']}\n")
                else:
                    f.write(f"{roundedKey}-{round(key + step, stepDecimal)}, {value['total']}, {rolling_count}, {value['belowPvalue']}, {value['abovePvalue']}\n")
                rolling_count -= ratioCount[key]["total"]
            
        pass
