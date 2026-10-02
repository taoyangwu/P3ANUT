"""
Abundance Ranking Plot data logic.

The Tkinter interface this module used to carry has been replaced by the PyQt
node graph in src/p3anut_ui. What remains is the ranking logic, which the
Ranking Plot block, CLI_rankingPlot.py and any other caller share.
"""

import bisect
import os

import numpy as np


class supportingLogic:


    def  csvComparision(fileA, fileB):
        
        
        fileA_index = []
    
        with open(fileA, 'r') as f:
            for i, line in enumerate(f):
                line_split = line.strip().split(',')
                if(line_split[0] in ["sequence" or "NORMALIZED_ONE_COUNT"]  ):
                    continue
                fileA_index.append((line_split[0], float(line_split[1])))
                
        fileA_index.sort(key=lambda x: x[1], reverse=True)
        fileA_dict = {seq: [rank + 1, freq] for rank, (seq, freq) in enumerate(fileA_index)}
        
                
        fileB_index = []
    
        with open(fileB, 'r') as f:
            for i, line in enumerate(f):
                line_split = line.strip().split(',')
                if(line_split[0] in ["sequence" or "NORMALIZED_ONE_COUNT"]  ):
                    continue
                fileB_index.append((line_split[0], float(line_split[1])))
                
        fileB_index.sort(key=lambda x: x[1], reverse=True)
        fileB_Dictionary = {seq: [rank + 1, freq] for rank, (seq, freq) in enumerate(fileB_index)}
        
        return {"fileA_index": fileA_index, "fileB_index": fileB_index,
                "fileA_Dictionary" : fileA_dict, "fileB_Dictionary": fileB_Dictionary}
        
    
    def gatherScatterData(data, fileA = True, fileB = False, percentOrCount = "%", maskValues = 100):
        
        points = [[], []]
        seqs = []
        
        if(percentOrCount == "%"):
            if(fileA and not fileB):
                x = data["fileA_index"]
                    
                min_i = 0
                for i, (seq, freq) in enumerate(x):
                    if(freq < maskValues / 100):
                        min_i = i
                        break
                    else:
                        min_i += 1
                
                for i, (seq, _) in enumerate(data["fileA_index"][:min_i]):
                    y = data["fileB_Dictionary"].get(seq, -1)
                    if(y == -1):
                        continue
                    
                    points[0].append(i + 1)
                    points[1].append(y[0])
                    seqs.append(seq)
            elif(not fileA and fileB):
                y = data["fileB_index"]
                    
                min_i = 0
                for i, (seq, freq) in enumerate(y):
                    if(freq < maskValues / 100):
                        min_i = i
                        break
                    else:
                        min_i += 1
                
                for i, (seq, _) in enumerate(data["fileB_index"][:min_i]):
                    x = data["fileA_Dictionary"].get(seq, -1)
                    if(x == -1):
                        continue
                    
                    points[0].append(x[0])
                    points[1].append(i + 1)
                    seqs.append(seq)
            elif(fileA and fileB):
                y = data["fileB_index"]
                    
                x_i = 0
                for i, (seq, freq) in enumerate(y):
                    if(freq < maskValues / 100):
                        x_i = i
                        break
                    else:
                        x_i += 1
                        
                x = data["fileA_index"]
                    
                y_i = 0
                for i, (seq, freq) in enumerate(x):
                    if(freq < maskValues / 100):
                        y_i = i
                        break
                    else:
                        y_i += 1
                
                top_x = data["fileA_index"][:x_i]
                top_y = data["fileB_index"][:y_i]
                
                top_x_sequences = set([seq for seq, _ in top_x])
                top_y_sequences = set([seq for seq, _ in top_y])
                
                combined_sequences = top_x_sequences.union(top_y_sequences)
                for seq in combined_sequences:
                    x = data["fileA_Dictionary"].get(seq, -1)
                    y = data["fileB_Dictionary"].get(seq, -1)
                    if(x == -1 or y == -1):
                        continue
                    points[0].append(x[0])
                    points[1].append(y[0])
                    seqs.append(seq)
            else:
                raise ValueError('Invalid Option must include at least one file')
                
        elif(percentOrCount == "#"):
        # else:
            if(fileA and not fileB):
                
                #Get the top N sequences from file A
                print(maskValues)
                x = data["fileA_index"][:int(maskValues)]
                
                #Loop through and get the corresponding value from file B
                for i, (seq, _) in enumerate(x):
                    
                    #Safe guard against missing sequences
                    y = data["fileB_Dictionary"].get(seq, -1)
                    if(y == -1):
                        continue
                    
                    #Append the values to the points list
                    points[0].append(i + 1)
                    points[1].append(y[0])
                    seqs.append(seq)
                
            elif(not fileA and fileB):
                
                #Get the top N sequences from file B
                y = data["fileB_index"][:int(maskValues)]
                
                #Loop through and get the corresponding value from file A
                for i, (seq, _) in enumerate(y):
                    
                    #Safe guard against missing sequences
                    x = data["fileA_Dictionary"].get(seq, -1)
                    if(x == -1):
                        continue
                    
                    #Append the values to the points list
                    points[0].append(x[0])
                    points[1].append(i + 1)
                    seqs.append(seq)
                    
            elif(fileA and fileB):
                
                
                
                top_x = data["fileA_index"][:int(maskValues)]
                top_y = data["fileB_index"][:int(maskValues)]
                
                top_x_sequences = set([seq for seq, _ in top_x])
                top_y_sequences = set([seq for seq, _ in top_y])
                
                combined_sequences = top_x_sequences.union(top_y_sequences)
                for seq in combined_sequences:
                    x = data["fileA_Dictionary"].get(seq, -1)
                    y = data["fileB_Dictionary"].get(seq, -1)
                    if(x == -1 or y == -1):
                        continue
                    points[0].append(x[0])
                    points[1].append(y[0])
                    seqs.append(seq)
                
            else:
                raise ValueError('Invalid Option must include at least one file')
        
        print("Points gathered: ", len(points[0]))
        return points,seqs

    def aboveBelowCounts(x, y, slope, b):
        above, bellow = 0, 0
        
        x_arr = np.array(x)
        y_arr = np.array(y)
        
        y_line = slope * x_arr + b
        
        y_sign = np.sign(y_arr-y_line)
        
        above = np.sum(y_sign > 0) + np.count_nonzero(y_sign == 0)
        bellow = np.sum(y_sign < 0)
        
        return above, bellow
        

    def returnQuadrant(df, slope, b, above = True):
        
        x = np.array(df['FileA_Index'])
        y = x * slope + b
        
        if(above):
            return 0
        elif(not above):
            return 1
        else:
            raise ValueError('Invalid Option must be above (True) or bellow (False)')


    def exportIndices(fileA, indies, exportName, file1_indices=None, file2_indices=None):
        
        lines = []
        writeCounter = 0
        header = ""
        
        with open(fileA, 'r') as f:
            header = f.readline().rstrip('\n')
            for i, line in enumerate(f):
                line_split = line.strip().split(',')
                if(line_split[0] in ["sequence" or "NORMALIZED_ONE_COUNT"]  ):
                    lines.append((line.rstrip('\n'), None, None))
                elif(i in indies):
                    pos = indies.index(i)
                    f1 = file1_indices[pos] if file1_indices is not None else ""
                    f2 = file2_indices[pos] if file2_indices is not None else ""
                    lines.append((line.rstrip('\n'), f1, f2))
                    writeCounter += 1
                
                if(writeCounter >= len(indies)):
                    break
                
        with open(exportName, 'w') as f:
            f.write(f"{header},file1_index,file2_index\n")
            for line, f1, f2 in lines:
                if f1 is None:
                    f.write(f"{line},,\n")
                else:
                    f.write(f"{line},{f1},{f2}\n")
                
        
        
        
    
    def smallest_greater_index(data, value, ascending=True):
        """
        Return index of the smallest element in sorted 'data' that is > value.
        If no element is greater, returns len(data).
        Works with list/tuple/numpy array/pandas Series.
        ascending=True for ascending-sorted data (default).
        """
        # normalize to a Python sequence for bisect
        seq = data if isinstance(data, (list, tuple)) else list(data)

        if not seq:
            return 0

        if ascending:
            # bisect_right gives insertion point after any equals -> first element > value
            return bisect.bisect_right(seq, value)
        else:
            # for descending order, flip sign to reuse bisect on ascending data
            neg_seq = [-x for x in seq]
            return bisect.bisect_right(neg_seq, -value)
