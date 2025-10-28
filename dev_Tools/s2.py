


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
        else:
            raise ValueError('Invalid Option must include at least one file')
            
    # elif(percentOrCount == "#"):
    else:
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
            
        else:
            raise ValueError('Invalid Option must include at least one file')
    
    print("Points gathered: ", len(points[0]))
    return points

data = csvComparision("dev_Tools/P4.merge.fa.csv", "dev_Tools/P1.merge.fa.csv")

points = gatherScatterData(data, fileA=True, fileB=False, percentOrCount="%", maskValues=0.1)