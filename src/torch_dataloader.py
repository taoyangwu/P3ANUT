import torch
import re


class _P3ANUT_Dataset(torch.utils.data.Dataset):
    def __init__(self, data, **kwargs):
        self.data = data

        #Set up conversion arrays that will convert the bases into ints. Additionally it provides additional protection against invalid bases
        self.baseConversionArray = torch.full(ord('T') + 1, 0, dtype=torch.uint8)
        self.baseConversionArray[ord('T')] = 1
        self.baseConversionArray[ord('G')] = 2
        self.baseConversionArray[ord('C')] = 3
        self.baseConversionArray[ord('A')] = 4
        
        #Seperate conversion array is needed when the sequence is flipped
        self.flipBaseConversionArray = torch.full(ord('T') + 1, 0, dtype=torch.uint8)
        self.flipBaseConversionArray[ord('T')] = 4
        self.flipBaseConversionArray[ord('G')] = 3
        self.flipBaseConversionArray[ord('C')] = 2
        self.flipBaseConversionArray[ord('A')] = 1

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        return self.data[idx][0], self._convert_sequence(self.data[idx][1], self.data[idx][2], flip=False), \
               self._convert_sequence(self.data[idx][3], self.data[idx][4], flip=True)
    
    def _convert_sequence(self, p_sequence, p_score, flip=False):
        
        #Prevent edge conditions where the sequence and the score are not the same length
        if(len(p_sequence) != len(p_score)):
            #mismatchCount += 1
            return None, None, None
        
        #Convert the sequence into a numpy array of ints
        sequnceByteArray = bytearray(p_sequence, 'utf-8')
        
        
        #Convert the scores into a numpy array of ints
        qualityByteArray = bytearray(p_score, 'utf-8')
        
        
        if(flip):
            sequence = torch.frombuffer(reversed(sequnceByteArray), dtype=torch.uint8)
            qualityScore = torch.frombuffer(reversed(qualityByteArray), dtype=torch.uint8)
            torch.take(self.flipBaseConversionArray, sequence, out=sequence)
        else:
            sequence = torch.frombuffer(sequnceByteArray, dtype=torch.uint8)
            qualityScore = torch.frombuffer(qualityByteArray, dtype=torch.uint8)
            torch.take(self.baseConversionArray, sequence, out=sequence)
        # sequence = codonConversionFunction(sequence).astype(np.uint8)

        
        #The quality score is offset by 33 in the fastq file, so we need to subtract it to get the actual score and allow for compression  
        torch.subtract(qualityScore, 33, out=qualityScore)
        
        #Force the scores of the N bases to be 0
        NbaseIndex = torch.where(sequence == 0)
        torch.put(qualityScore, NbaseIndex, 0)
        
        #Modify the quality score instead of a new array to make slightly more memeory efficent
        torch.multiply(qualityScore, 5, out=qualityScore)
        torch.add(qualityScore, sequence, out=qualityScore)
        
        return qualityScore
    
def _read_FASTQ(self, filePath, **kwargs):

    cull_maxlength = kwargs.get("cull_maxlength",90)
    cull_minlength = kwargs.get("cull_minlength",20)
    includeN = kwargs.get("includeN",False)
    fileSeperator = kwargs.get("fileSeperator",r"\+")

    with open(filePath, "r") as file:
    #Example Entry after regex as a tuple
    #('NB501061:163:HVYLLAFX3:1:11101:1980:1063' - Run ID, 
    # '1' - Run Number, 
    # ':N:0:GTATTATCT+CATATCGTT' - Additional Run Information, 
    # 'TGTAGACTATTCTCACTCTTCTTGTCTGGTTCCTCCGCGTCCGACGTGTGGTGGAGGTTCGGTCGACG', - DNA Sequence 
    # 'AAAAAEE<EE<EEEEEEEEEEAEEEE/EEAEEEEA//AA/EAAEEEEEEEEAEEEEEE/EA<EEEEE6' - DNA Quality Score)
        regexExpression = r"@([A-Z0-9:]+)(?:\s)(\d)([A-Z0-9:+]+)(?:\s*)([CODONS]{cull_minlength,cull_maxlength})(?:\s+SPLIT\s*)([!-I]{cull_minlength,cull_maxlength})".replace("cull_minlength", str(cull_minlength)).replace("cull_maxlength", str(cull_maxlength))
        regexExpression = regexExpression.replace("CODONS", "ATGCN" if includeN else "ATGC")
        
        cleanedFileSeparator = fileSeperator.replace("\\\\", "\\")
        regexExpression = regexExpression.replace("SPLIT", f"[{cleanedFileSeparator}]?")
        entries = re.findall(regexExpression, file.read())

    return {entry[0]: (entry[3], entry[4], len(entry[3])) for entry in entries if len(entry[3]) ==  len(entry[4])}

def _split_datasets(forwardData, ReverseData, **kwargs):
    
    forward_keys = set(forwardData.keys())
    reverse_keys = set(ReverseData.keys())

    # Find common keys in both dictionaries
    common_keys = forward_keys.intersection(reverse_keys)

    # keys only in one file
    forward_only_keys = forward_keys - common_keys
    reverse_only_keys = reverse_keys - common_keys

    matched_data, indel_data, forward_only_data, reverse_only_data = [], [], [], []

    insert_lambda = lambda data, key: (key, data[key][0], data[key][1])
    forward_only_data = [insert_lambda(forwardData, key) for key in forward_only_keys]
    reverse_only_data = [insert_lambda(ReverseData, key) for key in reverse_only_keys]

    for key in common_keys:
        forward_entry = forwardData[key]
        reverse_entry = ReverseData[key]

        if forward_entry[2] == reverse_entry[2]:
            matched_data.append((key, forward_entry[0], forward_entry[1], reverse_entry[0], reverse_entry[1]))
        else:
            indel_data.append((key, forward_entry[0], forward_entry[1], reverse_entry[0], reverse_entry[1]))

    return matched_data, indel_data, forward_only_data, reverse_only_data

        
    
def return_dataloader(forwardPath, revesePath, batch_size=32, **kwargs):

    forwardData = _read_FASTQ(forwardPath, **kwargs)
    reverseData = _read_FASTQ(revesePath, **kwargs)

    matched_data, indel_data, forward_only_data, reverse_only_data = _split_datasets(forwardData, reverseData, **kwargs)

    create_dataloader= lambda data, batchsize: torch.utils.data.DataLoader(_P3ANUT_Dataset(data), batch_size=batch_size)
    
    return create_dataloader(matched_data, batch_size), \
           create_dataloader(indel_data, batch_size), \
           create_dataloader(forward_only_data, batch_size), \
           create_dataloader(reverse_only_data, batch_size)

if __name__ == "__main__":
    # Example usage
    forwardPath = "path/to/forward.fastq"
    reversePath = "path/to/reverse.fastq"
    
    dataloader, indel_dataloader, forward_only_dataloader, reverse_only_dataloader = return_dataloader(forwardPath, reversePath, batch_size=32)
    
    for data in dataloader:
        print(data)