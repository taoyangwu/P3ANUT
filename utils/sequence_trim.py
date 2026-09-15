import os
import csv
import shutil


def trip_count_file(filePath, start_bases = 18,end_bases =0, output = "trimmed.csv"):
    with open(filePath, "r") as f:
        reader = csv.reader(f)
        lines = list(reader)

    header = lines[0]
    one_count = lines[1]

    file_counts = {}

    for line in lines[2:]:
        seq, count, std = line

        trimmedSeq = seq[start_bases:-end_bases] if end_bases > 0 else seq[start_bases:]
        if(trimmedSeq in file_counts):
            file_counts[trimmedSeq] += float(count)
        else:
            file_counts[trimmedSeq] = float(count)

    with open(output, "w") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerow(one_count)

        for seq, count in file_counts.items():
            writer.writerow([seq, count, 0])

if(__name__ == "__main__"):
    files = [f"data/Rhau/Reverse/R{i}_Rhau18_12aa_R/Rhau18_12aa_R.csv" for i in [3,4,5,6,7]]
    files_new = [f"data/Rhau/Reverse/R{i}_Rhau18_12aa_R/R{i}_Rhau18_12aa_R.csv" for i in [3,4,5,6,7]]
    files2 = [f"data/Rhau/Forward/R{i}_Rhau18_12aa_F/Rhau18_12aa_F.csv" for i in [3,4,5,6,7]]
    files2_new = [f"data/Rhau/Forward/R{i}_Rhau18_12aa_F/R{i}_Rhau18_12aa_F.csv" for i in [3,4,5,6,7]]
    f3 = [f"data/Rhau/matched/R{i}_Rhau18_12aa/Amino_30.csv" for i in [3,4,5,6,7]]
    f3_new = [f"data/Rhau/matched/R{i}_Rhau18_12aa/R{i}_A30.csv" for i in [3,4,5,6,7]]


    # for file, new_file in zip(files, files_new):
    #     shutil.copy(file, new_file)

    # # for file, new_file in zip(files2, files2_new):
    # #     shutil.copy(file, new_file) 

    # # for file, new_file in zip(f3, f3_new):
    # #     shutil.copy(file, new_file)

    # for file in files_new:
    #     output = file[:-4] + "_trimmed.csv"
    #     trip_count_file(file, start_bases = 18, output = output)    

    # for file in files2_new:
    #     output = file[:-4] + "_trimmed.csv"
    #     trip_count_file(file, start_bases = 18, output = output)   

    for file in f3_new:
        #output = file[:-4] + "_trimmed.csv"
        trip_count_file(file, start_bases=3, end_bases=3, output = file) 
