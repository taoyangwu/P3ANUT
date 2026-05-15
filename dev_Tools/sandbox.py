import json

file_lengths_1 = {
}

with open("file_lengths.json", 'r') as f:
    file_lengths_1 = json.load(f)

file_lengths_2 = {
}

with open("rebollo_file_lengths.json", 'r') as f:
    file_lengths_2 = json.load(f)

overlap = set(file_lengths_1.keys()) & set(file_lengths_2.keys())

joined_lengths =  file_lengths_2 | file_lengths_1

with open("joined_file_lengths.json", 'w') as f:
    json.dump(joined_lengths, f, indent=4)