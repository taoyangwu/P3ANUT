import matplotlib.pyplot as plt
import numpy as np
import json

def load_distribution_data(file_path):
    #Ratio_Range, Total_Count, Rolling_Count, Count_Below_Pvalue, Count_Above_Pvalue
    ratio_ranges = []
    counts = []
    counts_below = []
    counts_below_pvalue = []
    counts_above_pvalue = []
    
    midRatio = lambda r: (float(r.split('-')[0]) + (float(r.split('-')[1]))) / 2
    
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split(',')
            ratio_ranges.append(midRatio(parts[0]) if '+' not in parts[0] else float(parts[0].split('-')[0]) + 0.05)
            counts.append(int(parts[1]))
            counts_below.append(int(parts[2]))
            counts_below_pvalue.append(int(parts[3]))
            counts_above_pvalue.append(int(parts[4]))
    
    return ratio_ranges, counts, counts_below, counts_below_pvalue, counts_above_pvalue

def plot_distribution(ratio_ranges, counts, counts_below, output_path, counts_below_pvalue, counts_above_pvalue):
    plt.figure(figsize=(10, 6))
    
    plt.plot(ratio_ranges, counts_below, color='red', marker='o', label='Cumulative Count Below')
    
    bar_dict = {
        'Count Below P-Value': counts_below_pvalue,
        'Count Above P-Value': counts_above_pvalue
    }
    bottom = np.zeros(len(ratio_ranges))
    for label, bar_counts in bar_dict.items():
        plt.bar(ratio_ranges, bar_counts, width=0.08, bottom=bottom, alpha=0.7, label=label, 
                color = "yellow" if label == 'Count Above P-Value' else "red")
        bottom += np.array(bar_counts)
     
    xticks = [str(x) for x in np.arange(0, 10.1, 0.5)]
    xticks[-1] = '10+'
    plt.xticks(xticks)
    
    plt.xlabel('A1/A2 Ratios within 0.1 intervals')
    plt.ylabel("Number of sequences")
    plt.yscale('log')
    plt.title('Distribution of Average A1/A2 Ratios')
    plt.legend()
    
    plt.grid(True)
    plt.savefig(output_path)
    
    plt.show()
    
def plot_distribution_2(ratio_ranges, counts, counts_below, output_path, counts_below_pvalue, counts_above_pvalue):
    fig = plt.figure(figsize=(10, 6))
    ax1 = fig.add_subplot(111)
    
    
    ax1.bar(ratio_ranges, counts_below_pvalue, width=0.08, label="Count Below P-Value", color='pink')
    rolling_sum = np.array([sum(counts_below_pvalue[i:]) for i in range(len(counts_below_pvalue))])
    
    
    
    s1 = np.array(counts_above_pvalue) + np.array(counts_below_pvalue)
    
    rolling_s1 = np.array([sum(s1[i:]) for i in range(len(s1))])
    s1 =(rolling_sum / rolling_s1)
    
    ax1.plot(ratio_ranges, rolling_sum, color='red', marker='o', label='Cumulative Count', alpha=0.5)
    ax1.legend()
    
    max_s1 = max(s1)
    max_s1 = round(max_s1, 1)
    
    ax2 = ax1.twinx()

    ax2.plot(ratio_ranges, s1 * 100, color='darkblue', marker='o', label='Ratio of squence count below P-Value to sequence count greater than A1vA2 Ratio')
    ax2.set_ylabel('Percentage of Cumulative Count below P-value', color='darkblue')
    ax2.tick_params(axis='y')
    ax2.set_ylim(bottom=0, top=max_s1 * 100)
    #ax2.legend()
    
    
   
    
     
    xticks = [x for x in np.arange(0, 10.1, 0.5)]
    ax1.xaxis.set_ticks(xticks)
    
    #ax1.xaxis.xlabel('A1/A2 Ratios within 0.1 intervals')
    ax1.set_xlabel('A1/A2 Ratios within 0.1 intervals')
    ax1.set_ylabel("Number of sequences")
    #plt.yscale('log')
    ax1.set_title('Distribution of Average A1/A2 Ratios')
    
    ax1.set_ylim(bottom=1)
    
    ax1.grid(True)
    plt.savefig(output_path)
    
    plt.show()
    
def graph_levenstein_distribution(data_file, output_path):
    data = np.load(data_file, allow_pickle=True)
    data = np.unique(data, return_counts=True)
    print(max(data[0]))
    counts = data[1] / np.sum(data[1])  # Normalize to frequency
    plt.bar(data[0], counts, width=1.0, alpha=0.7)
    plt.xlabel('Number of Insert/Delete Operations')
    plt.ylabel('Frequency')
    plt.title('Distribution of Levenstein Insert/Delete Operations \n in GAL-CA Merged Reads')
    plt.xticks(range(1, max(data[0]) + 1, 2))
    plt.savefig(output_path)
    plt.show()
    
def graph_sequence_length_distribution(data_path, output_path):
    
    data = np.load(data_path, allow_pickle=True)
    
    has_values = np.where(data > 0)[0]
    
    total = np.sum(data)
    
    percentages = (data / total) * 100
    
    fig = plt.figure(figsize=(8, 8))
    
    ax1 = fig.add_subplot(111)
    
    ax1.bar(np.arange(len(data)) + 1, data, width=1.0, alpha=0.7)
    ax1.set_xlabel('Sequence Length (base pairs)')
    ax1.set_xlim(min(has_values)-5, max(has_values)+5)
    ax1.set_ylabel('Count')
    ax1.set_yscale('log')
    
    ax2 = ax1.twinx()
    ax2.plot(np.arange(len(data)) + 1, percentages, color='xkcd:burnt orange', marker='o', label='Percentage of Total Sequences')
    ax2.set_ylabel('Percentage of Total Sequences (%)', color='xkcd:burnt orange')
    ax2.tick_params(axis='y')
    ax2.set_ylim(bottom=0, top=max(percentages[has_values]) + 5)
    
    
    plt.title(f'Distribution of Sequence Lengths - Target length 68 base pairs \n Percentage at Target Length: {percentages[68]:.2f}, Percetage within 1bp: {np.sum(percentages[67:69]):.2f}%')
    plt.savefig(output_path)
    plt.show()
    
def graph_delta_score(json_file):
    data = {}
    with open(json_file, 'r') as f:
        data = json.load(f)

    extracted_data = {d_length : values["both_counts"] / values["total_reads"] for d_length, values in data.items()}
    delta_lengths = list(extracted_data.keys())
    ratios = list(extracted_data.values())
    plt.figure(figsize=(10, 6))
    plt.plot(delta_lengths, ratios, color='blue', marker='o', label='Total Reads to Final Reads Ratio')
    plt.xlabel('Threshold of Quailty Score Difference (δ)')
    plt.ylabel('Percentage of Reads with Fixed Regions')
    # plt.title('Effect of Delta Threshold on Read Retention')
    plt.xticks(delta_lengths)
    plt.grid(True)
    plt.savefig("dev_Tools/delta_score_effect.png")
    plt.show()

def graph_improvements(json_file):
    data = {}
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    forawrds = []
    reverses = []
    mergeds = []
    
    for file, value in data.items():
        forawrds.append(value["forward"]["both_counts"] / value["forward"]["total_reads"])
        reverses.append(value["reverse"]["both_counts"] / value["reverse"]["total_reads"])
        mergeds.append(value["merged"]["both_counts"] / value["merged"]["total_reads"])
        
    data = [forawrds, reverses, mergeds]
    
    forawrds_npy = np.array(forawrds) * 100
    reverses_npy = np.array(reverses) * 100
    mergeds_npy = np.array(mergeds) * 100
    
    print("Forward Reads - Mean: ", np.mean(forawrds_npy), " Std Dev: ", np.std(forawrds_npy))
    print("Reverse Reads - Mean: ", np.mean(reverses_npy), " Std Dev: ", np.std(reverses_npy))
    print("Merged Reads - Mean: ", np.mean(mergeds_npy), " Std Dev: ", np.std(mergeds_npy))
        
    plt.figure(figsize=(10, 6))
    plt.boxplot(data, labels=['Forward Reads', 'Reverse Reads', 'Merged Reads'])
    plt.ylabel('Percentage of Reads with Fixed Regions')
    #plt.title('Improvement in Read Retention After Merging')
    plt.grid(True)
    plt.savefig("dev_Tools/merging_improvement.png")
    plt.show()
    
    
        

if __name__ == "__main__":
    input_file = "dev_Tools/demo.txt"
    output_image = "dev_Tools/volcano_plot_distribution.png"
    
    #Count_Below_Pvalue - Blue and Red seperated by a1va2 ration
    #Count_Above_Pvalue - Yellow and Green seperated by a1va2 ration
    # ratio_ranges, counts, counts_below, counts_below_pvalue, counts_above_pvalue = load_distribution_data(input_file)
    # plot_distribution_2(ratio_ranges, counts, counts_below, output_image, counts_below_pvalue, counts_above_pvalue)
    
    # input_file = "leven_idel_Count.npy"
    # output_image = "dev_Tools/levenstein_distribution.png"
    # graph_levenstein_distribution(input_file, output_image)
    
    # input_file = "sequence_length_distribution.npy"
    # output_image = "dev_Tools/sequence_len_distribution.png"
    # graph_sequence_length_distribution(input_file, output_image)
    
    # input_file = "dev_Tools/p3anut_delta_evaluation.json"
    # graph_delta_score(input_file)
    
    input_file = "dev_Tools/p3anut_evaluation_1.json"
    graph_improvements(input_file)