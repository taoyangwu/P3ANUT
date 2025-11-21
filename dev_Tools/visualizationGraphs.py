import matplotlib.pyplot as plt
import numpy as np

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
    s1 =(rolling_sum / rolling_s1) * rolling_sum[0]
    
    ax1.plot(ratio_ranges, rolling_sum, color='red', marker='o', label='Cumulative Count', alpha=0.5)
    ax1.legend()
    
    ax2 = ax1.twinx()

    ax2.plot(ratio_ranges, s1, color='lightskyblue', marker='o', label='Ratio of Counts above and below P-Value')
    ax2.set_ylabel('Percentage of Counts above and below P-Value (%)', color='lightskyblue')
    ax2.tick_params(axis='y')
    ax2.set_ylim(bottom=0, top=100)
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
    counts = data[1] / np.sum(data[1])  # Normalize to frequency
    plt.bar(data[0], counts, width=1.0, alpha=0.7)
    plt.xlabel('Number of Insert/Delete Operations')
    plt.ylabel('Frequency')
    plt.title('Distribution of Levenstein Insert/Delete Operations \n in GAL-BSA Merged Reads')
    plt.xticks(range(1, max(data[0]) + 1, 2))
    plt.savefig(output_path)
    plt.show()
        

if __name__ == "__main__":
    input_file = "dev_Tools/demo.txt"
    output_image = "dev_Tools/volcano_plot_distribution.png"
    
    ratio_ranges, counts, counts_below, counts_below_pvalue, counts_above_pvalue = load_distribution_data(input_file)
    plot_distribution_2(ratio_ranges, counts, counts_below, output_image, counts_below_pvalue, counts_above_pvalue)
    
    # input_file = "testMerge.npy"
    # output_image = "dev_Tools/levenstein_distribution.png"
    # graph_levenstein_distribution(input_file, output_image)