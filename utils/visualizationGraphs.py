import matplotlib.pyplot as plt
import numpy as np

def load_distribution_data(file_path):
    ratio_ranges = []
    counts = []
    counts_below = []
    
    midRatio = lambda r: (float(r.split('-')[0]) + (float(r.split('-')[1]))) / 2
    
    with open(file_path, 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split(',')
            ratio_ranges.append(midRatio(parts[0]))
            counts.append(int(parts[1]))
            counts_below.append(int(parts[2]))
    
    return ratio_ranges, counts, counts_below

def plot_distribution(ratio_ranges, counts, counts_below, output_path):
    plt.figure(figsize=(10, 6))
    plt.bar(ratio_ranges, counts, width=0.08, alpha=0.7, label='Count per Range')
    plt.plot(ratio_ranges, counts_below, color='red', marker='o', label='Cumulative Count Below')
    
    plt.xlabel('File1/File2 Ratios within 0.1 intervals')
    plt.ylabel("Number of sequences")
    plt.yscale('log')
    plt.title('Distribution of Average File1/File2 Ratios')
    plt.legend()
    
    plt.grid(True)
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
    # input_file = "/Users/ethankoland/Desktop/FOR ETHAN/demo.txt"
    # output_image = "dev_Tools/volcano_plot_distribution.png"
    
    # ratio_ranges, counts, counts_below = load_distribution_data(input_file)
    # plot_distribution(ratio_ranges, counts, counts_below, output_image)
    
    input_file = "testMerge.npy"
    output_image = "dev_Tools/levenstein_distribution.png"
    graph_levenstein_distribution(input_file, output_image)