import matplotlib.pyplot as plt

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
    
    plt.xlabel('Average B/A Ratio Ranges')
    plt.ylabel('Counts')
    plt.yscale('log')
    plt.title('Distribution of Average B/A Ratios')
    plt.legend()
    
    plt.grid(True)
    
    plt.show()
    
if __name__ == "__main__":
    input_file = "/Users/ethankoland/Desktop/FOR ETHAN/demo.txt"
    output_image = "dev_Tools/volcano_plot_distribution.png"
    
    ratio_ranges, counts, counts_below = load_distribution_data(input_file)
    plot_distribution(ratio_ranges, counts, counts_below, output_image)