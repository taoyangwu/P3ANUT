import json
import matplotlib.pyplot as plt
import numpy as np

def runtime_comparison(metrics_path, file_sizes):

    with open(metrics_path, 'r') as f:
        data = json.load(f)

    with open(file_sizes, 'r') as f:
        size_data = json.load(f)

    flash_times = []
    casper_times = []
    p3anut_times = []
    run_names = []

    pass

def avg_rebollo(avg_metrics):

    rebollo_avg = avg_metrics["rebollo"]
    rebollo_reverse_avg = avg_metrics["rebollo_reverse"]

    combined_avg = {}
    for metric in rebollo_avg.keys():
        combined_avg[metric] = (rebollo_avg[metric] + rebollo_reverse_avg[metric]) / 2

    return combined_avg


def calc_averages(metrics_path, print_results = True):

    with open(metrics_path, 'r') as f:
        data = json.load(f)

    base = {
        "time_taken": 0,
        "retention_rate": 0,
        "tau_score": 0,
        "upsilon_score": 0,
        "phi_score": 0,
        "sequence_length_score": 0
    }

    p3anut_avg = base.copy()
    flash_avg = base.copy()
    casper_avg = base.copy()
    rebollo_avg = base.copy()
    rebollo_reverse_avg = base.copy()
    forward_avg = base.copy()
    reverse_avg = base.copy()

    for run_name, metrics in data.items():
        for tool in ["p3anut", "flash", "casper", "rebollo", "rebollo_reverse", "forward", "reverse"]:
            for metric in base.keys():
                if tool in metrics and metric in metrics[tool]:
                    if metrics[tool][metric] is not None:
                        if tool == "p3anut":
                            p3anut_avg[metric] += metrics[tool][metric]
                        elif tool == "flash":
                            flash_avg[metric] += metrics[tool][metric]
                        elif tool == "casper":
                            casper_avg[metric] += metrics[tool][metric]
                        elif tool == "rebollo":
                            rebollo_avg[metric] += metrics[tool][metric]
                        elif tool == "forward":
                            forward_avg[metric] += metrics[tool][metric]
                        elif tool == "reverse":
                            reverse_avg[metric] += metrics[tool][metric]
                        elif tool == "rebollo_reverse":
                            rebollo_reverse_avg[metric] += metrics[tool][metric]

    num_runs = len(data)
    for metric in base.keys():
        p3anut_avg[metric] /= num_runs
        flash_avg[metric] /= num_runs
        casper_avg[metric] /= num_runs
        rebollo_avg[metric] /= num_runs
        forward_avg[metric] /= num_runs
        reverse_avg[metric] /= num_runs
        rebollo_reverse_avg[metric] /= num_runs

    total_score = lambda x: x["retention_rate"] * x["tau_score"] *  x["sequence_length_score"] 

    if print_results:

        print("Average Metrics:")
        print("P3ANUT:", p3anut_avg)
        print("FLASH:", flash_avg)
        print("CASPER:", casper_avg)
        print("REBOLLO:", rebollo_avg)
        print("REBOLLO REVERSE:", rebollo_reverse_avg)
        print("FORWARD:", forward_avg)
        print("REVERSE:", reverse_avg)  

        print("\nAverage Total Scores:")
        print("P3ANUT:", total_score(p3anut_avg))
        print("FLASH:", total_score(flash_avg))
        print("CASPER:", total_score(casper_avg))
        print("REBOLLO:", total_score(rebollo_avg))
        print("REBOLLO REVERSE:", total_score(rebollo_reverse_avg))
        print("FORWARD:", total_score(forward_avg))
        print("REVERSE:", total_score(reverse_avg))

    return {
        "p3anut": p3anut_avg,
        "flash": flash_avg,
        "casper": casper_avg,
        "rebollo": rebollo_avg,
        "rebollo_reverse": rebollo_reverse_avg,
        "forward": forward_avg,
        "reverse": reverse_avg
    }

def compare_table(avg_metrics, non_included_metrics = [], inverse_metrics = set(["time_taken"])):

    def _permute_ranking(sorted_values):
        offset = 1
        adjusted_rankings = np.zeros_like(sorted_values)
        for i in range(len(sorted_values)):
            if i > 0 and sorted_values[i] == sorted_values[i-1]:
                adjusted_rankings[i] = adjusted_rankings[i-1]
            else:
                adjusted_rankings[i] = offset
            offset += 1

        

    tools = list(avg_metrics.keys())
    t1 = tools[0]
    t2 = avg_metrics[tools[0]]
    metircs = list(avg_metrics[tools[0]].keys())

    for metric in metircs:
        if metric in non_included_metrics:
            print(f"removing {metric} from comparison")
            metircs.remove(metric)

    score_table = np.zeros((len(tools), len(metircs)))
    rankings = np.zeros((len(tools), len(metircs)), dtype=int)
    percentage_table = np.zeros((len(tools), len(metircs)))
    
    for j, metric in enumerate(metircs):
        metric_values = np.array([avg_metrics[tool][metric] for tool in tools]) 
        sorted_indices = np.argsort(metric_values * (-1 if metric in inverse_metrics else 1))
        sorted_metric_values = metric_values[sorted_indices]

        offset = 1
        adjusted_rankings = np.zeros_like(sorted_indices)
        for i in range(len(sorted_metric_values)):
            if i > 0 and sorted_metric_values[i] == sorted_metric_values[i-1]:
                adjusted_rankings[i] = adjusted_rankings[i-1]
            else:
                adjusted_rankings[i] = offset
            offset += 1

        rankings[sorted_indices,j] = adjusted_rankings
        score_table[:, j] = metric_values
        percentage_table[:, j] = metric_values / np.max(metric_values)
        
        
    total_scores = np.sum(rankings, axis=1)
    percentage_total_scores = np.sum(percentage_table, axis=1)
    percentage_total_scores /= np.max(percentage_total_scores)
    total_scores_rankings = np.argsort(total_scores)

    score_table = np.round(np.hstack((score_table, total_scores[:, np.newaxis])), 3)
    rankings = np.hstack((rankings, total_scores_rankings[:, np.newaxis]))
    percentage_table = np.round(np.hstack((percentage_table, percentage_total_scores[:, np.newaxis])),3)



    plt.figure(figsize=(10, 10))
    plt.imshow(rankings, cmap="viridis", aspect="auto")
    plt.colorbar(label="Rank")

    #Plot the rank text on the heatmap
    for i in range(rankings.shape[0]):
        for j in range(rankings.shape[1]):
            plt.text(j, i, f"{rankings[i,j]}\n{score_table[i,j]}\n{percentage_table[i,j]}", ha="center", va="center", color="white")

    metircs.append("total_score")

    plt.xticks(ticks=range(len(metircs)), labels=metircs, rotation=45, ha="right")
    plt.yticks(ticks=range(len(tools)), labels=tools)
    plt.title("Tool Rankings Across Metrics")
    plt.tight_layout()
    plt.show()


def main():
    avgs = calc_averages("joined_data.json")
    compare_table(avgs)



if __name__ == "__main__":
    main()



