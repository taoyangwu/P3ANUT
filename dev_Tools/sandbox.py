import json
import matplotlib.pyplot as plt
import numpy as np
import os

def runtime_comparison(metrics_data, file_sizes, save_path=None, show_plot=True):

    with open(file_sizes, 'r') as f:
        size_data = json.load(f)

    flash_times = []
    casper_times = []
    p3anut_times = []
    rebollo_times = []
    sequence_lengths = []

    for run_name, metrics in metrics_data.items():
        if "flash" in metrics and "time_taken" in metrics["flash"]:
            flash_times.append(metrics["flash"]["time_taken"])
        else:
            flash_times.append(None)

        if "casper" in metrics and "time_taken" in metrics["casper"]:
            casper_times.append(metrics["casper"]["time_taken"])
        else:
            casper_times.append(None)

        if "p3anut" in metrics and "time_taken" in metrics["p3anut"]:
            p3anut_times.append(metrics["p3anut"]["time_taken"])
        else:
            p3anut_times.append(None)

        if "rebollo" in metrics and "time_taken" in metrics["rebollo"]:
            rebollo_times.append(metrics["rebollo"]["time_taken"] * 2)
        else:
            rebollo_times.append(None)

        if run_name in size_data:
            sequence_lengths.append(size_data[run_name]["total_length"])
        else:            
            sequence_lengths.append(None)

    plt.figure(figsize=(10, 6))
    plt.scatter(sequence_lengths, flash_times, label="FLASH", color='blue', marker='o')
    plt.scatter(sequence_lengths, casper_times, label="CASPER", color='orange', marker='s')
    plt.scatter(sequence_lengths, p3anut_times, label="P3ANUT", color='green', marker='^')
    plt.scatter(sequence_lengths, rebollo_times, label="REBOLLO", color='red', marker='x')
    plt.xlabel("Sequence Length")
    plt.ylabel("Time Taken (seconds)")
    plt.title("Runtime Comparison of Tools")
    plt.legend()
    plt.xscale('log')
    #plt.yscale('log')
    plt.grid(True, which="both", ls="--", linewidth=0.5)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300)
        print(f"Runtime comparison plot saved to {save_path}")
    if show_plot:
        plt.show()

    plt.cla()  
   

    pass

def avg_rebollo(avg_metrics):

    rebollo_avg = avg_metrics["rebollo"]
    rebollo_reverse_avg = avg_metrics["rebollo_reverse"]

    combined_avg = {}
    for metric in rebollo_avg.keys():
        combined_avg[metric] = (rebollo_avg[metric] + rebollo_reverse_avg[metric]) / 2

    return combined_avg

def load_data(metrics_path, included_fileNames = None):

    with open(metrics_path, 'r') as f:
        data = json.load(f)

    if included_fileNames is not None:
        filtered_data = {k: v for k, v in data.items() if k in included_fileNames}
        return filtered_data

    return data


def calc_averages(metrics_data, print_results = True):

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

    for run_name, metrics in metrics_data.items():
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

    num_runs = len(metrics_data)
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

def calculate_total_scores(avg_metrics):
    total_scores = {}
    for tool, metrics in avg_metrics.items():
        total_scores[tool] = metrics["retention_rate"] * metrics["tau_score"] *  metrics["sequence_length_score"]

    for tool, score in total_scores.items():
        avg_metrics[tool]["total_score"] = score

def compare_table(avg_metrics, non_included_metrics = [], non_included_tools = [], inverse_metrics = set(["time_taken"]), include_ranking_sum = False):

    def _permute_ranking(values, invert=False):

        sorted_indices = np.argsort(values * (-1 if invert else 1))
        sorted_metric_values = metric_values[sorted_indices]

        offset = 1
        adjusted_rankings = np.zeros_like(values, dtype=int)
        for i in range(len(sorted_metric_values)):
            if i > 0 and sorted_metric_values[i] == sorted_metric_values[i-1]:
                adjusted_rankings[i] = adjusted_rankings[i-1]
            else:
                adjusted_rankings[i] = offset
            offset += 1

        return adjusted_rankings, sorted_indices

    tools = [t for t in avg_metrics if t not in non_included_tools]
    metircs = [m for m in avg_metrics[tools[0]] if m not in non_included_metrics]

    for tool in tools:
        if tool in non_included_tools:
            print(f"removing {tool} from comparison")
            tools.remove(tool)

    score_table = np.zeros((len(tools), len(metircs)))
    rankings = np.zeros((len(tools), len(metircs)), dtype=int)
    percentage_table = np.zeros((len(tools), len(metircs)))
    
    for j, metric in enumerate(metircs):
        metric_values = np.array([avg_metrics[tool][metric] for tool in tools]) 
        adjusted_rankings, sorted_indices = _permute_ranking(metric_values, invert=metric in inverse_metrics)
        rankings[sorted_indices, j] = adjusted_rankings
        score_table[:, j] = metric_values
        percentage_table[:, j] = metric_values / np.max(metric_values)
        
    if include_ranking_sum:

        metircs.append("sum_of_ranks")
        
        total_scores = np.sum(rankings, axis=1)
        total_scores_rankings, sorted_indices = _permute_ranking(total_scores, invert=True)

        score_table = np.hstack((score_table, total_scores[:, np.newaxis]))
        rankings = np.hstack((rankings, total_scores_rankings[:, np.newaxis]))

        percentage_total_scores = np.sum(percentage_table, axis=1)
        percentage_total_scores /= np.max(percentage_total_scores)

        percentage_table = np.hstack((percentage_table, percentage_total_scores[:, np.newaxis]))

    score_table = np.round(score_table, 3)
    percentage_table = np.round(percentage_table,3)



    plt.figure(figsize=(10, 10))
    plt.imshow(rankings, cmap="viridis", aspect="auto")

    #Plot the rank text on the heatmap
    for i in range(rankings.shape[0]):
        for j in range(rankings.shape[1]):
            plt.text(j, i, f"{rankings[i,j]}\n{score_table[i,j]}\n{percentage_table[i,j]}", ha="center", va="center", color="white")

    plt.xticks(ticks=range(len(metircs)), labels=metircs, rotation=45, ha="right")
    plt.yticks(ticks=range(len(tools)), labels=tools)
    plt.title("Tool Rankings Across Metrics")
    plt.tight_layout()
    plt.savefig("tool_comparison_heatmap_3.png")
    # plt.show()


def plot_metric_barplot(metrics_data, metric_name, excluded_tools=None, excluded_metrics=None, 
                        output_path=None, figsize=(14, 6), show_plot=True):
    
    if excluded_tools is None:
        excluded_tools = []
    
    if excluded_metrics is None:
        excluded_metrics = []
    
    # Check if metric exists in data
    metric_found = False
    for file_name, tools_dict in metrics_data.items():
        for tool_name, metrics_dict in tools_dict.items():
            if metric_name in metrics_dict:
                metric_found = True
                break
        if metric_found:
            break
    
    if not metric_found:
        raise ValueError(f"Metric '{metric_name}' not found in any file or tool data")
    
    # Extract file names and collect available tools
    file_names = list(metrics_data.keys())
    all_tools = set()
    
    # Gather all available tools
    for file_name, tools_dict in metrics_data.items():
        for tool_name in tools_dict.keys():
            if tool_name not in excluded_tools:
                all_tools.add(tool_name)
    
    all_tools = sorted(list(all_tools))
    
    if len(all_tools) == 0:
        raise ValueError("No tools available after applying exclusions")
    
    # Prepare data for plotting
    x = np.arange(len(file_names))
    width = 0.8 / len(all_tools)  # Width of each bar
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Define colors for tools
    colors = plt.cm.Set3(np.linspace(0, 1, len(all_tools)))
    
    # Plot bars for each tool
    for i, tool in enumerate(all_tools):
        values = []
        for file_name in file_names:
            if tool in metrics_data[file_name] and metric_name in metrics_data[file_name][tool]:
                val = metrics_data[file_name][tool][metric_name]
                values.append(val if val is not None else 0)
            else:
                values.append(0)
        
        ax.bar(x + i * width, values, width, label=tool, color=colors[i])
    
    # Customize plot
    ax.set_xlabel('File Name', fontsize=12, fontweight='bold')
    ax.set_ylabel(metric_name.replace('_', ' ').title(), fontsize=12, fontweight='bold')
    ax.set_title(f'{metric_name.replace("_", " ").title()} Comparison Across Files and Tools', 
                 fontsize=14, fontweight='bold')
    
    # Set x-axis ticks and labels with 45 degree rotation
    ax.set_xticks(x + width * (len(all_tools) - 1) / 2)
    ax.set_xticklabels(file_names, rotation=45, ha='right')
    
    # Add legend and grid
    ax.legend(title='Tools', loc='upper left', fontsize=10)
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    
    # Save figure if path is provided
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {output_path}")
    
    if show_plot:
        plt.show()
    
    return fig, ax


def plot_metric_stacked_bar_with_error(metrics_data, metric_names=None, excluded_tools=None,
                                       excluded_metrics=None, output_path=None,
                                       figsize=(14, 8), show_plot=True, error_type="sem",
                                       include_total_score=True):
    """
    Create a vertical grouped bar plot for one or more metrics with error bars.

    The x-axis shows the metrics being compared. For each metric, each tool has a separate
    bar. If the input is raw run data, the function computes the mean value per tool and
    metric, then adds error bars using either standard deviation or standard error of the mean.

    Parameters:
    1. metrics_data: dict
        Either raw run data shaped like {file_name: {tool: {metric: value}}}
        or aggregated data shaped like {tool: {metric: value}}.
    2. metric_names: str | list[str] | None
        Metric name or metric names to plot. If None, all metrics are plotted.
    3. excluded_tools: list[str] | None
        Tools to exclude from the plot.
    4. excluded_metrics: list[str] | None
        Metrics to exclude from the plot.
    5. output_path: str | None
        Optional path to save the figure.
    6. figsize: tuple
        Figure size.
    7. show_plot: bool
        Whether to show the figure.
    8. error_type: str
        "sem" for standard error of the mean or "std" for standard deviation.
    9. include_total_score: bool
        If True, include a derived metric:
        total_score = retention_rate * tau_score * sequence_length_score

    Returns:
    fig, ax
        Matplotlib figure and axes.
    """

    if excluded_tools is None:
        excluded_tools = []

    if excluded_metrics is None:
        excluded_metrics = []

    def _looks_like_raw_data(data):
        first_value = next(iter(data.values()), None)
        if not isinstance(first_value, dict) or len(first_value) == 0:
            return False
        first_nested_value = next(iter(first_value.values()), None)
        return isinstance(first_nested_value, dict)

    def _collect_metrics(data, raw_data):
        metric_set = set()
        if raw_data:
            for run_metrics in data.values():
                for tool_metrics in run_metrics.values():
                    if isinstance(tool_metrics, dict):
                        metric_set.update(tool_metrics.keys())
        else:
            for tool_metrics in data.values():
                if isinstance(tool_metrics, dict):
                    metric_set.update(tool_metrics.keys())
        return sorted(metric_set)

    def _collect_tools(data, raw_data):
        tool_set = set()
        if raw_data:
            for run_metrics in data.values():
                tool_set.update(run_metrics.keys())
        else:
            tool_set.update(data.keys())
        return sorted([tool for tool in tool_set if tool not in excluded_tools])

    def _build_summary_from_raw(data, selected_metrics, selected_tools):
        summary = {metric: {} for metric in selected_metrics}

        if include_total_score and "total_score" in selected_metrics:
            for tool in selected_tools:
                total_values = []
                for run_metrics in data.values():
                    tool_metrics = run_metrics.get(tool, {})
                    r = tool_metrics.get("retention_rate", None)
                    t = tool_metrics.get("tau_score", None)
                    s = tool_metrics.get("sequence_length_score", None)
                    if r is not None and t is not None and s is not None:
                        total_values.append(r * t * s)

                if len(total_values) == 0:
                    summary["total_score"][tool] = {"value": None, "error": 0.0}
                else:
                    total_values = np.asarray(total_values, dtype=float)
                    if error_type == "std":
                        total_error = float(np.std(total_values, ddof=1)) if len(total_values) > 1 else 0.0
                    else:
                        total_error = float(np.std(total_values, ddof=1) / np.sqrt(len(total_values))) if len(total_values) > 1 else 0.0

                    summary["total_score"][tool] = {
                        "value": float(np.mean(total_values)),
                        "error": total_error,
                    }

        for metric in selected_metrics:
            if metric == "total_score":
                continue
            for tool in selected_tools:
                values = []
                for run_metrics in data.values():
                    tool_metrics = run_metrics.get(tool, {})
                    value = tool_metrics.get(metric, None)
                    if value is not None:
                        values.append(value)

                if len(values) == 0:
                    summary[metric][tool] = {"value": None, "error": 0.0}
                    continue

                values = np.asarray(values, dtype=float)
                if error_type == "std":
                    error_value = float(np.std(values, ddof=1)) if len(values) > 1 else 0.0
                else:
                    error_value = float(np.std(values, ddof=1) / np.sqrt(len(values))) if len(values) > 1 else 0.0

                summary[metric][tool] = {
                    "value": float(np.mean(values)),
                    "error": error_value,
                }
        return summary

    def _build_summary_from_aggregated(data, selected_metrics, selected_tools):
        summary = {metric: {} for metric in selected_metrics}
        for metric in selected_metrics:
            for tool in selected_tools:
                if metric == "total_score" and include_total_score:
                    tool_metrics = data.get(tool, {})
                    r = tool_metrics.get("retention_rate", None)
                    t = tool_metrics.get("tau_score", None)
                    s = tool_metrics.get("sequence_length_score", None)
                    value = (r * t * s) if (r is not None and t is not None and s is not None) else None
                else:
                    value = data.get(tool, {}).get(metric, None)

                if value is None:
                    summary[metric][tool] = {"value": None, "error": 0.0}
                else:
                    summary[metric][tool] = {"value": float(value), "error": 0.0}
        return summary

    raw_data = _looks_like_raw_data(metrics_data)

    if metric_names is None:
        selected_metrics = _collect_metrics(metrics_data, raw_data)
    elif isinstance(metric_names, str):
        selected_metrics = [metric_names]
    else:
        selected_metrics = list(metric_names)

    selected_metrics = [metric for metric in selected_metrics if metric not in excluded_metrics]

    if include_total_score and "total_score" not in selected_metrics and "total_score" not in excluded_metrics:
        selected_metrics.append("total_score")

    if len(selected_metrics) == 0:
        raise ValueError("No metrics remain after applying exclusions")

    selected_tools = _collect_tools(metrics_data, raw_data)
    if len(selected_tools) == 0:
        raise ValueError("No tools remain after applying exclusions")

    if raw_data:
        summary = _build_summary_from_raw(metrics_data, selected_metrics, selected_tools)
    else:
        summary = _build_summary_from_aggregated(metrics_data, selected_metrics, selected_tools)

    fig, ax = plt.subplots(figsize=figsize)
    x_positions = np.arange(len(selected_metrics))
    width = 0.8 / len(selected_tools)
    colors = plt.cm.Set3(np.linspace(0, 1, len(selected_tools)))

    for tool_index, tool in enumerate(selected_tools):
        values = []
        errors = []

        for metric in selected_metrics:
            entry = summary[metric][tool]
            value = entry["value"]
            if value is None:
                values.append(0.0)
                errors.append(0.0)
            else:
                values.append(value)
                errors.append(entry["error"])

        values = np.asarray(values, dtype=float)
        errors = np.asarray(errors, dtype=float)

        bar_x = x_positions + tool_index * width
        ax.bar(bar_x, values, width=width, label=tool, color=colors[tool_index], yerr=errors,
               capsize=3, ecolor='black', error_kw={"linewidth": 1})

    ax.set_xticks(x_positions + width * (len(selected_tools) - 1) / 2)
    ax.set_xticklabels([metric.replace('_', ' ').title() for metric in selected_metrics],
                       rotation=45, ha='right')
    ax.set_xlabel("Metric")
    ax.set_ylabel("Metric Value")
    ax.set_title("Metric Comparison Across Tools")
    ax.legend(title="Tools", loc="upper left", bbox_to_anchor=(1.02, 1.0))
    ax.grid(axis='y', alpha=0.3, linestyle='--')
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to {output_path}")

    if show_plot:
        plt.show()

    return fig, ax


def main():

    full_37_data = load_data("joined_data.json")
    paper_7_data = load_data("joined_data.json", 
                             included_fileNames = ["PID-1309-GAL-CA-1-PC_S105", "PID-1309-GAL-CA-2-PC_S106", "PID-1309-GAL-CA3-PC_S91",
                                                    "PID-1309-M7-MON-BSA-1_S87", "PID-1309-M7-MON-BSA-2_S88", "PID-1309-M7-MON-CONA-1_S85", "PID-1309-M7-MON-CONA-2_S86"])

    # runtime_comparison(full_37_data, "file_lengths.json", save_path="data/final_comparisons/full_37/runtime_comparison.png", show_plot=False)
    # runtime_comparison(paper_7_data, "file_lengths.json", save_path="data/final_comparisons/paper_7/runtime_comparison.png", show_plot=False)

    plot_metric_stacked_bar_with_error(full_37_data, metric_names=["retention_rate", "tau_score", "sequence_length_score"],
                                       excluded_tools=['forward', 'reverse', 'meta'], error_type="sem",
                                       output_path="data/final_comparisons/full_37/stacked_bar_comparison.png", show_plot=False)
    plot_metric_stacked_bar_with_error(paper_7_data, metric_names=["retention_rate", "tau_score", "sequence_length_score"],
                                       excluded_tools=['forward', 'reverse', 'meta'], error_type="sem",
                                       output_path="data/final_comparisons/paper_7/stacked_bar_comparison.png", show_plot=False)


    avgs = calc_averages(full_37_data, print_results=False)
    avgs["avg_rebollo"] = avg_rebollo(avgs)
    calculate_total_scores(avgs)
    # compare_table(avgs, include_ranking_sum=False,
    #               non_included_metrics=["upsilon_score", "phi_score"], 
    #               non_included_tools=["forward", "reverse", "rebollo_reverse", "rebollo"])
    
    # Example usage of plot_metric_barplot:
    # Plot time_seconds excluding rebollo tools
    # root_dir = "data/final_comparisons"
    # for data, name in zip([full_37_data, paper_7_data], ["full_37", "paper_7"]):
    #     print(f"Plotting time_taken for {name} dataset")
    #     save_dir = os.path.join(root_dir, name)
    #     os.makedirs(save_dir, exist_ok=True)
    #     for metric in ["time_taken", "retention_rate",  "sequence_length_score", "tau_score"]:
    #         output_path = os.path.join(save_dir, f"{name}_{metric}_comparison.png")

    #         _, _ = plot_metric_barplot(paper_7_data, metric_name=metric, 
    #                                     excluded_tools=['forward', 'reverse', 'meta'],
    #                                     output_path=output_path, show_plot=False)
            
    #         plt.cla()  # Clear the current axes for the next plot



if __name__ == "__main__":
    main()



