import json
import matplotlib.pyplot as plt

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


def calc_averages(metrics_path):

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


if __name__ == "__main__":
    calc_averages("joined_data.json")



