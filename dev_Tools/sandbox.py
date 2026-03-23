import json
import os
import csv

def load_metrics(metrics_path):

    data = {}

    with open(metrics_path, 'r') as f:
        data = json.load(f)

    parsed_data = {}

    for run_name, metrics in data.items():

        sanity_lambda = lambda x, default = 0: default if x is None else x
        
        parsed_data[run_name] = {
            "flash": {
                "retention_rate": sanity_lambda(data[run_name].get("flash", {}).get("total_reads", 0)) / sanity_lambda(data[run_name].get("flash", {}).get("total_pairs"), 1),
                "time_taken": sanity_lambda(data[run_name].get("flash", {}).get("time_seconds")),
                "tau_score": sanity_lambda(data[run_name].get("flash", {}).get("both_counts",0)) / sanity_lambda(data[run_name].get("flash", {}).get("total_reads"), 1),
            },
            "casper": {
                "retention_rate": sanity_lambda(data[run_name].get("casper", {}).get("total_reads", 0)) / sanity_lambda(data[run_name].get("casper", {}).get("total_pairs"), 1),
                "time_taken": sanity_lambda(data[run_name].get("casper", {}).get("time_seconds")),
                "tau_score": sanity_lambda(data[run_name].get("casper", {}).get("both_counts", 0))/ sanity_lambda(data[run_name].get("casper", {}).get("total_reads"),1)
            }
        }

    return parsed_data

def load_p3anut_evaluation(evaluation_path, p3anut_log_path):

    p3anut_log = load_p3anut_log(p3anut_log_path)

    data = {}

    with open(evaluation_path, 'r') as f:
        inp_data = json.load(f)
    

    for run_name, metrics in inp_data.items():
        file_name = os.path.basename(run_name).split("_R1_")[0]

        data[file_name] = {
            "forward": {
                "retention_rate": 1,
                "tau_score": inp_data[run_name].get("forward", {}).get("both_counts")/ inp_data[run_name].get("forward", {}).get("total_reads"),
            },
            "reverse": {
                "retention_rate": 1,
                "tau_score": inp_data[run_name].get("reverse", {}).get("both_counts")/ inp_data[run_name].get("reverse", {}).get("total_reads")
            },
            "p3anut": {
                "retention_rate": p3anut_log.get(file_name, {}).get("retention_rate"),
                "time_taken": p3anut_log.get(file_name, {}).get("time_taken"),
                "tau_score": inp_data[run_name].get("merged", {}).get("both_counts") / p3anut_log.get(file_name, {}).get("final_count")
            }
        }
    
    return data

def load_p3anut_log(log_path):

    data = []
    header = {}
    cleaned_data = {}

    with open(log_path, 'r') as f:
        data = csv.reader(f)
            
        t = next(data)
        for i, v in enumerate(t):
            header[v] = i



        for line in list(data)[20:]:
            file_name = os.path.basename(line[header["forwardFile"]]).split("_R1_")[0]

            cleaned_data[file_name] = {
                "retention_rate": float(line[header["finalCount"]]) / (float(line[header["finalCount"]]) + float(line[header["droppedCount"]])),
                "time_taken": float(line[header["totalTime"]]),
                "final_count": int(line[header["finalCount"]])
            }
    


    return cleaned_data

def join(p3anut_evaluation, flash_casper_metrics):

    joined_data = {}

    for run_name in p3anut_evaluation.keys():
        joined_data[run_name] = {
            "forward": p3anut_evaluation[run_name]["forward"],
            "reverse": p3anut_evaluation[run_name]["reverse"],
            "p3anut": p3anut_evaluation[run_name]["p3anut"],
            "flash": flash_casper_metrics.get(run_name, {}).get("flash", {}),
            "casper": flash_casper_metrics.get(run_name, {}).get("casper", {})
        }

    return joined_data


if __name__ == "__main__":

    p3anut_evaluation = load_p3anut_evaluation("p3anut_evaluation.json", "log.csv")
    flash_casper_metrics = load_metrics("metrics.json")

    joined_data = join(p3anut_evaluation, flash_casper_metrics)

    with open("joined_data.json", 'w') as f:
        json.dump(joined_data, f, indent=4)