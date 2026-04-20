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
                "upsilon_score": sanity_lambda(data[run_name].get("flash", {}).get("upsilon_score", 0)) / sanity_lambda(data[run_name].get("flash", {}).get("total_reads"), 1),
                "phi_score": sanity_lambda(data[run_name].get("flash", {}).get("phi_score", 0)) / sanity_lambda(data[run_name].get("flash", {}).get("total_reads"), 1),
                "sequence_length_score": sanity_lambda(data[run_name].get("flash", {}).get("sequence_length_score"))
            },
            "casper": {
                "retention_rate": sanity_lambda(data[run_name].get("casper", {}).get("total_reads", 0)) / sanity_lambda(data[run_name].get("casper", {}).get("total_pairs"), 1),
                "time_taken": sanity_lambda(data[run_name].get("casper", {}).get("time_seconds")),
                "tau_score": sanity_lambda(data[run_name].get("casper", {}).get("both_counts", 0))/ sanity_lambda(data[run_name].get("casper", {}).get("total_reads"),1),
                "upsilon_score": sanity_lambda(data[run_name].get("casper", {}).get("upsilon_score", 0)) / sanity_lambda(data[run_name].get("casper", {}).get("total_reads"), 1),
                "phi_score": sanity_lambda(data[run_name].get("casper", {}).get("phi_score", 0)) / sanity_lambda(data[run_name].get("casper", {}).get("total_reads"), 1),
                "sequence_length_score": sanity_lambda(data[run_name].get("casper", {}).get("sequence_length_score"))
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
                "time_taken": 0,
                "retention_rate": inp_data[run_name].get("forward", {}).get("retention_rate", 0),
                "tau_score": inp_data[run_name].get("forward", {}).get("both_counts")/ inp_data[run_name].get("forward", {}).get("total_reads"),
                "upsilon_score": inp_data[run_name].get("forward", {}).get("upsilon_score") / inp_data[run_name].get("forward", {}).get("total_reads"),
                "phi_score": inp_data[run_name].get("forward", {}).get("phi_score") / inp_data[run_name].get("forward", {}).get("total_reads"),
                "sequence_length_score": inp_data[run_name].get("forward", {}).get("sequence_length_score")
            },
            "reverse": {
                "time_taken": 0,
                "retention_rate": inp_data[run_name].get("reverse", {}).get("retention_rate", 0),
                "tau_score": inp_data[run_name].get("reverse", {}).get("both_counts")/ inp_data[run_name].get("reverse", {}).get("total_reads"),
                "upsilon_score": inp_data[run_name].get("reverse", {}).get("upsilon_score") / inp_data[run_name].get("reverse", {}).get("total_reads"),
                "phi_score": inp_data[run_name].get("reverse", {}).get("phi_score") / inp_data[run_name].get("reverse", {}).get("total_reads"),
                "sequence_length_score": inp_data[run_name].get("reverse", {}).get("sequence_length_score")
            },
            "p3anut": {
                "retention_rate": p3anut_log.get(file_name, {}).get("retention_rate"),
                "time_taken": p3anut_log.get(file_name, {}).get("time_taken"),
                "tau_score": inp_data[run_name].get("merged", {}).get("both_counts") / p3anut_log.get(file_name, {}).get("final_count"),
                "upsilon_score": inp_data[run_name].get("merged", {}).get("upsilon_score") / p3anut_log.get(file_name, {}).get("final_count"),
                "phi_score": inp_data[run_name].get("merged", {}).get("phi_score") / p3anut_log.get(file_name, {}).get("final_count"),
                "sequence_length_score": inp_data[run_name].get("merged", {}).get("sequence_length_score")
            }
        }
    
    return data

def load_rebollo(rebollo_path, sequence_lengths):

    data = {}
    forward_data = {}
    reverse_data = {}
    
    with open(rebollo_path, 'r') as f:
        data = json.load(f)

    run_name_lambda = lambda x, n = 1: os.path.basename(x).split(f"_R{n}_")[0]

    for raw in data:

        if raw["status"] != "ok":
           continue 

        if "R1" in raw["file_fastq"]:
            local_run_name = run_name_lambda(raw["file_fastq"])
            file_length = sequence_lengths.get(local_run_name, {}).get("forward_length", 1)
            pass
        elif "R2" in raw["file_fastq"]:
            local_run_name = run_name_lambda(raw["file_fastq"], n=2)
            file_length = sequence_lengths.get(local_run_name, {}).get("reverse_length", 1)
            pass

        log_entry = {
            "time_taken": raw.get("time_taken"),
            "retention_rate": raw.get("retention_rate") / file_length,
            "tau_score": raw.get("tau_score") / file_length,
            "upsilon_score": raw.get("upsilon_score") / file_length,
            "phi_score": raw.get("phi_score") / file_length,
            "sequence_length_score": raw.get("sequence_length_score") / file_length
        }

        if "R1" in raw["file_fastq"]:
            forward_data[local_run_name] = log_entry
        elif "R2" in raw["file_fastq"]:
            reverse_data[local_run_name] = log_entry

        

        


    return forward_data, reverse_data

def _load_lengths(sequence_length_path):

    sequence_lengths = {}
    with open(sequence_length_path, 'r') as f:
        sequence_lengths = json.load(f)

    return sequence_lengths


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
                "final_count": int(line[header["finalCount"]]),
            }
    


    return cleaned_data

def join(p3anut_evaluation, flash_casper_metrics, rebollo, rebollo_reverse, sequence_lengths):

    joined_data = {}

    for run_name in p3anut_evaluation.keys():
        joined_data[run_name] = {
            "forward": p3anut_evaluation[run_name]["forward"],
            "reverse": p3anut_evaluation[run_name]["reverse"],
            "p3anut": p3anut_evaluation[run_name]["p3anut"],
            "flash": flash_casper_metrics.get(run_name, {}).get("flash", {}),
            "casper": flash_casper_metrics.get(run_name, {}).get("casper", {}), 
            "rebollo": rebollo.get(run_name, {}),
            "rebollo_reverse": rebollo_reverse.get(run_name, {}),
            "meta": sequence_lengths.get(run_name, {})
        }

    return joined_data


if __name__ == "__main__":

    sequence_lengths = _load_lengths("file_lengths.json")
    p3anut_evaluation = load_p3anut_evaluation("p3anut_evaluation_ups_phi.json", "log.csv")
    flash_casper_metrics = load_metrics("metrics_ups_phi.json")
    rebollo_forward, rebollo_reverse = load_rebollo("rebollo_output_ups_phi_trimmed.json", sequence_lengths)

    joined_data = join(p3anut_evaluation, flash_casper_metrics, rebollo_forward, rebollo_reverse, sequence_lengths)

    with open("joined_data_r_trimmed.json", 'w') as f:
        json.dump(joined_data, f, indent=4)