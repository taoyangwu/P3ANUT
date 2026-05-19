import json 
import os
import re


def analyse_file(file, total_count, forward = True):

    def _calculate_tau_score(seq):
        if forward:
            return "TCT" in seq and "GGTGGAGGT" in seq
        else:
            return "AGA" in seq and "ACCTCCACC" in seq
    sequence_lambda = lambda x: len(x) == 39

    with open(file, 'r') as f:
        lines = f.readlines()

    count = 0
    tau_count = 0
    sequence_count = 0

    for line in lines:
        splits = re.split(r'\s+', line.strip())
        count += int(splits[1])
        tau_count += int(splits[1]) if _calculate_tau_score(splits[2]) else 0
        sequence_count += int(splits[1]) if sequence_lambda(splits[2]) else 0


    return count/total_count, tau_count/count, sequence_count/count

def main(root_path, file_lengths, method = "max"):

    def _manage_method(current, new):
        if method == "max":
            return max(current, new)
        elif method == "min":
            return min(current, new)
        elif method == "rebollo_post_a":
            return current + new
        else:
            raise ValueError(f"Invalid method: {method}")
        
    def _calculate_average(data):
        forward_keys = set(data["rebollo_post_f"].keys())
        reverse_keys = set(data["rebollo_post_r"].keys())
        union_keys = forward_keys.union(reverse_keys)

        temp_data = {}

        for key in union_keys:
            forward_metrics = data["rebollo_post_f"].get(key, {})
            reverse_metrics = data["rebollo_post_r"].get(key, {})

            average_metrics = {}
            for metric in set(forward_metrics.keys()).union(set(reverse_metrics.keys())):
                forward_value = forward_metrics.get(metric)
                reverse_value = reverse_metrics.get(metric)

                if forward_value is not None and reverse_value is not None:
                    average_metrics[metric] = (forward_value + reverse_value) / 2
                elif forward_value is not None:
                    average_metrics[metric] = forward_value
                elif reverse_value is not None:
                    average_metrics[metric] = reverse_value

            temp_data[key] = average_metrics

        return temp_data


    main_data = {
        "rebollo_post_f": {},
        "rebollo_post_r": {},
        "rebollo_post_a": {},
    }


    for folder in os.listdir(root_path):

        run_name = re.split(r"_R\d_", folder)[0]

        foward_file = "_R1_" in folder

        run_count = 0

        run_retention_rate = 0
        run_tau_score = 0
        run_sequence_length_score = 0

        folder_path = os.path.join(root_path, folder)
        for folder_2 in os.listdir(folder_path):

            if "txt" in folder_2:
                continue

            folder_2_path = os.path.join(folder_path, folder_2)
            for file in os.listdir(folder_2_path):

                

                if "Translated_stats.txt" == file:
                    continue

                file_path = os.path.join(folder_2_path, file)

                #Skipping files if it is a directory
                if os.path.isdir(file_path):
                    continue

                if "union_length"  in file_lengths.get(run_name, {}):
                    total_count = file_lengths.get(run_name).get("union_length", 0)
                else:
                    total_count = file_lengths.get(run_name).get("total_length", 0)

                try:
                    run_count += 1
                    retention_rate, tau_score, sequence_length_score = analyse_file(file_path, total_count, forward = foward_file)

                    run_retention_rate = _manage_method(run_retention_rate, retention_rate)
                    run_tau_score = _manage_method(run_tau_score, tau_score)
                    run_sequence_length_score = _manage_method(run_sequence_length_score, sequence_length_score)

                    print(f"{file}: Retention Rate: {retention_rate}, Tau Score: {tau_score}, Sequence Length Score: {sequence_length_score}")
                except Exception as e:
                    print(f"Error processing {file_path}: {e}")
                    pass

        if run_count > 0:
            if method == "rebollo_post_a":
                run_retention_rate /= run_count
                run_tau_score /= run_count
                run_sequence_length_score /= run_count

            print(f"Run: {run_name}, Retention Rate: {run_retention_rate}, Tau Score: {run_tau_score}, Sequence Length Score: {run_sequence_length_score}")

        main_data["rebollo_post_f" if foward_file else "rebollo_post_r"][run_name] = {
            "retention_rate": run_retention_rate,
            "tau_score": run_tau_score,
            "sequence_length_score": run_sequence_length_score
        }

    main_data["rebollo_post_a"] = _calculate_average(main_data)

    with open("rebollo_postprocessing.json", 'w') as f:
        json.dump(main_data, f, indent=4)
        


if __name__ == "__main__":
    root_path = "/home/proxima/Desktop/Side_Projects/P3ANUT/data/Rebollo_Filtered"
    with open("joined_file_lengths.json", 'r') as f:
        file_lengths = json.load(f)

    main(root_path, file_lengths)

    