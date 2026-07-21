import os
import sys
import csv

time_domain = ["Trim",
        "BFS first",
        "Second parent",
        "Queue",
        "Bridge",
        "Total"]

def write_list_to_csv(result_list, output_file):

    with open(output_file, "w") as csv_file:
        csv_write = csv.writer(csv_file, delimiter = ',')

        csv_write.writerow(["graph"] + time_domain)

        for line in result_list:
            csv_write.writerow(line)


def get_ours_time(input_folder):

    file_list = os.listdir(input_folder)
    result_list = []

#    print(file_list)
    for file_one in file_list:
        time_list = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        if file_one == "ours_time.csv":
            continue
        input_file = os.path.join(input_folder, file_one)

        with open(input_file, "r") as fp:
            for line in csv.reader(fp):
                if len(line) > 1:
                    for index in range(len(time_domain) - 1):
                        if line[0] == time_domain[index]:
                            time_list[index] = float(line[1])
                            break

        for index in range(len(time_list) - 1):
            time_list[len(time_list) - 1] += time_list[index]

        result_list.append([file_one[:-4]] + time_list)

    return result_list

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print "Usage: python collect_runtime_dfs.py <input_folder>"
        exit(-1)
    input_folder = sys.argv[1]
    print("Output file ours_time.csv will be stored under ", input_folder)
    result_list = get_ours_time(input_folder)
    write_list_to_csv(result_list, os.path.join(input_folder, "ours_time.csv"))
