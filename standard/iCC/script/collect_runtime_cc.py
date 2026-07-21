import os
import sys
import csv


def write_list_to_csv(result_list, output_file):

    with open(output_file, "w") as csv_file:
        csv_write = csv.writer(csv_file, delimiter = ',')
        for line in result_list:
            csv_write.writerow(line)


def get_dfs_time(input_folder):

    file_list = os.listdir(input_folder)
    result_list = []
#    print(file_list)
    for file_one in file_list:
        if file_one == "cc_time.csv" or "histogram" in file_one:
            continue
        input_file = os.path.join(input_folder, file_one)
        time_dfs = 0.0
        num = 0

        with open(input_file, "r") as fp:
            for line in csv.reader(fp):
                if(len(line) > 1):
                    if(line[0] == "wcc computation time"):
#                    if(line[0] == "DFS Time"):
#                    if(line[0] == "Time(s)"):
                        time_dfs += float(line[-1][:-4])
                        num += 1
        if num != 0:
            result_list.append([file_one[:-4], time_dfs / num])

    return result_list
#    print result_list

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print "Usage: python collect_runtime_cc.py <input_folder>"
        exit(-1)
    input_folder = sys.argv[1]
    print("Output file cc_time.csv will be stored under ", input_folder)
    result_list = get_dfs_time(input_folder)
    write_list_to_csv(result_list, os.path.join(input_folder, "cc_time.csv"))
