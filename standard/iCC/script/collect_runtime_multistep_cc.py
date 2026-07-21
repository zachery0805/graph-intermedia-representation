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
        if "time" in file_one or "histogram" in file_one:
            continue
        input_file = os.path.join(input_folder, file_one)
        time_multistep = 0.0

        with open(input_file, "r") as fp:
            for line in csv.reader(fp, delimiter = ":"):
                if(len(line) > 1):
                    if(line[0] == "TOTAL"):
#                    if(line[0] == "Time(s)"):
#                        print line[-1]
                        time_dfs = float(line[-1]) * 1000

        result_list.append([file_one[:-4], time_dfs])

    return result_list
#    print result_list

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print "Usage: python collect_runtime_multistep_cc.py <input_folder>"
        exit(-1)
    input_folder = sys.argv[1]
    print("Output file multistep_cc_time.csv will be stored under ", input_folder)
    result_list = get_dfs_time(input_folder)
    write_list_to_csv(result_list, os.path.join(input_folder, "multistep_cc_time.csv"))
