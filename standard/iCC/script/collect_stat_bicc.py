import os
import sys
import csv

stat_domain = ["num_ap",
        "bicc_number",
        "num_check_trim",
        "num_ap_trim",
        "num_bicc",
        "num_second_parent",
        "Queue"]

def write_list_to_csv(result_list, output_file):

    with open(output_file, "w") as csv_file:
        csv_write = csv.writer(csv_file, delimiter = ',')

        csv_write.writerow(["graph"] + stat_domain)

        for line in result_list:
            csv_write.writerow(line)


def get_ours_stat(input_folder):

    file_list = os.listdir(input_folder)
    result_list = []

#    print(file_list)
    for file_one in file_list:
        stat_list = [0, 0, 0, 0, 0, 0, 0]
        if file_one == "ours_stat.csv":
            continue
        input_file = os.path.join(input_folder, file_one)

        with open(input_file, "r") as fp:
            for line in csv.reader(fp):
                if len(line) > 1:
                    for index in range(len(stat_domain)):
                        if line[0].strip() == stat_domain[index] and stat_list[index] == 0:
                            #print line
                            if index != len(stat_domain) - 1:
                                stat_list[index] = int(line[1])
                            else:
                                stat_list[index] = int(line[-1])

#        for index in range(len(stat_list) - 1):
#            stat_list[len(stat_list) - 1] += stat_list[index]

        result_list.append([file_one[:-4]] + stat_list)

    return result_list

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print "Usage: python collect_runstat_dfs.py <input_folder>"
        exit(-1)
    input_folder = sys.argv[1]
    print("Output file ours_stat.csv will be stored under ", input_folder)
    result_list = get_ours_stat(input_folder)
    write_list_to_csv(result_list, os.path.join(input_folder, "ours_stat.csv"))
