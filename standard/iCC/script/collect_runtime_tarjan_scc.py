import os
import sys
import csv


def write_list_to_csv(result_list, output_file):

    with open(output_file, "w") as csv_file:
        csv_write = csv.writer(csv_file, delimiter = ',')
        for line in result_list:
            csv_write.writerow(line)


def get_scc_time(input_folder):

    file_list = os.listdir(input_folder)
    result_list = []
#    print(file_list)
    for file_one in file_list:
        if "time" in file_one or "histogram" in file_one:
            continue
        input_file = os.path.join(input_folder, file_one)
        time_multistep = 0.0

        time_scc = []
        with open(input_file, "r") as fp:
            for line in csv.reader(fp):
                if(len(line) > 1):
                    if("running_time" in line[0]):
            #            print line
                        time_scc.append(float(line[-1]))

        if len(time_scc) > 0:
            avg_time = sum(time_scc) / len(time_scc)
            result_list.append([file_one[:-4], avg_time])

    return result_list
#    print result_list

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print "Usage: python collect_runtime_multistep_cc.py <input_folder>"
        exit(-1)
    input_folder = sys.argv[1]
    print("Output file multistep_cc_time.csv will be stored under ", input_folder)
    result_list = get_scc_time(input_folder)
    write_list_to_csv(result_list, os.path.join(input_folder, "scc_tarjan_time.csv"))
