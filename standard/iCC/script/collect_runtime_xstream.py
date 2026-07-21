import os
import sys
import csv

time_domain = ["IO_WAIT", "SETUP", "WALL", "Compute"]

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
        time_list = [0.0, 0.0, 0.0, 0.0]
        if "xstream" in file_one:
            continue
        print file_one
        input_file = os.path.join(input_folder, file_one)

        with open(input_file, "r") as fp:
            for line in csv.reader(fp, delimiter = " "):
                if len(line) == 5:
                    if line[2].strip() == "CORE::TIME::IO_WAIT":
#                        print time_list[0] = float(line[-1].strip())
                        time_list[0] = float(line[-2].strip()) * 1000
                    if line[2].strip() == "CORE::TIME::SETUP":
                        time_list[1] = float(line[-2].strip()) * 1000
                    if line[2].strip() == "CORE::TIME::WALL":
                        time_list[2] = float(line[-2].strip()) * 1000

        time_list[3] = time_list[2] - time_list[0] - time_list[1]

#
#                    for index in range(len(time_domain)):
#                        if line[0].strip() == time_domain[index]:
##                            print line[0]
#                            time_list[index] = float(line[1])

#        for index in range(len(time_list)):
#            time_list[len(time_list) - 1] += time_list[index]

        result_list.append([file_one[:-4]] + time_list)

    return result_list

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print "Usage: python collect_runtime_Galois_async.py <input_folder>"
        exit(-1)
    input_folder = sys.argv[1]
    print("Output file galois_async_time.csv will be stored under ", input_folder)
    result_list = get_ours_time(input_folder)
    write_list_to_csv(result_list, os.path.join(input_folder, "xstream_time.csv"))
