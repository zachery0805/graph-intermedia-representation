import os
import os.path
import sys
import csv

file_all = "file_name.txt"

def print_head(out_file):
    fp = open(out_file, "w")
    content = "Graph,Total Time"
    fp.writelines(content + "\n")
    fp.close()

def is_result_valid(input_file):

    with open(input_file, "r") as fp:
        for line in fp:
            if (line.strip()).startswith("avg. elapsed"):
                return True

    return False

def get_time(input_file):

    time_list = []
    with open(input_file, "r") as fp:
        for line in fp:
            if (line.strip()).startswith("avg. elapsed"):
                line_list = (line.strip()).split()
                run_time = float(line_list[-2])
                time_list.append(run_time)
    return time_list

def dump_list_to_csv(list_all, output_file):

    with open(output_file, "w") as fp:
        csv_write = csv.writer(fp)
        for line_one in list_all:
            csv_write.writerow(line_one)

    print "The list is saved in %s" %(output_file)

def deal_graph(file_name, graph, out_file):

#    graph_name = (file_name.split('result_')[1]).split('.')[0]
    time = 0.0
    #step_list = ['Pviot time', 'BFS time', 'Queue time', 'Color time', 'wcc computation time']
 #   step_list = ['WCC time']

    count = 0;
    with open(file_name, "r") as fp:
#        for line in csv.reader(fp, delimiter = ' '):
        for line in fp:
            if len(line) > 1 and (line.strip()).startswith("avg. elapsed"):
                line_list = (line.strip()).split()
#                print line_list
                time += float(line_list[-2])
                count += 1
#                time_list.append(run_time)

    print graph, count, time
    if count == 0:
        return
    content = graph
#    for i in range(len(time_list)):
    content += ',' + str(round(time / count, 3))

    with open(out_file, "a") as fp:
        fp.write(content + '\n')

def main(file_name, file_folder):
    fp = open(file_name, "r")
    data = fp.readlines()
    out_file = os.path.join(file_folder, "result_time_summary_gunrock.csv")
#    print data
    if len(data) < 1:
        print "No csv result files\n"
    else:
        print_head(out_file)
        for one in data:
            file_in = one.split(',')[0].strip()
            graph = one.split(',')[1].strip()

            deal_graph(os.path.join(file_folder, file_in), graph, out_file)

#    time_all = []
#    result_list = os.listdir(result_folder)
#    for result_one in result_list:
#        input_file = os.path.join(result_folder, result_one)
#        if result_one.endswith(".result") and is_result_valid(input_file):
#            time_list = [result_one]
#            time_list += get_time(input_file)
#            time_all.append(time_list)
#
#    dump_list_to_csv(time_all, output_file)

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "\nUsage: python summarize_result.py <input_file_of_all> <file_folder>\n"
        exit()

    print len(sys.argv), sys.argv[0], sys.argv[1], sys.argv[2]
    main(sys.argv[1], sys.argv[2])

