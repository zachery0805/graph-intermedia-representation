import os
import os.path
import sys
import csv

file_all = "file_name.txt"

def print_head(out_file):
    fp = open(out_file, "w")
    content = "Graph,Pivot,Trim,Large WCC,Graph Contraction,Small WCCs,Total Time"
    fp.writelines(content + "\n")
    fp.close()

def deal_graph(file_name, graph, out_file):

#    graph_name = (file_name.split('result_')[1]).split('.')[0]
    time_list = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    #step_list = ['Pviot time', 'BFS time', 'Queue time', 'Color time', 'wcc computation time']
    step_list = ['Pivot time', 'Trim time', 'BFS time', 'Queue time', 'Color time', 'WCC time']

    count = 0;
    with open(file_name, "r") as fp:
        for line in csv.reader(fp):
            if len(line) > 1:
                for i in range(len(step_list)):
                    if line[0].startswith(step_list[i]):
                        if len(line[1].strip().split(' ')) > 1:
                            time_list[i] += float(line[1].strip().split(' ')[0])
                        else:
                            time_list[i] += float(line[1].strip())
                        if i == 4:
                            count += 1
                        break

    print graph, count
    if count == 0:
        return
    content = graph
    for i in range(len(time_list)):
        content += ',' + str(round(time_list[i] / count, 3))

#    content += str(round(time_list[-1] / count, 3))

#    print content
    fp = open(out_file, "a")
    fp.writelines(content + "\n")
    fp.close()

def deal_all(file_name, file_folder):
    fp = open(file_name, "r")
    data = fp.readlines()
    out_file = os.path.join(file_folder, "result_time_summary_gpu_trim2.csv")
#    print data
    if len(data) < 1:
        print "No csv result files\n"
    else:
        print_head(out_file)
        for one in data:
            file_in = one.split(',')[0].strip()
            graph = one.split(',')[1].strip()
#            print file_in
            deal_graph(os.path.join(file_folder, file_in), graph, out_file)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "\nUsage: python avg.py <input_file_of_all> <file_folder>\n"
        exit()

    print len(sys.argv), sys.argv[0], sys.argv[1], sys.argv[2]
#    cmd = "ls " + sys.argv[1] +  "result_*.csv > " + file_all;
#    os.system(cmd)
    deal_all(sys.argv[1], sys.argv[2])

