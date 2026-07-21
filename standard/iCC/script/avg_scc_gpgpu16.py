import os
import os.path
import sys
import csv

def print_head(out_file, method):
    fp = open(out_file, "w")
    if "hong" in method.lower():
        content = "Graph,Pivot,Trim,BFS,Update,WCC,Total Time"
        fp.writelines(content + "\n")
        fp.close()
    elif "slota" in method.lower():
        content = "Graph,Pivot,Trim,BFS,Update,P,B,Total Time"
        fp.writelines(content + "\n")
        fp.close()

def deal_graph_slota(file_name, graph, out_file):
#    graph_name = (file_name.split('result_')[1]).split('.')[0]
    time_list = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    #step_list = ['Pviot time', 'BFS time', 'Queue time', 'Color time', 'wcc computation time']
    step_list = ['trim_', 'pivot time', 'bfs time', 'update time', 'p time', 'b time', 'Time']

    count = 0;
    with open(file_name, "r") as fp:
        for line in csv.reader(fp):
            if len(line) > 0:
                for i in range(len(step_list)):
                    if line[0].startswith(step_list[i]):
                        if ":" in line[0]:
                            time_list[i] += float(line[0].strip().split(' ')[-1])
                            if i == len(step_list) - 1:
                                count += 1
                        else:
                            if len(line[1].strip().split(' ')) > 1:
                                time_list[i] += float(line[-1].strip().split(' ')[0])
                            else:
                                time_list[i] += float(line[1].strip())
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


def deal_graph_hong(file_name, graph, out_file):
#    graph_name = (file_name.split('result_')[1]).split('.')[0]
    time_list = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    #step_list = ['Pviot time', 'BFS time', 'Queue time', 'Color time', 'wcc computation time']
    step_list = ['trim_', 'pivot time', 'bfs time', 'update time', 'wcc time', 'Time']

    count = 0;
    with open(file_name, "r") as fp:
        for line in csv.reader(fp):
            if len(line) > 0:
                for i in range(len(step_list)):
                    if line[0].startswith(step_list[i]):
                        if ":" in line[0]:
                            time_list[i] += float(line[0].strip().split(' ')[-1])
                            if i == len(step_list) - 1:
                                count += 1
                        else:
                            if len(line[1].strip().split(' ')) > 1:
                                time_list[i] += float(line[-1].strip().split(' ')[0])
                            else:
                                time_list[i] += float(line[1].strip())
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

def deal_all(file_name, file_folder, method):
    fp = open(file_name, "r")
    data = fp.readlines()
    out_file = os.path.join(file_folder, "result_time_summary_{}.csv".format(method))
#    print data
    if len(data) < 1:
        print "No csv result files\n"
    else:
        print_head(out_file, method)
        for one in data:
            file_in = one.split(',')[0].strip()
            graph = one.split(',')[1].strip()
#            print file_in
            if "hong" in method.lower():
                deal_graph_hong(os.path.join(file_folder, file_in), graph, out_file)
            elif "slota" in method.lower():
                deal_graph_slota(os.path.join(file_folder, file_in), graph, out_file)



if __name__ == "__main__":
    if len(sys.argv) < 4:
        print "\nUsage: python avg.py <input_file_of_all> <file_folder> <vHong, wHong, vSlota, wSlota>\n"
        exit()

    print len(sys.argv), sys.argv[0], sys.argv[1], sys.argv[2], sys.argv[3]
#    cmd = "ls " + sys.argv[1] +  "result_*.csv > " + file_all;
#    os.system(cmd)
    deal_all(sys.argv[1], sys.argv[2], sys.argv[3])

