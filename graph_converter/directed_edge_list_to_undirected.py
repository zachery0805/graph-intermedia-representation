import os
import sys
from sets import Set


def get_undirected_adj_list(input_file):

    undirected_edge_list = []
    edge_set = Set()
    is_head = True
    n = 0
    m = 0
    index = 0

    adj_dict = {}
    with open(input_file, 'r') as fp:
        for raw_line in fp:
            line_list = raw_line.strip().split()
            if is_head:
                n = int(line_list[0])
#                m = int(line_list[1])
                is_head = False
            else:
                u = int(line_list[0])
                v = int(line_list[1])

                if u not in adj_dict:
                    adj_dict[u] = {}
                if v not in adj_dict[u]:
                    adj_dict[u][v] = 1
                    m += 1


                if v not in adj_dict:
                    adj_dict[v] = {}

                if u not in adj_dict[v]:
                    adj_dict[v][u] = 1
                    m += 1

#                if (u, v) not in edge_set:
#                    undirected_edge_list.append((u, v))
#                    edge_set.add((u, v))
#
#                if (v, u) not in edge_set:
#                    undirected_edge_list.append((v, u))
#                    edge_set.add((v, u))
            index += 1
            if index % 1000000 == 0:
                print index

    return adj_dict, n, m
#    return undirected_edge_list, n

def dump_adj_dict_to_file(adj_dict, n, m, output_file):
    with open(output_file, "w") as fp:
        fp.write('%{} {}\n'.format(n, m))
        for v_adj in sorted(adj_dict.keys()):
            for w in sorted(adj_dict[v_adj].keys()):
                fp.write('{} {}\n'.format(v_adj, w))


    print("new edge list is saved to {}".format(output_file))

def directed_to_undirected(input_file, output_file):

    adj_dict, n, m = get_undirected_adj_list(input_file)

    dump_adj_dict_to_file(adj_dict, n, m, output_file)

#    fp_out = open(output_file, 'w')
#    is_head = True
#    with open(input_file, 'r') as fp:
#        for raw_line in fp:
#            line_list = raw_line.strip().split()
#            if is_head:
#                fp_out.write('%s %s %s\n' %(line_list[0], line_list[0], line_list[1]))
#                is_head = False
#            else:
#                fp_out.write('%s %s\n' %(int(line_list[0]) + 1, int(line_list[1]) + 1))
#    fp_out.close()

if __name__ == "__main__":

    if len(sys.argv) != 3:
        print("Usage: python <directed_edge_list_to_undirected.py> <edge_list_file> <output_file>\n")
        exit(-1)

    directed_to_undirected(sys.argv[1], sys.argv[2])

