#test_wikipedia_link_en:$(exe)
make clean
make

thread_num=24
#graph=random
graph=livejournal
#graph=friendster
#graph=wikipedia_link_en
#graph=pokec
#graph=random_large_cc
#graph=toy
#graph=facebook

./scc ~/data/data/$graph/undirected_begin.bin ~/data/data/$graph/undirected_adjacent.bin $thread_num 1

#./bicc ~/data/data/friendster/undirected_begin.bin ~/data/data/friendster/undirected_adjacent.bin 24 1
#./bicc ~/data/data/facebook/undirected_begin.bin ~/data/data/facebook/undirected_adjacent.bin 24 1
#./bicc ~/data/data/livejournal/undirected_begin.bin ~/data/data/livejournal/undirected_adjacent.bin 24 1
#./bicc ~/data/data/epinions/undirected_begin.bin ~/data/data/epinions/undirected_adjacent.bin 24 1
#./bicc ~/data/data/toy/undirected_begin.bin ~/data/data/toy/undirected_adjacent.bin 56 1

