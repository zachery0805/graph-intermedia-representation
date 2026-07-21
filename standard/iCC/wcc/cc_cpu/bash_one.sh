#test_wikipedia_link_en:$(exe)
make clean
make

thread_num=64
#graph=twitter_mpi
#graph=random
graph=livejournal
#graph=friendster
#graph=wikipedia_link_en
#graph=pokec
#graph=RMAT
#graph=toy
#graph=facebook

export OMP_NUM_THREADS=$thread_num
#./cc_cpu ~/data/data/${graph}_large_cc/undirected_begin.bin ~/data/data/${graph}_large_cc/undirected_adjacent.bin $thread_num 1
./cc_cpu /home/cxl6029/data/${graph}/undirected_begin.bin /home/cxl6029/data/${graph}/undirected_adjacent.bin $thread_num 10 test.log

#./cc_cpu ~/data/data/friendster/undirected_begin.bin ~/data/data/friendster/undirected_adjacent.bin 24 1
#./cc_cpu ~/data/data/facebook/undirected_begin.bin ~/data/data/facebook/undirected_adjacent.bin 24 1
#./cc_cpu ~/data/data/livejournal/undirected_begin.bin ~/data/data/livejournal/undirected_adjacent.bin 24 1
#./cc_cpu ~/data/data/epinions/undirected_begin.bin ~/data/data/epinions/undirected_adjacent.bin 24 1
#./cc_cpu ~/data/data/toy/undirected_begin.bin ~/data/data/toy/undirected_adjacent.bin 56 1

