#test_wikipedia_link_en:$(exe)
make clean
make

#thread_num=24
#graph=twitter_mpi
#graph=random
#graph=livejournal
#graph=friendster
#graph=wikipedia_link_en
#graph=pokec
graph=baidu
#graph=RMAT
#graph=toy
#graph=facebook


deal () {

    export OMP_NUM_THREADS=$1
#./cc_cpu ~/data/data/${graph}_large_cc/undirected_begin.bin ~/data/data/${graph}_large_cc/undirected_adjacent.bin $thread_num 1
    #./cc_cpu ~/data/data/${graph}/undirected_begin.bin ~/data/data/${graph}/undirected_adjacent.bin $thread_num 1 test.log
    echo ./cc_cpu ~/data/data/${graph}/undirected_begin.bin ~/data/data/${graph}/undirected_adjacent.bin $1 1 test.log
    ./cc_cpu ~/data/data/${graph}/undirected_begin.bin ~/data/data/${graph}/undirected_adjacent.bin $1 1 test.log
}
#./cc_cpu ~/data/data/friendster/undirected_begin.bin ~/data/data/friendster/undirected_adjacent.bin 24 1
#./cc_cpu ~/data/data/facebook/undirected_begin.bin ~/data/data/facebook/undirected_adjacent.bin 24 1
#./cc_cpu ~/data/data/livejournal/undirected_begin.bin ~/data/data/livejournal/undirected_adjacent.bin 24 1
#./cc_cpu ~/data/data/epinions/undirected_begin.bin ~/data/data/epinions/undirected_adjacent.bin 24 1
#./cc_cpu ~/data/data/toy/undirected_begin.bin ~/data/data/toy/undirected_adjacent.bin 56 1

for scale in 1 2 4 8 16 32 64
do
    deal $scale
done
