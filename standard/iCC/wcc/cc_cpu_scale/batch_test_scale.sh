#alpha=30
##beta=200
#gamma=10
#theta=0.01
#async_limit=1
thread_num=24
run_times=10
result_folder=result_scale

#file[1]="livejournal"
#file[2]="flickr-growth"
#file[3]="wikipedia_en"
#file[4]="random"

file[1]="livejournal"
file[2]="baidu"
file[3]="pokec"
file[4]="wikipedia_link_en"
file[5]="wikipedia_en"
file[6]="RMAT"
file[7]="facebook"
file[8]="random"
file[9]="friendster"
file[10]="twitter_www"
file[11]="twitter_mpi"
#file[6]="dbpedia"

#file[3]="flickr-growth"
#file[4]="zhishi-hudong-relatedpages"
#file[10]="twitter_www"
#file[11]="wiki_talk_en"
#file[12]="wiki_communication"
#file[14]="us_patent"
#file[15]="citeseer"
#file[16]="scale25"
#file[10]="twitter_mpi"

deal () {

#    mkdir /mnt/raid0_huge/yuedeji/data/scc_result_scalability/${file[$1]}
#    ./scc_cpu /mnt/raid0_huge/yuedeji/data/${file[$1]}/fw_begin.bin /mnt/raid0_huge/yuedeji/data/${file[$1]}/fw_adjacent.bin /mnt/raid0_huge/yuedeji/data/${file[$1]}/bw_begin.bin /mnt/raid0_huge/yuedeji/data/${file[$1]}/bw_adjacent.bin $2 $alpha $beta $gamma $theta $run_times > /mnt/raid0_huge/yuedeji/scc_result/scc_time_result/${file[$1]}.csv
#    echo $2

    echo "./cc_cpu /home/yuedeji/data/data/${file[$1]}/undirected_begin.bin /home/yuedeji/data/data/${file[$1]}/undirected_adjacent.bin $2 $run_times $result_folder/scale_$2/histogram_${file[$1]}.csv > $result_folder/scale_$2/${file[$1]}.csv"
    ./cc_cpu /home/yuedeji/data/data/${file[$1]}/undirected_begin.bin /home/yuedeji/data/data/${file[$1]}/undirected_adjacent.bin $2 $run_times $result_folder/scale_$2/histogram_${file[$1]}.csv > $result_folder/scale_$2/${file[$1]}.csv

}

deal_scale () {
    export OMP_NUM_THREADS=$1
    for index in `seq 6 7`;
    do
        echo $index
        echo ${file[$index]}
        deal $index $1
    done
}

make clean
make

for scale in 1 2 4 8 16 32 64
do
    deal_scale $scale
done

