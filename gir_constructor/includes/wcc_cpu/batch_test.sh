alpha=30
#beta=200
gamma=10
theta=0.01
thread_count=56
run_times=20
async_limit=1
thread_num=56

#file[1]="livejournal"
#file[2]="flickr-growth"
#file[3]="wikipedia_en"
#file[4]="random"

file[1]="livejournal"
file[2]="baidu"
file[3]="flickr-growth"
file[4]="zhishi-hudong-relatedpages"
file[5]="pokec"
file[6]="wikipedia_link_en"
file[7]="wikipedia_en"
file[8]="dbpedia"
file[9]="facebook"
file[10]="twitter_www"
file[11]="wiki_talk_en"
file[12]="wiki_communication"
file[13]="random"
file[14]="RMAT"
file[15]="twitter_mpi"
file[16]="friendster"
#file[14]="us_patent"
#file[15]="citeseer"
#file[16]="scale25"

deal () {

#    mkdir /mnt/raid0_huge/yuede/data/scc_result_scalability/${file[$1]}
#    ./scc_cpu /mnt/raid0_huge/yuede/data/${file[$1]}/fw_begin.bin /mnt/raid0_huge/yuede/data/${file[$1]}/fw_adjacent.bin /mnt/raid0_huge/yuede/data/${file[$1]}/bw_begin.bin /mnt/raid0_huge/yuede/data/${file[$1]}/bw_adjacent.bin $2 $alpha $beta $gamma $theta $run_times > /mnt/raid0_huge/yuede/scc_result/scc_time_result/${file[$1]}.csv
#    echo $2
    echo "./wcc_cpu /mnt/raid0_huge/yuede/data/${file[$1]}/fw_begin.bin /mnt/raid0_huge/yuede/data/${file[$1]}/fw_adjacent.bin /mnt/raid0_huge/yuede/data/${file[$1]}/bw_begin.bin /mnt/raid0_huge/yuede/data/${file[$1]}/bw_adjacent.bin $thread_count $run_times /home/yuede/gpu_cc/gpu_cc/iCC/wcc/result/result_wcc_${file[$1]}.csv /home/yuede/gpu_cc/gpu_cc/iCC/wcc/result/result_time_${file[$1]}.csv"
    ./wcc_cpu /mnt/raid0_huge/yuede/data/${file[$1]}/fw_begin.bin /mnt/raid0_huge/yuede/data/${file[$1]}/fw_adjacent.bin /mnt/raid0_huge/yuede/data/${file[$1]}/bw_begin.bin /mnt/raid0_huge/yuede/data/${file[$1]}/bw_adjacent.bin $thread_count $run_times /home/yuede/gpu_cc/gpu_cc/iCC/wcc/result/result_wcc_${file[$1]}.csv /home/yuede/gpu_cc/gpu_cc/iCC/wcc/result/result_time_${file[$1]}.csv > /home/yuede/gpu_cc/gpu_cc/iCC/wcc/result/result_time_${file[$1]}.csv
#    ./scc_cpu /mnt/raid0_huge/yuede/data/${file[$1]}/fw_begin.bin /mnt/raid0_huge/yuede/data/${file[$1]}/fw_adjacent.bin /mnt/raid0_huge/yuede/data/${file[$1]}/bw_begin.bin /mnt/raid0_huge/yuede/data/${file[$1]}/bw_adjacent.bin $thread_count $alpha $2 $gamma $theta $run_times > ../result/test_beta/result_${file[$1]}_beta_${2}.csv

}
make clean
make
for index in `seq 1 16`;
do
    echo $index
    echo ${file[$index]}
    deal $index
#    sleep 2 
#    for alpha in 0 25 50 75 100 125 150 175 200 225 250 275 300;
#    do
#        echo $thread_num
#        deal $index
#        sleep 2 
#    done
#    for a in `seq 5 30`;
#    do
#        echo $a
#        deal $index $a
#    done
done
#./scc_cpu /mnt/raid0_huge/yuede/data/livejournal/fw_begin.bin /mnt/raid0_huge/yuede/data/livejournal/fw_adjacent.bin /mnt/raid0_huge/yuede/data/livejournal/bw_begin.bin /mnt/raid0_huge/yuede/data/livejournal/bw_adjacent.bin $alpha $beta $thread_count


