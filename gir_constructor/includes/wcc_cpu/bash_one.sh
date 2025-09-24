#test_wikipedia_link_en:$(exe)
make clean
make
# ./wcc_cpu /mnt/raid0_huge/yuede/data/livejournal/fw_begin.bin /mnt/raid0_huge/yuede/data/livejournal/fw_adjacent.bin /mnt/raid0_huge/yuede/data/livejournal/bw_begin.bin /mnt/raid0_huge/yuede/data/livejournal/bw_adjacent.bin 56 10 /home/yuede/gpu_cc/gpu_cc/iCC/wcc/result/result_wcc_livejournal.csv /home/yuede/gpu_cc/gpu_cc/iCC/wcc/result/result_time_livejournal.csv

#./wcc_cpu /mnt/raid0_huge/yuede/data/twitter_www/fw_begin.bin /mnt/raid0_huge/yuede/data/twitter_www/fw_adjacent.bin /mnt/raid0_huge/yuede/data/twitter_www/bw_begin.bin /mnt/raid0_huge/yuede/data/twitter_www/bw_adjacent.bin 56 10 /home/yuede/gpu_cc/gpu_cc/iCC/wcc/result/result_wcc_twitter_www.csv /home/yuede/gpu_cc/gpu_cc/iCC/wcc/result/result_time_twitter_www.csv

#./wcc_cpu /mnt/raid0_huge/yuede/data/zhishi-hudong-relatedpages/fw_begin.bin /mnt/raid0_huge/yuede/data/zhishi-hudong-relatedpages/fw_adjacent.bin /mnt/raid0_huge/yuede/data/zhishi-hudong-relatedpages/bw_begin.bin /mnt/raid0_huge/yuede/data/zhishi-hudong-relatedpages/bw_adjacent.bin 56 10 /home/yuede/gpu_cc/gpu_cc/iCC/wcc/result/result_wcc_zhishi-hudong-relatedpages.csv /home/yuede/gpu_cc/gpu_cc/iCC/wcc/result/result_time_zhishi-hudong-relatedpages.csv

./wcc_cpu /mnt/raid0_huge/yuede/data/baidu/fw_begin.bin /mnt/raid0_huge/yuede/data/baidu/fw_adjacent.bin /mnt/raid0_huge/yuede/data/baidu/bw_begin.bin /mnt/raid0_huge/yuede/data/baidu/bw_adjacent.bin 56 10 /home/yuede/gpu_cc/gpu_cc/iCC/wcc/result/result_wcc_baidu.csv /home/yuede/gpu_cc/gpu_cc/iCC/wcc/result/result_time_baidu.csv

