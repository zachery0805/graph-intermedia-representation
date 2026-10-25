# Aspen driver, built like Aspen's own tools (homegrown scheduler; NUM_THREADS=<p> at run time).
ASPEN ?= /home/cxl6029/projects/aspen/code
aspen_driver: aspen_driver.cpp ../common/protocol.h
	g++ -DEDGELONG -mcx16 -O3 -std=c++17 -march=native -Wno-subobject-linkage -DUSEMALLOC -DNDEBUG -DHOMEGROWN -pthread \
	    -I$(ASPEN) -o $@ $< -ldl
clean:
	rm -f aspen_driver
