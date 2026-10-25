# Terrace driver. Needs Terrace built with CILK=1 (OpenCilk) in $(TERRACE).
TERRACE ?= /home/cxl6029/projects/terrace
CXX = /home/cxl6029/environment/opencilk/OpenCilk-12.0.0-Linux/bin/clang++
FLAGS = -std=c++17 -Ofast -DNDEBUG -march=native -DOPENMP=0 -DCILK=1 -fopencilk -I$(TERRACE) -I$(TERRACE)/include
OBJS = $(TERRACE)/obj/util.o $(TERRACE)/obj/PMA.o $(TERRACE)/obj/PMA_Lock.o $(TERRACE)/obj/partitioned_counter.o
LIBS = -lm -lpthread -lssl -lcrypto -fopencilk

all: terrace_driver terrace_sssp
terrace_driver: terrace_driver.cc ../common/protocol.h
	$(CXX) $(FLAGS) -o $@ $< $(OBJS) $(LIBS)
terrace_sssp: terrace_driver.cc ../common/protocol.h
	$(CXX) $(FLAGS) -DTERRACE_SSSP -o $@ $< $(OBJS) $(LIBS)
clean:
	rm -f terrace_driver terrace_sssp
