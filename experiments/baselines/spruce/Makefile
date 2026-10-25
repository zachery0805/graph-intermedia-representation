# Spruce driver. Needs Spruce's Release build (libspruce.a) and Junction.
SPRUCE ?= /home/cxl6029/projects/Spruce
JUNCTION_INC ?= /home/cxl6029/junction-install/include
JUNCTION_LIB ?= /home/cxl6029/projects/junction/build
spruce_driver: spruce_driver.cpp ../common/protocol.h
	g++ -O3 -DNDEBUG -std=gnu++20 -fopenmp -I$(SPRUCE)/src -I$(JUNCTION_INC) -o $@ $< \
	    $(SPRUCE)/build/libspruce.a $(JUNCTION_LIB)/libjunction.a $(JUNCTION_LIB)/turf/libturf.a -lpthread
clean:
	rm -f spruce_driver
