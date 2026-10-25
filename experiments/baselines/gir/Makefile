GIRC ?= ../../../gir_constructor
JEMALLOC ?= $(firstword $(wildcard /lib/x86_64-linux-gnu/libjemalloc.so.2 /usr/lib/x86_64-linux-gnu/libjemalloc.so.2))
CORE = $(GIRC)/gir.o $(GIRC)/gir_query.o $(GIRC)/graph.o $(GIRC)/graph_undirected.o

# The core objects are built by gir_constructor's own Makefile (its flags:
# -O3 -fopenmp); never by make's built-in rule from here.
gir_stream: gir_stream.cpp ../common/protocol.h core
	g++ -std=c++17 -O3 -fopenmp -I$(GIRC) -I$(GIRC)/includes -o $@ $< $(CORE) $(JEMALLOC)
core:
	$(MAKE) -C $(GIRC) gir.o gir_query.o graph.o graph_undirected.o
clean:
	rm -f gir_stream
.PHONY: core clean
