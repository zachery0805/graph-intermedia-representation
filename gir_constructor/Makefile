# Compiler
CXX = g++
CXXFLAGS = -std=c++17 -O3 -fopenmp -I./includes -I./includes/scc_cpu

# Source files
SRCS = gir.cpp graph.cpp graph_undirected.cpp \
       includes/scc_cpu/scc_core.cpp \
       includes/scc_cpu/load_print.cpp

# Object files
OBJS = $(SRCS:.cpp=.o)

TARGET = gir

all: $(TARGET)

# Compile rule
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Link
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean
