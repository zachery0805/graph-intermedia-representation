# Compiler
CXX = g++
CXXFLAGS = -std=c++17 -O3 -fopenmp -I./includes -I./includes/scc_cpu

# 可选编译开关
# 用法：
#   make DEBUG=1        # 开启DEBUG日志
#   make OUTPUT_TIME=1  # 开启时间输出
#   make DEBUG=1 OUTPUT_TIME=1 # 同时开启
ifeq ($(DEBUG),1)
    CXXFLAGS += -DDEBUG
endif
ifeq ($(OUTPUT_TIME),1)
    CXXFLAGS += -DOUTPUT_TIME
endif

# Source files
SRCS = main.cpp gir.cpp graph.cpp graph_undirected.cpp \
       includes/scc_cpu/scc_core.cpp \
       includes/scc_cpu/load_print.cpp \
       includes/wcc_cpu/wcc_core.cpp


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
