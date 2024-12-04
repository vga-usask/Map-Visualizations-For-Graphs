# Compilers
NVCC := nvcc
GXX := g++

# Directories
L_SRC_DIR := ../../lib
SRC_DIR := ../../src
OBJ_DIR := obj

# Debug compilation, disable by default
DEBUG ?= 0
ifeq ($(DEBUG), 1)
    NVCCFLAGS := -g 
    GXXFLAGS := -g
    O_LVL := -O0
else
    NVCCFLAGS := -O3 
    GXXFLAGS := -O3
endif

# C preprocessor flags
CPPFLAGS := -DNO_FREETYPE

# C/C++ compiler flags
GXXFLAGS := $(GXXFLAGS) -std=c++17 -lpthread `pkg-config --cflags cairo` -lboost_system -lboost_filesystem -lgvc -lcgraph 
NVCCFLAGS := $(NVCCFLAGS) -rdc=true --std=c++17 -lpthread `pkg-config --cflags cairo` -lboost_system -lboost_filesystem -lgvc -lcgraph 

# Linker flags - Add CGAL and potentially GMP and MPFR
LDFLAGS := -lc -lm -lstdc++ -lgmp -lmpfr `pkg-config --libs cairo` -lboost_system -lboost_filesystem -lgvc -lcgraph 

# Source and object files
GRAPH_VIEWER_SRCS := $(SRC_DIR)/graph_viewer.cpp
GRAPH_VIEWER_OBJS := $(OBJ_DIR)/graph_viewer.o

OTHER_SRCS := $(filter-out $(GRAPH_VIEWER_SRCS), $(wildcard $(SRC_DIR)/*.cpp))
OTHER_OBJS := $(OTHER_SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

PNGWRITER_SRCS := $(wildcard $(L_SRC_DIR)/*.cc)
PNGWRITER_OBJS := $(PNGWRITER_SRCS:$(L_SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)

CUDA_SRCS := $(wildcard $(SRC_DIR)/*.cu)
CUDA_OBJS := $(CUDA_SRCS:$(SRC_DIR)/%.cu=$(OBJ_DIR)/%.o)
CUDA_DEPS := $(wildcard $(SRC_DIR)/*.cuh)

OBJECTS := $(GRAPH_VIEWER_OBJS) $(OTHER_OBJS) $(PNGWRITER_OBJS) $(CUDA_OBJS)

depend: ddepend
ddepend: $(GRAPH_VIEWER_SRCS) $(OTHER_SRCS) $(PNGWRITER_SRCS) $(CUDA_SRCS)
	rm -f ddepend
	$(GXX) $(GXXFLAGS) $(CPPFLAGS) -M $(OTHER_SRCS) $(PNGWRITER_SRCS) >> ddepend
	$(NVCC) $(NVCCFLAGS) $(CPPFLAGS) -M $(GRAPH_VIEWER_SRCS) $(CUDA_SRCS) >> ddepend

include ./ddepend

graph_viewer: $(OBJECTS)
	$(NVCC) $(OBJECTS) $(LDFLAGS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	mkdir -p $(@D)
	$(GXX) -c $(CPPFLAGS) $(GXXFLAGS) -o $@ $<

$(OBJ_DIR)/graph_viewer.o: $(SRC_DIR)/graph_viewer.cpp
	mkdir -p $(@D)
	$(NVCC) -x cu -c $(CPPFLAGS) $(NVCCFLAGS) -o $@ $<

$(OBJ_DIR)/%.o: $(L_SRC_DIR)/%.cc
	mkdir -p $(@D)
	$(GXX) -c $(CPPFLAGS) $(GXXFLAGS) -o $@ $<

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cu
	mkdir -p $(@D)
	$(NVCC) -c $(CPPFLAGS) $(NVCCFLAGS) -o $@ $<

clear: clean

clean:
	rm -f graph_viewer $(OBJ_DIR)/*.o ddepend

.PHONY: all clear clean depend

