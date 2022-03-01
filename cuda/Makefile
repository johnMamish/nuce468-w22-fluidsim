
include ../common.make

PROJ_BASE := ..
MODULE_BASE := .
# CUDA_LOC := /usr/local/cuda-11.0/targets/x86_64-linux

# Thanks to https://yuukidach.github.io/2019/08/05/makefile-learning/

TARGET = fluidsim-gpu

LNCC := nvcc
CC := g++
CUC := nvcc
CFLAGS := -g -Wall -Werror
CUFLAGS := -g -Xcompiler -Wall -Xcompiler -Werror -Werror all-warnings

OUT = ${PROJ_BASE}/out
BUILD = ${MODULE_BASE}/build
INCLUDE= ${MODULE_BASE}/include

MODULES := src src/sim src/io

INCS = $(wildcard *.h $(foreach fd, $(shell find ${INCLUDE} -type d), $(fd)/*.h))
LIB_INCS = 
C_SRCS = $(wildcard *.c $(foreach fd, $(MODULES), $(fd)/*.c))
CU_SRCS = $(wildcard *.c $(foreach fd, $(MODULES), $(fd)/*.cu))
NODIR_C_SRC = $(notdir $(C_SRCS))
NODIR_CU_SRC = $(notdir $(CU_SRCS))
OBJS = $(addprefix $(BUILD)/, $(C_SRCS:.c=_cc.o)) $(addprefix $(BUILD)/, $(CU_SRCS:.cu=_cu.o)) # obj/xxx.o obj/folder/xxx .o
LIBS =
LIB_DIRS =

CFLAGS += $(foreach incl, ${INCLUDE}, -I ${incl}) $(foreach incl, ${LIB_INCS}, -I ${incl})
CUFLAGS += $(foreach incl, ${INCLUDE}, -I ${incl}) $(foreach incl, ${LIB_INCS}, -I ${incl})

all: ${TARGET}

${BUILD}:
	-mkdir -p ${BUILD}

${OUT}:
	-mkdir -p ${OUT}

PHONY := $(TARGET)
$(TARGET): $(OBJS)
	echo "${TC_LINK} $(LNCC) -o $(OUT)/$@ $(OBJS) $(LIB_DIRS) $(LIBS)${TC_END}"
	$(LNCC) -o $(OUT)/$@ $(OBJS) $(LIB_DIRS) $(LIBS)

$(BUILD)/%_cc.o: %.c $(INCS)
	echo "${TC_BUILD}$(CC) -o $@ $(CFLAGS) -c $<${TC_END}"
	mkdir -p $(@D)
	$(CC) -o $@ $(CFLAGS) -c $<

$(BUILD)/%_cu.o: %.cu $(INCS)
	echo "${TC_BUILD_CU}$(CUC) -o $@ $(CUFLAGS) -c $<${TC_END}"
	mkdir -p $(@D)
	$(CUC) -o $@ $(CUFLAGS) -c $<

PHONY += clean
clean:
	echo "${TC_CLEAN} fluidsim-gpu${TC_END}"
	-rm -r ${BUILD} ${OUT}/${TARGET}

PHONY += echoes
echoes:
	@echo "INC files: $(INCS)"
	@echo "SRC files: $(C_SRCS)"
	@echo "OBJ files: $(OBJS)"
	@echo "LIB files: $(LIBS)"
	@echo "INC DIR: $(INCLUDE)"
	@echo "LIB DIR: $(LIB_DIRS)"

.PHONY: $(PHONY)