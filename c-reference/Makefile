include ../common.make

SOURCES:= main.c
#SOURCES+= ../external/cJSON/cJSON.c
SOURCES+= simulation.c

INCLUDE := ../lib/build/include/json-c

all: ../out/fluidsim

../out/fluidsim: $(SOURCES)
	echo "${TC_BUILD}gcc -g -Ofast -ffast-math -Wall -std=gnu99 -I${INCLUDE} $(SOURCES) -lm -o ../out/fluidsim${TC_END}"
	-mkdir -p ../out
	gcc -g -Ofast -ffast-math -Wall -std=gnu99 -I${INCLUDE} $(SOURCES) -lm -o ../out/fluidsim

clean:
	-rm fluidsim

.PHONY: clean
