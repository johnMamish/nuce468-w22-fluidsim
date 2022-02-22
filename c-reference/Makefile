SOURCES:= main.c
#SOURCES+= ../external/cJSON/cJSON.c
SOURCES+= simulation.c

all: $(SOURCES)
	gcc -g -Ofast -ffast-math -Wall -std=gnu99 $(SOURCES) -lm -o fluidsim

clean:
	rm fluidsim

.PHONY: clean
