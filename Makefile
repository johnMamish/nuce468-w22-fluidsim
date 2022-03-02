
include common.make

.PHONY: c-reference cuda clean libs all

all: c-reference cuda

libs:
	echo "${TC_NOTIF}Building libraries${TC_END}"
	make all -C lib

c-reference: libs
	echo "${TC_NOTIF}Building module 'c-reference'${TC_END}"
	make all -C c-reference

cuda: libs
	echo "${TC_NOTIF}Building module 'cuda'${TC_END}"
	make all -C cuda

test: test-c-reference test-cuda

test-cuda: cuda
	echo "${TC_NOTIF}Testing CUDA implementation${TC_END}"
	-mkdir -p out/results
	echo "${TC_RUN} ./out/fluidsim-gpu -o ./out/results/sim-gpu.results${TC_END}"
	./out/fluidsim-gpu -o ./out/results/sim-gpu.results -l 120 -d "(400,200)" -b "(0.1,0)"
	echo "${TC_RUN} python3 ./tools/render_fluid.py ./out/results/sim-gpu.results ./out/results/sim-gpu.avi${TC_END}"
	python3 ./tools/render_fluid.py ./out/results/sim-gpu.results ./out/results/sim-gpu.avi

test-c-reference: c-reference
	echo "${TC_NOTIF}Testing C implementation${TC_END}"
	-mkdir -p out/results
	echo "${TC_RUN} ./out/fluidsim ./out/results/sim-cpu.results${TC_END}"
	./out/fluidsim ./out/results/sim-cpu.results
	echo "${TC_RUN} python3 ./tools/render_fluid.py ./out/results/sim-cpu.results ./out/results/sim-cpu.avi${TC_END}"
	python3 ./tools/render_fluid.py ./out/results/sim-cpu.results ./out/results/sim-cpu.avi

clean:
	echo "${TC_CLEAN}Cleaning modules${TC_END}"
	make clean -C c-reference
	make clean -C cuda
	-rm -r out

libs-clean:
	echo "${TC_CLEAN}Cleaning libraries${TC_END}"
	cd lib && make clean

full-clean: clean libs-clean