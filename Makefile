.SILENT:
.PHONY: c-reference cuda clean libs all

TC_NOTIF=\033[38;5;207m
TC_CLEAN=\033[38;5;006m
TC_BUILD=\033[35;5;123m
TC_RUN=\033[35;5;202mRUN:
TC_END=\033[0m

all: c-reference cuda

libs:
	echo "${TC_NOTIF}Building libraries${TC_END}"
	cd lib && make all

c-reference: libs
	echo "${TC_NOTIF}Building module 'c-reference'${TC_END}"
	cd c-reference && make all

cuda: libs
	echo "${TC_NOTIF}Building module 'cuda'${TC_END}"
	echo "not yet implemented"

test: test-c-reference test-cuda

test-cuda: cuda
	echo "${TC_NOTIF}Testing CUDA implementation${TC_END}"
	echo "not yet implemented"

test-c-reference: c-reference
	echo "${TC_NOTIF}Testing C implementation${TC_END}"
	-mkdir -p out/results
	echo "${TC_RUN} ./out/fluidsim ./out/results/sim.results${TC_END}"
	./out/fluidsim ./out/results/sim.results
	echo "${TC_RUN} python3 ./tools/render_fluid.py ./out/results/sim.results ./out/results/sim.avi${TC_END}"
	python3 ./tools/render_fluid.py ./out/results/sim.results ./out/results/sim.avi

clean:
	echo "${TC_CLEAN}Cleaning modules${TC_END}"
	cd c-reference && make clean
	-rm -r out

libs-clean:
	echo "${TC_CLEAN}Cleaning libraries${TC_END}"
	cd lib && make clean

full-clean: clean libs-clean