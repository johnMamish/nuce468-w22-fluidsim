.SILENT:
TC_NOTIF=\033[38;5;207m
TC_CLEAN=\033[38;5;006m
TC_BUILD=\033[35;5;123mBUILD:
TC_RUN=\033[35;5;202mRUN:
TC_END=\033[0m

all: build/tmp/json-c

build/tmp/json-c:
	echo "  -> ${TC_BUILD} json-c${TC_END}"
	-rm -r build/tmp/json-c
	mkdir -p build/tmp/json-c
	cd build/tmp/json-c && cmake -DCMAKE_INSTALL_PREFIX=../.. -DBUILD_SHARED_LIBS=false ../../../json-c-0.15 && make -j16 && make install

clean:
	-rm -r build