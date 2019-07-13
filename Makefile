CXX        = g++

CXXFLAGS    = -O3 -std=c++0x -Wall -D NDEBUG -pthread -pipe
TARGET_DIR  = bin/Functions


primary:
	make -C icSHAPE
	make -C sliding_SHAPE

create_dir:
	mkdir -p ${TARGET_DIR}
	mkdir -p ${TARGET_DIR}/bin

install:
	cp sliding_SHAPE/sam2tab ${TARGET_DIR}
	cp sliding_SHAPE/calc_sliding_shape ${TARGET_DIR}
	cp sliding_SHAPE/countRT ${TARGET_DIR}

clean:
	rm ${TARGET_DIR}/sam2tab || true
	rm ${TARGET_DIR}/calc_sliding_shape || true
	rm ${TARGET_DIR}/countRT || true
	make -C icSHAPE clean
	make -C sliding_SHAPE clean

