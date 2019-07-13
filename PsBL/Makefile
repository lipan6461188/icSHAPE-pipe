CC        = g++

CPPFLAGS    = -O3 -std=c++0x -Wall -D NDEBUG -pthread -pipe
TARGET_DIR    = PsBL


primary: create_dir
	cd RNA_Structure_Class; make
	cp RNA_Structure_Class/*.h ${TARGET_DIR}/include/RNA_Structure_Class
	cp RNA_Structure_Class/src/*.h ${TARGET_DIR}/include/RNA_Structure_Class/src
	cp RNA_Structure_Class/libhybrid.a ${TARGET_DIR}/lib
	cd src; make
	cp src/*.h ${TARGET_DIR}/include
	cp src/libPsBL.a ${TARGET_DIR}/lib
	cd Bin_Src; make


create_dir:
	mkdir -p ${TARGET_DIR}
	mkdir -p ${TARGET_DIR}/include ${TARGET_DIR}/lib ${TARGET_DIR}/bin
	mkdir -p ${TARGET_DIR}/include/RNA_Structure_Class
	mkdir -p ${TARGET_DIR}/include/RNA_Structure_Class/src

install:
	cp Bin_Src/call_interaction ${TARGET_DIR}/bin
	cp Bin_Src/faformat ${TARGET_DIR}/bin
	cp Bin_Src/matrix_homo_summary ${TARGET_DIR}/bin
	cp Bin_Src/matrix_summary ${TARGET_DIR}/bin
	cp Bin_Src/paris_background ${TARGET_DIR}/bin
	cp Bin_Src/paris_prepare ${TARGET_DIR}/bin
	cp Bin_Src/sam2dg ${TARGET_DIR}/bin
	cp Bin_Src/sam2fq ${TARGET_DIR}/bin
	cp Bin_Src/sam2matrix ${TARGET_DIR}/bin
	cp Bin_Src/sam_group ${TARGET_DIR}/bin
	cp Bin_Src/sam_group_trim ${TARGET_DIR}/bin
	cp Bin_Src/sam_mismatch ${TARGET_DIR}/bin
	cp Bin_Src/sam_trim ${TARGET_DIR}/bin
	cp Bin_Src/sample_fold_params ${TARGET_DIR}/bin
	cp Bin_Src/sto2fa ${TARGET_DIR}/bin
	cp Bin_Src/dg_cluster ${TARGET_DIR}/bin

clean:
	rm -r ${TARGET_DIR} || true
	cd RNA_Structure_Class; make clean
	cd src; make clean
	cd Bin_Src; make clean

