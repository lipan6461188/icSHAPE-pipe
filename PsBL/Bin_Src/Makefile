CC        = g++

TARGET_DIR    = PsBL-1.0.0

CPPFLAGS    = -O3 -std=c++0x -Wall -D NDEBUG -pthread -pipe -D NO_QT

HYBRIDINC = -I../RNA_Structure_Class -L../RNA_Structure_Class -lhybrid
PsBLINC = -I../src  -L../src -lPsBL -lhts

TARGET_OBJ = call_interaction faformat matrix_homo_summary matrix_summary paris_background \
	paris_prepare sam2dg sam2fq sam2matrix sam_group sam_group_trim sam_mismatch sam_trim \
	sample_fold_params sto2fa dg_cluster

PsBL: $(TARGET_OBJ)

call_interaction: call_interaction.cpp
	$(CC) call_interaction.cpp $(CPPFLAGS) $(PsBLINC)  -o call_interaction 
faformat: faformat.cpp
	$(CC) faformat.cpp $(CPPFLAGS) $(PsBLINC)  -o faformat 
matrix_homo_summary: matrix_homo_summary.cpp
	$(CC) matrix_homo_summary.cpp $(CPPFLAGS) $(PsBLINC)  -o matrix_homo_summary 
matrix_summary: matrix_summary.cpp
	$(CC) matrix_summary.cpp $(CPPFLAGS) $(PsBLINC)  -o matrix_summary 
paris_background: paris_background.cpp
	$(CC) paris_background.cpp $(CPPFLAGS) $(PsBLINC)  -o paris_background 
paris_prepare: paris_prepare.cpp
	$(CC) paris_prepare.cpp $(CPPFLAGS) $(PsBLINC)  -o paris_prepare 
sam2dg: sam2dg.cpp
	$(CC) sam2dg.cpp $(CPPFLAGS) $(PsBLINC)  -o sam2dg 
sam2fq: sam2fq.cpp
	$(CC) sam2fq.cpp $(CPPFLAGS) $(PsBLINC)  -o sam2fq 
sam2matrix: sam2matrix.cpp
	$(CC) sam2matrix.cpp $(CPPFLAGS) $(PsBLINC)  -o sam2matrix 
sam_group: sam_group.cpp
	$(CC) sam_group.cpp $(CPPFLAGS) $(PsBLINC)  -o sam_group 
sam_group_trim: sam_group_trim.cpp
	$(CC) sam_group_trim.cpp $(CPPFLAGS) $(PsBLINC)  -o sam_group_trim 
sam_mismatch: sam_mismatch.cpp
	$(CC) sam_mismatch.cpp $(CPPFLAGS) $(PsBLINC)  -o sam_mismatch 
sam_trim: sam_trim.cpp
	$(CC) sam_trim.cpp $(CPPFLAGS) $(PsBLINC)  -o sam_trim 
sample_fold_params: sample_fold_params.cpp
	$(CC) sample_fold_params.cpp $(CPPFLAGS) $(PsBLINC) $(HYBRIDINC)  -o sample_fold_params 
sto2fa: sto2fa.cpp
	$(CC) sto2fa.cpp $(CPPFLAGS) $(PsBLINC)  -o sto2fa 
dg_cluster: dg_cluster.cpp
	$(CC) dg_cluster.cpp $(CPPFLAGS) $(PsBLINC) $(HYBRIDINC)  -o dg_cluster 

clean:
	rm *.o *.a || true
	rm call_interaction || true
	rm faformat || true
	rm matrix_homo_summary || true
	rm matrix_summary || true
	rm paris_background || true
	rm paris_prepare || true
	rm sam2dg || true
	rm sam2fq || true
	rm sam2matrix || true
	rm sam_group || true
	rm sam_group_trim || true
	rm sam_mismatch || true
	rm sam_trim || true
	rm sample_fold_params || true
	rm sto2fa || true
	rm dg_cluster || true
