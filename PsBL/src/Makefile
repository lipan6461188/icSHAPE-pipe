CC        = g++

TARGET_DIR    = PsBL-1.0.0

CPPFLAGS    = -O3 -std=c++0x -Wall -D NDEBUG -pthread -pipe

HYBRIDINC = -I../RNA_Structure_Class

TARGET_OBJ = align.o fasta.o fold.o pan_type.o param.o \
	paris_plot.o paris.o sam.o shape.o sstructure.o string_split.o htslib.o

libPsBL.a: $(TARGET_OBJ)
	ar rcs libPsBL.a $(TARGET_OBJ)

align.o: align.cpp
	$(CC) $(CPPFLAGS)  -c -o align.o align.cpp
fasta.o: fasta.cpp
	$(CC) $(CPPFLAGS)  -c -o fasta.o fasta.cpp
fold.o: fold.cpp
	$(CC) $(CPPFLAGS) $(HYBRIDINC) -c -o fold.o fold.cpp
pan_type.o: pan_type.cpp
	$(CC) $(CPPFLAGS)  -c -o pan_type.o pan_type.cpp
param.o: param.cpp
	$(CC) $(CPPFLAGS)  -c -o param.o param.cpp
paris_plot.o: paris_plot.cpp
	$(CC) $(CPPFLAGS) -D NO_QT -c -o paris_plot.o paris_plot.cpp
paris.o: paris.cpp
	$(CC) $(CPPFLAGS)  -c -o paris.o paris.cpp
sam.o: sam.cpp
	$(CC) $(CPPFLAGS)  -c -o sam.o sam.cpp
shape.o: shape.cpp
	$(CC) $(CPPFLAGS)  -c -o shape.o shape.cpp
sstructure.o: sstructure.cpp
	$(CC) $(CPPFLAGS)  -c -o sstructure.o sstructure.cpp
string_split.o: string_split.cpp
	$(CC) $(CPPFLAGS)  -c -o string_split.o string_split.cpp
htslib.o: htslib.cpp
	$(CC) $(CPPFLAGS)  -c -o htslib.o htslib.cpp


clean:
	rm *.o || true
	rm *.a || true
