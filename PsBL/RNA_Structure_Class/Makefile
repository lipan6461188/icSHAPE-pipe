CC        = g++

CPPFLAGS    = -O3 -std=c++0x -Wall -D NDEBUG -pthread -pipe

TARGET_OBJ = RNA.o thermodynamics.o design.o RsampleData.o src/algorithm.o src/alltrace.o src/DynProgArray.o \
	src/dotarray.o src/draw.o src/extended_double.o src/forceclass.o src/MaxExpect.o src/MaxExpectStack.o src/outputconstraints.o \
	src/pfunction.o src/probknot.o src/random.o src/rna_library.o src/stackclass.o src/stackstruct.o src/stochastic.o src/structure.o \
	src/TProgressDialog.o src/common_utils.o HybridRNA.o TwoRNA.o src/bimol.o src/ParseCommandLine.o

libhybrid.a: $(TARGET_OBJ)
	ar rcs libhybrid.a $(TARGET_OBJ)

RNA.o: RNA.cpp
	$(CC) $(CPPFLAGS)  -c -o RNA.o RNA.cpp
thermodynamics.o: thermodynamics.cpp
	$(CC) $(CPPFLAGS)  -c -o thermodynamics.o thermodynamics.cpp
design.o: design.cpp
	$(CC) $(CPPFLAGS)  -c -o design.o design.cpp
RsampleData.o: RsampleData.cpp
	$(CC) $(CPPFLAGS)  -c -o RsampleData.o RsampleData.cpp
src/algorithm.o: src/algorithm.cpp
	$(CC) $(CPPFLAGS)  -c -o src/algorithm.o src/algorithm.cpp
src/alltrace.o: src/alltrace.cpp
	$(CC) $(CPPFLAGS)  -c -o src/alltrace.o src/alltrace.cpp
src/DynProgArray.o: src/DynProgArray.cpp
	$(CC) $(CPPFLAGS)  -c -o src/DynProgArray.o src/DynProgArray.cpp
src/dotarray.o: src/dotarray.cpp
	$(CC) $(CPPFLAGS)  -c -o src/dotarray.o src/dotarray.cpp
src/draw.o: src/draw.cpp
	$(CC) $(CPPFLAGS)  -c -o src/draw.o src/draw.cpp
src/extended_double.o: src/extended_double.cpp
	$(CC) $(CPPFLAGS)  -c -o src/extended_double.o src/extended_double.cpp
src/forceclass.o: src/forceclass.cpp
	$(CC) $(CPPFLAGS)  -c -o src/forceclass.o src/forceclass.cpp
src/MaxExpect.o: src/MaxExpect.cpp
	$(CC) $(CPPFLAGS)  -c -o src/MaxExpect.o src/MaxExpect.cpp
src/MaxExpectStack.o: src/MaxExpectStack.cpp
	$(CC) $(CPPFLAGS)  -c -o src/MaxExpectStack.o src/MaxExpectStack.cpp
src/outputconstraints.o: src/outputconstraints.cpp
	$(CC) $(CPPFLAGS)  -c -o src/outputconstraints.o src/outputconstraints.cpp
src/pfunction.o: src/pfunction.cpp
	$(CC) $(CPPFLAGS)  -c -o src/pfunction.o src/pfunction.cpp
src/probknot.o: src/probknot.cpp
	$(CC) $(CPPFLAGS)  -c -o src/probknot.o src/probknot.cpp
src/random.o: src/random.cpp
	$(CC) $(CPPFLAGS)  -c -o src/random.o src/random.cpp
src/rna_library.o: src/rna_library.cpp
	$(CC) $(CPPFLAGS)  -c -o src/rna_library.o src/rna_library.cpp
src/stackclass.o: src/stackclass.cpp
	$(CC) $(CPPFLAGS)  -c -o src/stackclass.o src/stackclass.cpp
src/stackstruct.o: src/stackstruct.cpp
	$(CC) $(CPPFLAGS)  -c -o src/stackstruct.o src/stackstruct.cpp
src/stochastic.o: src/stochastic.cpp
	$(CC) $(CPPFLAGS)  -c -o src/stochastic.o src/stochastic.cpp
src/structure.o: src/structure.cpp
	$(CC) $(CPPFLAGS)  -c -o src/structure.o src/structure.cpp
src/TProgressDialog.o: src/TProgressDialog.cpp
	$(CC) $(CPPFLAGS)  -c -o src/TProgressDialog.o src/TProgressDialog.cpp
src/common_utils.o: src/common_utils.cpp
	$(CC) $(CPPFLAGS)  -c -o src/common_utils.o src/common_utils.cpp
HybridRNA.o: HybridRNA.cpp
	$(CC) $(CPPFLAGS)  -c -o HybridRNA.o HybridRNA.cpp
TwoRNA.o: TwoRNA.cpp
	$(CC) $(CPPFLAGS)  -c -o TwoRNA.o TwoRNA.cpp
src/bimol.o: src/bimol.cpp
	$(CC) $(CPPFLAGS)  -c -o src/bimol.o src/bimol.cpp
src/ParseCommandLine.o: src/ParseCommandLine.cpp
	$(CC) $(CPPFLAGS)  -c -o src/ParseCommandLine.o src/ParseCommandLine.cpp


TARGET_OBJ_SMP = RNA-smp.o thermodynamics-smp.o design-smp.o RsampleData-smp.o src/algorithm-smp.o src/alltrace-smp src/DynProgArray-smp \
	src/dotarray-smp.o src/draw-smp.o src/extended_double-smp.o src/forceclass-smp.o src/MaxExpect-smp.o src/MaxExpectStack-smp.o src/outputconstraints-smp.o \
	src/pfunction-smp.o src/probknot-smp.o src/random-smp.o src/rna_library-smp.o src/stackclass-smp.o src/stackstruct-smp.o src/stochastic-smp.o src/structure-smp.o \
	src/TProgressDialog-smp.o src/common_utils-smp.o HybridRNA-smp.o TwoRNA-smp.o src/bimol-smp.o src/ParseCommandLine-smp.o

libhybrid-smp.a: $(TARGET_OBJ_SMP)
	ar rcs libhybrid-smp.a $(TARGET_OBJ_SMP)

RNA-smp.o: RNA.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o RNA-smp.o RNA.cpp
thermodynamics-smp.o: thermodynamics.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o thermodynamics-smp.o thermodynamics.cpp
design-smp.o: design.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o design-smp.o design.cpp
RsampleData-smp.o: RsampleData.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o RsampleData-smp.o RsampleData.cpp
src/algorithm-smp.o: src/algorithm.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o src/algorithm-smp.o src/algorithm.cpp
src/alltrace-smp.o: src/alltrace.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o src/alltrace-smp.o src/alltrace.cpp
src/DynProgArray-smp.o: src/DynProgArray.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o src/DynProgArray-smp.o src/DynProgArray.cpp
src/dotarray-smp.o: src/dotarray.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o src/dotarray-smp.o src/dotarray.cpp
src/draw-smp.o: src/draw.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o src/draw-smp.o src/draw.cpp
src/extended_double-smp.o: src/extended_double.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o src/extended_double-smp.o src/extended_double.cpp
src/forceclass-smp.o: src/forceclass.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o src/forceclass-smp.o src/forceclass.cpp
src/MaxExpect-smp.o: src/MaxExpect.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o src/MaxExpect-smp.o src/MaxExpect.cpp
src/MaxExpectStack-smp.o: src/MaxExpectStack.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o src/MaxExpectStack-smp.o src/MaxExpectStack.cpp
src/outputconstraints-smp.o: src/outputconstraints.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o src/outputconstraints-smp.o src/outputconstraints.cpp
src/pfunction-smp.o: src/pfunction.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o src/pfunction-smp.o src/pfunction.cpp
src/probknot-smp.o: src/probknot.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o src/probknot-smp.o src/probknot.cpp
src/random-smp.o: src/random.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o src/random-smp.o src/random.cpp
src/rna_library-smp.o: src/rna_library.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o src/rna_library-smp.o src/rna_library.cpp
src/stackclass-smp.o: src/stackclass.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o src/stackclass-smp.o src/stackclass.cpp
src/stackstruct-smp.o: src/stackstruct.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o src/stackstruct-smp.o src/stackstruct.cpp
src/stochastic-smp.o: src/stochastic.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o src/stochastic-smp.o src/stochastic.cpp
src/structure-smp.o: src/structure.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o src/structure-smp.o src/structure.cpp
src/TProgressDialog-smp.o: src/TProgressDialog.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o src/TProgressDialog-smp.o src/TProgressDialog.cpp
src/common_utils-smp.o: src/common_utils.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o src/common_utils-smp.o src/common_utils.cpp
HybridRNA-smp.o: HybridRNA.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o HybridRNA-smp.o HybridRNA.cpp
TwoRNA-smp.o: TwoRNA.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o TwoRNA-smp.o TwoRNA.cpp
src/bimol-smp.o: src/bimol.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o src/bimol-smp.o src/bimol.cpp
src/ParseCommandLine-smp.o: src/ParseCommandLine.cpp
	$(CC) $(CPPFLAGS) -D SMP -fopenmp -c -o src/ParseCommandLine-smp.o src/ParseCommandLine.cpp

clean:
	rm *.o *.a || true
	rm src/*.o || true



