
g++ -o test_htslib test_htslib.cpp ../../src/sam.cpp ../../src/htslib.cpp ../../src/string_split.cpp ../../src/fasta.cpp -lhts
./test_htslib ../test_data/input.bam 10

