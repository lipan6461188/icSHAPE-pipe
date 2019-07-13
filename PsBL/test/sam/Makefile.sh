
g++ -o test_read_bamorsam test_read_bamorsam.cpp ../../src/sam.cpp ../../src/htslib.cpp ../../src/string_split.cpp ../../src/fasta.cpp -lhts

./test_read_bamorsam ../test_data/input.sam
./test_read_bamorsam ../test_data/input.bam


