
g++ -o test_Fasta_class ../../src/fasta.cpp test_Fasta_class.cpp ../../src/string_split.cpp
./test_Fasta_class ../test_data/mm10_transcriptome.fa

g++ -o test_qFasta_class ../../src/fasta.cpp test_qFasta_class.cpp ../../src/string_split.cpp 
./test_qFasta_class ../test_data/mm10_transcriptome.fa 

