
#include <gtest/gtest.h>
#include "treader.h"
#include <fstream>
#include <future>
#include <mutex>

std::vector<std::string> Files;
std::vector<std::string> reader_fasta_files1({"../tests/files/mm9_fasta/chr1.fa", "../tests/files/mm9_fasta/chr2.fa",
                                  "../tests/files/mm9_fasta/chr3.fa", "../tests/files/mm9_fasta/chr5.fa",
                                  "../tests/files/mm9_fasta/chr6.fa",
                                  "../tests/files/mm9_fasta/chr7.fa",
                                  "../tests/files/mm9_fasta/chr8.fa",
                                  "../tests/files/mm9_fasta/chr9.fa",
                                  "../tests/files/mm9_fasta/chr10.fa",
                                  "../tests/files/mm9_fasta/chr11.fa",
                                  "../tests/files/mm9_fasta/chr12.fa",
                                  "../tests/files/mm9_fasta/chr13.fa",
                                  "../tests/files/mm9_fasta/chr14.fa",
                                  "../tests/files/mm9_fasta/chr15.fa",
                                  "../tests/files/mm9_fasta/chr16.fa",
                                  "../tests/files/mm9_fasta/chr17.fa",
                                  "../tests/files/mm9_fasta/chr18.fa",
                                  "../tests/files/mm9_fasta/chr19.fa",
                                  "../tests/files/mm9_fasta/chrY.fa",
                                  "../tests/files/mm9_fasta/chrX.fa"});

std::string control_sequence1("TTCCCCAATAAAAAGACATAGACTAACAGACTGGCTACACAAACAGGACCCAACATTTTGCTGCTTAGAGGAAACCCATC");
std::string control_sequence2("TCCAATCCACCCTCTGACTATTCCACATCCCATACCTCCTCCTCATCCCCTTGTCTCCATGAGGATTTCTCCACCCCACC");
std::string control_sequence3("TCCAATCCACCCTCTGACTATTCCACATCCCATACCTCCTCCTCATCCCCTTGTCTCCATGAGGATTTCTCCACCCCACCTTCCCCAATAAAAAGACATAGACTAACAGACTGGCTACACAAACAGGACCCAACATTTTGCTGCTTAGAGGAAACCCATC");
static std::mutex read_mutex;

TEST(readingTest, parse) {
    for(size_t i(0); i < 4 ; ++i){
        Files.push_back("../tests/files/fasta_testfiles/seq"+ std::to_string(i+1)+ ".fa");
    }

    TReader reader(Files);
    reader.init_reading();


    EXPECT_EQ(79, reader.getMotifs(0)->getCounter());
    EXPECT_EQ(79, reader.getMotifs(1)->getCounter());
    EXPECT_EQ(79, reader.getMotifs(2)->getCounter());
    EXPECT_EQ(79, reader.getMotifs(3)->getCounter());

    EXPECT_EQ(control_sequence1, reader.getMotifs(0)->getEntireSequence());
    EXPECT_EQ(control_sequence1, reader.getMotifs(1)->getEntireSequence());
    EXPECT_EQ(control_sequence3, reader.getMotifs(2)->getEntireSequence());
    EXPECT_EQ(control_sequence2, reader.getMotifs(3)->getEntireSequence());

    EXPECT_EQ("chr11", reader.getMotifs(0)->getSequenceName());
    EXPECT_EQ("chr11", reader.getMotifs(1)->getSequenceName());
    EXPECT_EQ("chr11", reader.getMotifs(2)->getSequenceName());
    EXPECT_EQ("chr7", reader.getMotifs(3)->getSequenceName());
}





