#include <gtest/gtest.h>
#include "../src/pwm.h"
#include <deque>


vector<string> pwm_fasta_files1({"../tests/files/mm9_fasta/chr1.fa", "../tests/files/mm9_fasta/chr2.fa",
                                  "../tests/files/mm9_fasta/chr7.fa",
                                  "../tests/files/mm9_fasta/chr8.fa",
                                  "../tests/files/mm9_fasta/chr11.fa",});
vector<string> pwm_fasta_files2({"../tests/files/promoters.fasta"});


string bed_file1("../tests/files/bed3cols.bed");
string bed_file2("../tests/files/bed6cols.bed");
string bed_file3("../tests/files/fasta_bed.bed");


PWM pwm(6,pwm_fasta_files2, bed_file1);
PWM pwmc(5,pwm_fasta_files2, bed_file1);


deque<seq_t>seq4({{3,0,1,1,3,2,2,1,0},{3,1,1,0,1,2,0},{1,3,3,2,2},{1,3,3,2,2},{1,3,0,2,1},{0,0,3,2,2},{3,1,3,2,2}});

vector<array<double,4>> mat({{0.18125, 0.539583, 0.03125, 0.247917},
                            { 0.222917, 0.289583, 0.05625, 0.43125},  
                            { 0.197917, 0.164583, 0.08125, 0.55625},  
                            { 0.0729167, 0.122917, 0.747917, 0.05625},   
                            { 0.0979167, 0.222917, 0.622917, 0.05625}});


bool operator==(const gen_region& a, const gen_region& b) {
    bool cmp = true;
    cmp &= (a.beginning == b.beginning);
    cmp &= (a.end == b.end);
    cmp &= (a.seq_name == b.seq_name);
    cmp &= (a.s == b.s);
    return cmp;
}

TEST(pwmTest, readBEDLineWorksWithValidLine) {
    gen_region expected1 = {"chr1", 1218600000, 1224500000, PLUS};
    gen_region result1 = pwm.read_bed_line("chr1\t1218600000\t1224500000\tregion1\t12.5\t+");
    gen_region expected2 = {"chr20", 50, 113, BOTH};
    gen_region result2 = pwm.read_bed_line("chr20\t50\t113");
    gen_region expected3 = {"chr300", 1300, 1345, MINUS};
    gen_region result3 = pwm.read_bed_line("chr300\t1300\t1345\tA\t12.5\t-");
    EXPECT_EQ(expected1, result1);
    EXPECT_EQ(expected2, result2);
    EXPECT_EQ(expected3, result3);
}

TEST(pwmTest, readBEDWorksWith3ColBED) {
    PWM pwm1(6,pwm_fasta_files1, bed_file1);
    gen_region expected1 = {"chr1", 120, 125, BOTH};
    gen_region expected2 = {"chr2", 109276, 110001,  BOTH};
    gen_region expected3 = {"chr1", 525, 621,  BOTH};
    vector<gen_region> exp;
    exp.push_back(expected1);
    exp.push_back(expected2);
    exp.push_back(expected3);

    multimap<pair<std::string, unsigned int>, gen_region> exp_map;
    exp_map.insert({{expected1.seq_name, expected1.beginning}, expected1});
    exp_map.insert({{expected2.seq_name, expected2.beginning}, expected2});
    exp_map.insert({{expected3.seq_name, expected3.beginning}, expected3});

    multimap<pair<std::string, unsigned int>, gen_region> res;
    res = pwm1.getGenRegions();
    for(auto it = res.begin(); it != end(res); it++){
    }
    EXPECT_EQ(exp_map, res);
}

TEST(pwmTest, readBEDWorksWith6ColBED) {
    PWM pwm2(6, pwm_fasta_files1, bed_file2);
    gen_region expected1 = {"chr3", 1, 12, BOTH};
    gen_region expected2 = {"chr1", 525, 621, MINUS};
    gen_region expected3 = {"chr2", 109276, 110001, MINUS};
    gen_region expected4 = {"chr1", 120, 125, PLUS};

    multimap<pair<string, unsigned int>, gen_region> exp_map;
    exp_map.insert({{expected1.seq_name, expected1.beginning}, expected1});
    exp_map.insert({{expected2.seq_name, expected2.beginning}, expected2});
    exp_map.insert({{expected3.seq_name, expected3.beginning}, expected3});
    exp_map.insert({{expected4.seq_name, expected4.beginning}, expected4});
    multimap<pair<string, unsigned int>, gen_region> resu;
    resu = pwm2.getGenRegions();
    EXPECT_EQ(exp_map, resu);
}

TEST(pwmTest, computeMatrix) {
    vector<array<double, 4>> matrix1 = {{1.25,0.25,0.25,1.25},
                                        {0.25,0.25,1.25,1.25},
                                        {1.25,1.25,0.25,0.25},
                                        {0.25,0.25,0.25,2.25}};
    vector<int> s1({0,2,1,3});
    vector<int> s2({3,3,0,3});
    PWM pwm1(4,pwm_fasta_files1, bed_file1);
    pwm1.compute_matrix(s1);
    pwm1.compute_matrix(s2);
    EXPECT_EQ(pwm1.getMatrix(),matrix1);
}

TEST (pwmTest, divide_matrix) {
    vector<array<double, 4>> matrix2 = {{0.625,0.125,0.125,0.625},
                                        {0.125,0.125,0.625,0.625},
                                        {0.625,0.625,0.125,0.125},
                                        {0.125,0.125,0.125,1.125}};
    vector<int> s1({0,2,1,3});
    vector<int> s2({3,3,0,3});
    PWM pwm1(4,pwm_fasta_files1, bed_file1);
    pwm1.compute_matrix(s1);
    pwm1.compute_matrix(s2);
    pwm1.divide_matrix(2.0);
    EXPECT_EQ(pwm1.getMatrix(), matrix2);
}


TEST(pwmTest, readFasta){
    PWM pwm5(5, pwm_fasta_files2, bed_file3);
    pwm5.start_searching();
    deque<seq_t> seqs1( { { 1, 1, 0, 0, 3, 1, 1, 0, 1, 1, 1 }, { 3, 0, 0, 0, 0, 1, 0, 0, 2, 0, 1, 1 }, { 2, 0, 0, 0, 0 } });
    EXPECT_EQ(seqs1, pwm5.getSeqeunces());
}

TEST(pwmTest, getSequenceTest){
    ifstream is;
    string chr;
    seq_t ret;
    unsigned long int counter(0);

    is.open("../tests/files/promoters.fasta");
    getline(is, chr);
    ret = getSequence(is, 3, 23, counter);
    seq_t exp({ 0,0,3,1,1,0,1,1,1,3,1,3,2,0,1,3,0,3,3,1});
    EXPECT_EQ(exp, ret);


    seq_t exp2({0,2,3,1,3,1,3,3,2,0,2,0,2});
    ret = getSequence(is,122,135, counter);
    EXPECT_EQ(exp2, ret);

    seq_t exp3 ({3,0,1,1});
    while(is.good()){
        if(getline(is, chr).good() and chr[0]=='>'){
            counter = 0;
            ret = getSequence(is, 100,104, counter);
        }
    }
    EXPECT_EQ(exp3, ret);
}

TEST(pwmTest , calc){
    for(size_t i(0); i< seq4.size(); ++i){
        pwmc.compute_matrix(seq4[i]);
    }
    pwmc.divide_matrix(8);
    for(size_t i(0); i < mat.size(); ++i )
	{
		for(size_t j(0); j < mat.size(); ++j )
		{
		EXPECT_NEAR(mat[i][j], pwmc.getMatrix()[i][j],0.00001);
		}
	}
}







