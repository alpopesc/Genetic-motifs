#include <gtest/gtest.h>
#include <math.h>
#include "../src/extension.h"

/*! FILES */
vector<string> fasta_file({"../tests/files/promoters.fasta"});
string bed_file_1("../tests/files/bed3cols.bed");

/*! TEST 1 */
TEST(extensionTest, equalMatrixTrue) {
    Optimized_Matrix matrix(4, fasta_file, bed_file_1, pow(10, -4.0));
    matrix_t old_m = {{1,0,0,1},
                    {0,0,1,1},
                    {1,1,0,0},
                    {0,0,0,2}};
    matrix_t new_m = {{1,0,0,1},
                    {0,0,1,1},
                    {1,1,0,0},
                    {0,0,0,2}};
    matrix.set_old_matrix(old_m);
    matrix.set_new_matrix(new_m);
    EXPECT_TRUE(matrix.equal_matrix());
}

/*! TEST 2 */
TEST(extensionTest, equalMatrixFalse) {
    Optimized_Matrix matrix(4, fasta_file, bed_file_1,  pow(10, -4.0));
    matrix_t old_m = {{2,0,0,1},
                    {0,0,1,1},
                    {1,1,0,0},
                    {0,0,0,2}};
    matrix_t new_m = {{1,0,0,1},
                    {0,0,1,1},
                    {1,1,0,0},
                    {0,0,0,2}};
    matrix.set_old_matrix(old_m);
    matrix.set_new_matrix(new_m);
    EXPECT_FALSE(matrix.equal_matrix());
    
    old_m = {{2,0,0,1},
            {0,0,1,1},
            {0,0,0,2}};
    new_m = {{1,0,0,1},
            {0,0,1,1},
            {1,1,0,0},
            {0,0,0,2}};
    matrix.set_old_matrix(old_m);
    matrix.set_new_matrix(new_m);
    EXPECT_FALSE(matrix.equal_matrix());
}

/*! TEST 3 */
TEST(extensionTest, findBestSequence) {
    Optimized_Matrix matrix(4, fasta_file, bed_file_1, pow(10, -4.0));
    matrix_t old_m = {{1,0,0,1},
                    {0,0,1,1},
                    {1,1,0,0},
                    {0,0,0,2}};
    matrix_t new_m = {{1,0,0,1},
                    {0,0,1,1},
                    {1,1,0,0},
                    {0,0,0,2}};
    matrix.set_old_matrix(old_m);
    matrix.set_new_matrix(new_m);
    seq_t seq = {1,0,1,1,0,2,0,1,1,0,1,0};
    seq_t res = {0,2,0,1};
    EXPECT_EQ(matrix.find_best_score_sequences(seq), res);
}