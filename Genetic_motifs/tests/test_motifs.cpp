#include <gtest/gtest.h>
#include "../src/motifs.h"
#include "../src/reader.h"
#include <cmath>

/*! SEQUENCES */
std::deque <unsigned int> seq = {0, 3, 2};
std::deque <unsigned int> cseq = { 3, 0, 1};
std::deque <unsigned int> testSeq = {3, 2, 1};
std::deque <unsigned int> ctestSeq = {0, 1, 2};
std::string testseq = ("TGC");
std::string ctestseq = ("ACG");
std::string test_matrix = ("../tests/files/test_matrix.cvs");

/*! SCORES */
double score = 2+log2(0.1) +  2+log2(0.4) + 2+log2(0.3); 
double cscore = 2+log2(0.4) +  2+log2(0.1) + 2+log2(0.2); 


/*! CHARACTER */
char c('C');
char n('N');
char plus('+');
char minus('-');
char space(' ');
unsigned int one(1);
unsigned int five(5);


/*! MATRIX */
std::array<double, 4> vecteur_log = {2+log2(0.1),2+log2(0.2),2+log2(0.3),2+log2(0.4)};
Matrix matrix_log = { vecteur_log, vecteur_log, vecteur_log};

/*! MOTIF */
Motifs motif(cscore, "../tests/files/mm9_fasta/chr2.fa","../tests/files/test_matrix.cvs" );
Found_Motif found_motif = {" ", 0, '+', "ATG", score};
Found_Motif compl_motif = {" ", 0, '-', "CAT", cscore};

/*! FILES */
std::vector<std::string> fasta_files1({"../tests/files/mm9_fasta/chr2.fa"});
std::string mtx_file("../tests/files/test_matrix.cvs");
std::vector<std::string> fasta_files2({"../tests/files/mm9_fasta/chr1.fa", "../tests/files/mm9_fasta/chr2.fa",
                                  "../tests/files/mm9_fasta/chr3.fa", "../tests/files/mm9_fasta/chr5.fa",
                                  "../tests/files/mm9_fasta/chr6.fa",
                                  "../tests/files/mm9_fasta/chr7.fa",
                                  "../tests/files/mm9_fasta/chr8.fa",
                                  "../tests/files/mm9_fasta/chr9.fa"});

/*! TEST 1 */
	TEST(motifsTest, UpdateSequences) 
	{
		motif.reinitialize(seq, cseq, score, cscore);
		motif.update_sequences(one);
		EXPECT_EQ(motif.get_bufferseq(), testSeq);
		EXPECT_EQ(motif.get_cbufferseq(), ctestSeq);
	}



/*! TEST 2 */
	TEST(motifsTest, UpdateScore) 
	{
		motif.reinitialize(seq, cseq, score, cscore);
		motif.update_score();

		EXPECT_DOUBLE_EQ(motif.get_score(), 2+log2(0.1) + 2+log2(0.3) + 2+log2(0.4));
		EXPECT_DOUBLE_EQ(motif.get_cscore(), 2+log2(0.4) +  2+log2(0.1) + 2+log2(0.2));
	}


/*! TEST 3  */
	TEST(motifsTest, DecisionMaker) 
	{
		motif.reinitialize(seq, cseq, score, cscore);
		motif.decision_maker();
		EXPECT_EQ(motif.get_found_motif().score, found_motif.score);
		EXPECT_EQ(motif.get_found_motif().sequence, found_motif.sequence);
			
		motif.set_cscore(-100);
		motif.decision_maker();
		EXPECT_EQ(motif.get_found_motif().score, found_motif.score);
		EXPECT_EQ(motif.get_found_motif().sequence, found_motif.sequence); 
		
		motif.set_cscore(1);
		motif.set_score(-100);
		motif.decision_maker();
		EXPECT_EQ(motif.get_found_motif().score, 1);
		EXPECT_EQ(motif.get_found_motif().sequence, compl_motif.sequence);
	}
	
/*! TEST 4 */
	TEST(motifsTest, SetFoundMotifs) 
	{
		motif.reinitialize(seq, cseq, score, cscore);
		motif.set_found_motif(plus);
		EXPECT_EQ(motif.get_found_motif().chromosome, found_motif.chromosome);
		EXPECT_EQ(motif.get_found_motif().position, found_motif.position);
		EXPECT_EQ(motif.get_found_motif().coding, found_motif.coding);
		
		motif.set_found_motif(n);	
	} 


/*! TEST 5 */
	TEST(motifsTest, transformBufferseqToString) 
	{
		EXPECT_EQ(motif.transform_bufferseq_to_string(seq), "ATG");

	}

/*! TEST 6 */
	TEST(motifsTest, transformToComplementary) 
	{
		EXPECT_EQ(motif.transform_to_complementary(one), 2);
		EXPECT_EQ(motif.transform_to_complementary(five), 4);
	}

/*! TEST 7 	*/
	TEST(motifsTest, ReadMatrix) 
	{

		motif.reinitialize(seq, cseq, score, cscore);
		for(size_t i(0); i < matrix_log.size(); ++i )
		{
			EXPECT_DOUBLE_EQ(motif.get_matrix()[i][seq[i]], matrix_log[i][seq[i]]); 
		}

	}

	

