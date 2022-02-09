#ifndef TEAM8_EXTENSION_H
#define TEAM8_EXTENSION_H
#include "pwm.h"


    /**
     * @class Optimized_Matrix
     * Optimized_Matrix is a subclass of PWM
     * this class will serve to find the optimized matrix by iterating between calculating the score of best motifs and calculating a new matrix
     */
class Optimized_Matrix : public PWM 
{

public :
	/**
     * @brief Constructor with parameters
     * @param motif_lenght nb of lines in our matrix
     * @param fa_file fasta file name with sequences to analyse
     * @param bed_name bed file name
     * @param multi true if there the sequence (fasta) comes in multiple files
     * @param threshold limit at which two matrices are close enough to be considered equal
     */
	Optimized_Matrix(unsigned int motif_length, vector <string> fa_files, string bed_name, double threshold);

	/**
     * @brief Function that iterates between calculating new matrixes and finding best score sequences to recalculate the matrix
     */
	void iterate();
	/**
     * @brief Verifies if two matrix are equal
     * @return true if the two matrix are equal (close to a certain threshhold)
     */
	bool equal_matrix();
	/**
     * @brief From a given sequence finds the given lenght sequence with best score
     * @param sequence the given sequence from which we need to find the best score sequence
     * @return the smaller, given size sequence with the best score
     */
	seq_t find_best_score_sequences(seq_t sequence);
	
	/**
     * @brief Setters
     * @param m matrix we want to set as attribute
     */
	void set_old_matrix(matrix_t m);
	void set_new_matrix(matrix_t m);

private :
	matrix_t old_matrix;
	matrix_t new_matrix;
	deque<seq_t> best_score_sequences;
	double threshold_;
};
#endif

