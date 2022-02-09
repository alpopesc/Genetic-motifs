#ifndef TEAM8_MOTIFS_H
#define TEAM8_MOTIFS_H

#include <vector>
#include <deque>
#include <array>
#include <string>
#include <fstream>
#include <iostream>

typedef std::vector<std::array<double,4>> Matrix;

    /**
     * @struct Found_Motif
     * once a motif is found it will be saved in as found_motif to then be printed out
     * @param chromosome string that represents the chromosome the motif belongs to
     * @param position position of the forst nucleotide of the motif in the fasta file
     * @param coding shows if the motif is found on the coding or non-coding sequence
     * @param sequence string representing the motif
     * @param score score of the motif
     */
 
struct Found_Motif
{
	std::string chromosome;
	size_t position;
	char coding;
	std::string sequence;
	double score;
};

    /**
     * @class Motifs
     * a Motif is a set of nucleotids on a chromosome
     * this class will serve to find the motifs with highest score (affinity) from a PWM matrix and the genom fasta file
     */
 
class Motifs
{
	public:

	
    /**
     * @brief default constructor
     */
	Motifs();
	
	/**
     * @brief Constructor with parameters
     * @param  threshold maximum score threshold
     * @param  fa_filename Fasta file name
     * @param  mtx_filename tabulated file name
     */
	Motifs(double threshold, std::string fa_filename, std::string mtx_filename);
	
    /**
     * @brief Reads the tabulated file and initializes `matrix`
     * @param filename tabulated file name
     */
    void read_matrix(std::string filename);	
	
	 /**
     * @brief update function that can be called by reader that constently updates the called nuc and the position
     * @param c character representing the current read nucleotid in the fasta file
     * @param pos position of the current read nuc in the fasta file
     */
    void update(char c, unsigned long pos);
	
    /**
     * @brief Uptating the buffersequence so we can find a match and create a founf_motif
     * @param c character used to read fasta file char by char
     */
	void update_sequences(unsigned int c);

	/**
     * @brief Uptating the score for each buffersequence so we can know what sequence to save
     * @param no parameters
     * @return an int representing the column in the matrix
     */
    void update_score();
    
   	/**
     * @brief Once my_score is high enough our motif found on the coding strand is saved under found_motif and printed out
     */
    void decision_maker();
     

	void set_found_motif(char c);
	
	/**
     * @brief fonction that prints out the found motifs on a new tabulated file
     */

	void print_found_motifs();

	/**
     * @brief Transformes the deque of characters to string
     * @param d deque of character representing a sequence
     * @return a string equal to the deque of char
     */
	std::string transform_bufferseq_to_string(std::deque <unsigned int> d);
	
	/**
     * @brief Transformes the deque of characters representing a sequence to its complementary sequence
     * @param d deque of character representing a sequence
     * @return a deque of character 
     */
	unsigned int transform_to_complementary(unsigned int d);
	
	/**
     * @brief reinitialize the attributs bufferseq, cbufferseq, score and cscore of the motifs for the class motif's tests 
     * @param seq and cseq deque of unsigned int representing sequences, score and cscore double
     */
	void reinitialize(std::deque <unsigned int> seq, std::deque <unsigned int> cseq, double score, double cscore);
	
	 /**
     * @brief Getters of my attributs
     */
    std::deque <unsigned int> get_bufferseq();
    std::deque <unsigned int> get_cbufferseq();
	double get_score();
	double get_cscore();
	Found_Motif get_found_motif();
	Matrix get_matrix();
	
	/**
     * @brief Setters of my attributs
     */
	void set_length();
	void set_motif(std::vector<std::pair<std::deque <unsigned int>,size_t>>);
	void set_score(double score);
	void set_cscore(double cscore);
	void set_buffer_seq(std::deque <unsigned int> seq);
	void set_cbuffer_seq(std::deque <unsigned int> seq);
	void set_chromosome(std::string seq_name);
	std::string getFilename(){ return fasta_filename;};
	
    /**
     * @brief ressets the attribute my_buffersequence
     */	
	void reset_bufferseq();
	~Motifs();
	
	private:
	
		double score_threshold;
		std::string fasta_filename;
		std::deque <unsigned int> bufferseq;
		std::deque <unsigned int> cbufferseq;
		Matrix matrix;
		size_t position;
		std::string chromosome;
		Found_Motif found_motif;
		size_t length;
		double score;
		double score_c;
        std::ofstream outfile;

    std::string transform_bufferseq_to_string_comp(std::deque<unsigned int> d);
};

#endif
