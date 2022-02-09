#ifndef TEAM8_PWM_H
#define TEAM8_PWM_H

#include <vector>
#include <array>
#include <list>
#include <utility>
#include <iostream>
#include <string>
#include <map>
#include <future>
#include <mutex>
#include <fstream>
#include <deque>
using namespace std;

static mutex m1;
static mutex m2;

enum Strand {PLUS, MINUS, BOTH};

//===================================================
//**************** Types definitions ****************
//===================================================
typedef struct { 
    string seq_name;
    unsigned long int beginning;
    unsigned long int end;
    Strand s;
} gen_region;

typedef vector<array<double, 4>> matrix_t;
typedef vector<int> seq_t;

//===================================================
//************** Some useful functions **************
//===================================================
/**
* @brief verify if a certain string is a number
* @param s string we want to make sure is a number
* @return true if s is a number
*/
bool is_number(const string& s);

/**
* @brief converts a nucleotide to a number corresponding to its column in the matrix
* @param c character we want to convert in an int
* @return the number in int
*/
int to_int(char c);


/**
* @brief getter to a sequence from lower to upper given positions
* @param lower start position of the sequence we want to extract
* @param upper end position of the sequence we want to extract
* @param current_pos current position called will correspond to the current position in the fasta file reading
* @return the sequence wanted
*/
seq_t getSequence(ifstream& is, const unsigned long int& lower, const unsigned long int& upper, unsigned long int& current_pos);

/**
* @brief converts a sequence to it's complementary sequence
* @param seq sequence of which we desire the complement
*/
seq_t convert_to_complement(const seq_t& seq);

//===================================================

/**
* @class PWM
* A Position Weight Matrix is a matrix corresponding to a genomic sequence.
* The columms represent the nucleotides, respectively A, T, G and C, while the lignes 
* represent the postion of each nucleotides in the sequence.
* The values in the matrix show the probability of finding a nucleotide in that particular position.
* This class will serve to compute the matrix of a given motif.
*/
class PWM {

private:
    multimap<pair<string, unsigned int>, gen_region> gen_regions;
    vector<string> fasta_files; 
    string bed_file;
    unsigned long int motif_counter;
    vector<future<void>> fut_vec1;

protected:
    deque<seq_t> sequences;
    matrix_t matrix;    


public:
    /**
     * @brief Constructor with parameters
     * @param motif_length motif lenght
     * @param fa_files fasta file names (string)
     * @param bed_name bed file name
     * @param multi true if there are multiple chromosomes in the file
     */
    PWM(unsigned int motif_length, vector<string> fa_files, string bed_name);


    /**
     * @brief Reads the BED file and initializes `bed_info`
     */
    void read_bed();


    /**
     * @brief Computes the pwm from the fasta file
     * @param sequence sequence from which we want to compute the matrix, comes as a vector of int to facilitate this step
     */
    void compute_matrix(seq_t sequence);

    /**
     * @brief Divides each matrix entry by the number of motifs to get the average we want
     * To call only once, after reading the whole fasta file
     */
    void divide_matrix(long double nb_of_motifs);


    /**
     * @brief Extracts all info from one BED line
     * @param line  The line to read
     * @return A `gen_region` initialized with info from the line
     */
    gen_region read_bed_line(string line);


    /**
     * @brief Finds all the sequences on a given chromosome and updates the vector with sequences;
     * @param the index of the chromosome on which to look for the sequences
     */
    static void find_sequences_on_chromosome(string chr,deque<seq_t>* sequences, const multimap<pair<string, unsigned int>, gen_region>* gen_regions ,
                                             unsigned int length);

    /**
     * @brief finds all given sequences and updates the vector with them.
     * (uses find_sequence_on_chromosome(int chr) etc.
     */
    void start_searching();

    /**
     * @brief produces the output matrix as pwm file
     */
    void matrix_output();

    /**
     * @brief add a sequence to our deque of sequences
     * @param sequence the sequence to add
     */
    void addSequence(const seq_t& sequence);


	/**
     * @brief method that puts all other methodes together and computes the wanted matrix 
     * and calls the output for the user to see the result
     */
    void init_pwm();
	
	/**
     * @brief adds a nucleotide to the count (column j ligne i)
     * @param j will represent a nucleotide 
     * @param i position at wich this nucleotide is found
     * @param weight number we want to add to the matrix for a given nucleotide and position to later
     * compute probability
     */
    void matrix_counter(int j, size_t i, double weight);
	
	/**
     * @brief getter to a sequence
     */
    deque<seq_t> getSeqeunces();
    
    /**
     * @brief getter to a matrix
     */
    matrix_t getMatrix();
    
    /**
     * @brief getter to a GenRegion
     */
    multimap <pair<string, unsigned int>, gen_region> getGenRegions();


    void start_pwm();
};
#endif
