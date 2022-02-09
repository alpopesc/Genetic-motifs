#ifndef TEAM8_READER_H
#define TEAM8_READER_H


#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <deque>
#include <future>
#include <array>
#include <mutex>
#include <pthread.h>
#include <thread>
#include <memory>
#include "motifs.h"




class Reader {
public:
    /**
    * @brief Constructor
    * @param names/paths of the fasta files, name/path of the matrix file
    */
    Reader(std::vector<std::string> fa_files, std::string mtx_file, double threshold); 

    /**
    * @brief Initiates the reading of one file and sends the information of the corresponding set of Motifs (class Motifs)
     * each set of motifs contains the name of the file where the motifs where found and the current sequence name
    * @param A pointer on the motif and the corresponding filename
    * @return void
    */
    static void read(Motifs* a, std::string file);

    /**
    * @brief Initiates the reading of the all the fasta files and the search for motifs
     * @param no parameters
    * @return void
    */
    void init_reading();

    Motifs* getMotifs(int i){ return motifs_vec[i];}

private:
/*! @name motif_vec
  This is a vector of motifs of which each motif contains a pwm and the location(filename) to read from
  As the given fasta file is read, the corresponding motif is updated as a new nucleotide is read.
*/
    std::vector<Motifs*> motifs_vec;
    std::vector<std::future<void>> fut_vec;

/*! @name Fasta files
  In case the file names are not given by the motifs and the user specifies to read from fasta files different from
  the chromosomal specified ones.
 */
    std::vector<std::string> fasta_files;


};


#endif 

