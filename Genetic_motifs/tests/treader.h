#ifndef TEAM8_TREADER_H
#define TEAM8_TREADER_H



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
#include "../tests/testClassMotifs.h"



class TReader {
public:
    /**
    * @brief Constructor
    * @param the index of the chromosome where to look for the sequences
    */
    TReader(std::vector<std::string> fa_files);

    /**
    * @brief Initiates the reading of one file and sends the information of the corresponding set of Motifs (class Motifs)
     * each set of motifs contains the name of the file where the motifs where found and the current sequence name
    * @param A pointer on the motif and the corresponding filename
    * @return void
    */
    static void read(TMotifs* a, std::string file);

    /**
    * @brief Initiates the reading of the all the fasta files and the search for motifs
     * @param no parameters
    * @return void
    */
    void init_reading();

    TMotifs* getMotifs(int i){ return motifs_vec[i];}

private:
/*! @name motif_vec
  This is a vector of motifs of which each motif contains a pwm and the location(filename) to read from
  As the given fasta file is read, the corresponding motif is updated as a new nucleotide is read.
*/
    std::vector<TMotifs*> motifs_vec;
    std::vector<std::future<void>> fut_vec;

/*! @name Fasta files
  In case the file names are not given by the motifs and the user specifies to read from fasta files different from
  the chromosomal specified ones.
 */
  

/*! @name Filenames
  This is a text files with all the names of the fasta files to read from.
 */
    std::string file_filenames;

};


#endif 

