#include "reader.h"
 

Reader::Reader(std::vector<std::string> fasta_files, std::string mtx_file, double threshold){
    for(size_t i(0); i < fasta_files.size(); ++i){
        motifs_vec.push_back(new Motifs(threshold,fasta_files[i],mtx_file));
    }
}

void Reader::read(Motifs* element, std::string file) {

    char c('a');
    std::ifstream fin;
    fin.open( file );

    if (!fin) {
        std::cerr << "error opening" << file << std::endl;
    }

    long int j(0);
    std::string line;
    fin.seekg (0, fin.end);
    int length = fin.tellg();
    fin.seekg(0,fin.beg);

    if(std::getline(fin ,line).good() and line[0] == '>') {
        fin.seekg(0, fin.beg);
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        line.erase(0,1); 
        element -> set_chromosome(line);

    }else{
        std::cerr << "Warning your fasta file does not contain any header" << "\n";
    }

    while (fin.good() and fin.tellg() < length) {
        fin.get(c);
        if(c == '>' and fin.good()){
            getline(fin, line);
            unsigned long int i = fin.tellg();
            fin.seekg(i - 1, fin.beg);
            fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            element -> set_chromosome(line);
            j = 0;
        }else if(c != '\n' and not isspace(c)) {

        unsigned int nuc_nb;
            switch(toupper(c)) {
                case ('A'):
                    nuc_nb = 0;
                    break;

                case ('C'):
                    nuc_nb = 1;
                    break;

                case ('G'):
                    nuc_nb = 2;
                    break;

                case ('T'):
                    nuc_nb = 3;
                    break;

                default:
                    nuc_nb = 4;
                    break;
            }
            element->update(nuc_nb, j);
            ++j;
        }
    }
    fin.close();

}

void Reader::init_reading(){


#define ASYNCRON 1
#if ASYNCRON

    for(auto& motifs : motifs_vec) {
        fut_vec.push_back(std::async(std::launch::async, Reader::read,motifs,motifs->getFilename()));
    }

    for(size_t i(0);i < fut_vec.size(); ++i){
        fut_vec[i].wait();
    }

#else
    for(size_t i(0); i < motifs_vec.size(); ++i) {
        read( motifs_vec[i], motifs_vec[i] ->getFilename()); 
    }

#endif

}
