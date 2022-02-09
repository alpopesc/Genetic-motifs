

#include "treader.h"


TReader::TReader(std::vector<std::string> fasta_files ){
    for(size_t i(0); i < fasta_files.size(); ++i){
        motifs_vec.push_back(new TMotifs(fasta_files[i]));
    }
}

void TReader::read(TMotifs* element, std::string file) {
    char c('a');
    std::ifstream fin;
    fin.open(file, std::ios::in);

    if (!fin) {
        std::cerr << "error opening" << file << std::endl;
        throw std::string("File could not open. Please check if each file name corresponds to the name of the chromosome (Exemple Chromosome 1 -> seq1.fa");
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
        element -> setSequenceName(line);

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
            element -> setSequenceName(line);
            j = 0;
        }else if(c != '\n' and not isspace(c)){
            element->update(c, j);
            ++j;
        }
    }
    fin.close();
}

void TReader::init_reading(){


#define ASYNC 1
#if ASYNC

    for(auto& motifs : motifs_vec) {
        fut_vec.push_back(std::async(std::launch::async, TReader::read,motifs,motifs->getFileName()));
    }

    for(size_t i(0);i < fut_vec.size(); ++i){
        fut_vec[i].wait();
    }

#else
    for(auto& a : chromosomes) {
        read(a.first, a.second);
    }

#endif

}
