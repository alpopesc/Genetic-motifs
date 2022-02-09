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
#include "pwm.h"
#include "reader.h"
#include "../tests/treader.h"
#include <cmath>
#include "Timer.h"
#include "extension.h"

/*
class Timer{

public:
    Timer(){
        m_StartTimepoint = std::chrono::high_resolution_clock::now();
    }

    ~Timer(){
        Stop();
    }

    void Stop() {
        auto endTimepoint = std::chrono::high_resolution_clock::now();
        auto start = std::chrono::time_point_cast<std::chrono::microseconds>(m_StartTimepoint).time_since_epoch();
        auto end = std::chrono::time_point_cast<std::chrono::microseconds>(endTimepoint).time_since_epoch();

        auto duration = end - start;
        double ms = duration.count()*0.000001;

        std::cout << ms << "s";
    }

private:
    std::chrono::time_point<std::chrono::high_resolution_clock > m_StartTimepoint;

};
*/

static std::mutex read_mutex;

static void some_work(char c){
    std::cout << c;
}

class Queuetest {
public:
    void update(char c) {
        //std::cout<< buffer_sequence.size() << '\n';
        if (buffer_sequence.size() < 7) {
            buffer_sequence.push_back(c);
        } else {
            //print();
            buffer_sequence.pop_front();
            buffer_sequence.push_back(c);
            //std::cout << c << "\n";
        }
    }

    void print(){
        for (std::size_t i(0); i < buffer_sequence.size(); ++i) {
            std::cout << buffer_sequence[i];
        }
        std::cout << std:: endl;
    }

private:
    std::deque<char> buffer_sequence;
};

/*
static void read(Queuetest* a, std::string file) {
    char c('a');
    std::ifstream fin;
    fin.open(file, std::ios::in);

    if (!fin) {
        std::cerr << "error opening" << file << std::endl;
        throw -1;
    }


    long int j(0);
    bool sp(false);
    std::string line;
    */
    /*
    if(std::getline(fin ,line).good() and line[0] == '>') {
        fin.seekg(0, fin.beg);
        fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }*/
    /*
    while (!fin.eof()) {

        fin.get(c);
        a -> update(c);
        ++j;
    }

    fin.close();
}

int main() {
    //auto start = std::chrono::high_resolution_clock::now();
    Timer b;


    std::vector<std::future<void>> vec;
    std::vector<std::pair<Queuetest *, std::string>> chromosomes;
    for (int i(0); i < 4; i++) {
        chromosomes.push_back(
                std::pair<Queuetest *, std::string>(new Queuetest, "mm9_fasta/chr" + std::to_string(i + 1) + ".fa"));
    }

    //chromosomes.push_back(std::pair<Queuetest *, std::string>(new Queuetest, "mm9_fasta/chrX.fa"));
    //chromosomes.push_back(std::pair<Queuetest *, std::string>(new Queuetest, "mm9_fasta/chrY.fa"));


    for (auto &a : chromosomes) {
        read(a.first, a.second);
    }
*/

/*
#define ASYNC 1
#if ASYNC
    for(auto& a : chromosomes) {
        vec.push_back(std::async(std::launch::async, read,a.first, a.second));
    }
#else
    for(auto& a : chromosomes){
        read(a.first, a.second);
        }
    }
#endif





    return 0;
}
*/
/*
int main(){
    Timer t;
    std::ifstream fin("baba.fna", std::ios::in);
    char c;
    Queuetest q;
    while(!fin.eof()){
        c = fin.get();
        q.update(c);

    }
}
*/
/*So reading the fasta file and doing part of the update operation is extremely time intensive.
 * It is however the fastest way to find the motifs (maybe... Im not entirely sure). So For finding the
 * motifs a multi thread implementation can come in quite handy. At the moment it takes about 0.09 seconds to read the
 * dystrophin gene (~3mb) and 90(3.4gb) seconds for the human genome. With a multi thread implementation we could speed up this
 * process by a factor of 4-8 depending on how many cores are available.
 */
/*
int main() {
    Queuetest a;
    auto start = std::chrono::high_resolution_clock::now();
    char c('a');
    string line;
    vector<char> res;
    ifstream fin;
    //fin.open("baba.fna", std::ifstream::binary);
    fin.open("baba.fna", ios::in);
    if (!fin) {
        cerr << "error" << endl;
        return -1;
    }
    fin.seekg (0, fin.end);
    long int length = fin.tellg();
    for(int i(0); i < 2; ++i){
        long int x = 0.5*i;
        long int y = 0.5 + 0.5*i;
        string file;
        async(launch::async, read, x, y, file);
    }
    long int i(0);
    bool sp(false);
    //fin.seekg(0, fin.beg);
    while (!fin.eof()) {
        if(!sp) {
            fin.seekg(0, fin.beg);
            sp = true;
        }
        fin.get(c);
        //some_work(c);
        a.update(c);
        ++i;
    }
    cout << i << endl;
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> duration = end - start;
    std::cout << duration.count() << "s" <<  std::endl;
    return 0;
}
*/

/*
 * Reading a fasta file this way is very slow if we intend to find motifs. However it turns of to be very fast when
 * searching certain regions in a genome given the cooridantes of that genome. This can turn out useful when calculating
 * the pwm
 */

/*
int main()
{
    auto start = std::chrono::high_resolution_clock::now();
    std::ifstream is ("tata.fa", std::ifstream::binary);
    if (is) {
        // get length of file:
        is.seekg (0, is.end);
        int length = is.tellg();
        // allocate memory:
        char *buffer = new char[1];
        for(int i(0); i < length; ++i) {
            is.seekg(i, is.beg);
            // read data as a block:
            is.read(buffer, 1);
            // print content:
            //std::cout.write (buffer,length);
            //cout << buffer << endl;
        }
            is.close();
        // print content:
        //std::cout.write (buffer,length);
        //cout << buffer << endl;
        //delete[] buffer;
        }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> duration = end - start;
    std::cout << duration.count() << "s" <<  std::endl;
        return 0;
}
*/


/*Top-Down approach: sequence is known The down side is that that we cant do good error handling
 * it is however quite fast. As an example: To find a given sequence in a genome it takes the about 10s. So
 * if we have a list of about 30 possible motifs( motifs that surpass the threshold) it would take the algo about 5 min
 * to find all the corresponding positions in the chromosome. While this works fin for fa files with only one sequence
 * it can get quite messy with files that have more then one sequence. One idea is to use this algorithm only if the file
 * contains only one big sequence (for example a chromosome ) and if there are only a few sequences that surpass the
 * given threshold.
 */
/*
 int main() {
     Queuetest a;
    auto start = std::chrono::high_resolution_clock::now();
    std::ifstream is("baba.fna", std::ios::in);
    if (is) {
        is.seekg(0, is.end);
        long unsigned int length = is.tellg();
        std::cout << " length b " << length << std::endl;
        is.seekg(0, is.beg);
        size_t pos(0);
        std::string searched_sequence = "TTTTCACTAACAATTAACAACTCTAAAGGTGGATCTAATCCTTTTTATTCTTACAACCTCTTCAGCAGGTTTTGGATTCC";
        std::string txt;
        char *buffer = new char[length];
        is.read(buffer, length);
        txt = buffer;
        char g;

        for (size_t i(0); i < length; i++) {
            g = buffer[i];
            a.update(g);

        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> duration = end - start;
    std::cout << duration.count() << "s" <<  std::endl;
    return 0;
}*/
        /*
        std::size_t abs_position(0);
        while ( pos < length) {
            abs_position = txt.find(searched_sequence, pos + 1);
            if (abs_position != std::string::npos) {
                std::cout << "found position " << abs_position << std::endl;
                pos = abs_position;
                std::cout << pos << std::endl;
                //pos += 1;
            }else{
                std::cout << "no more positions" << std::endl;
                pos = length;
            }
        }
        delete[] buffer;
        std::cout << "length " << length << std::endl;
    }
    is.close();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> duration = end - start;
    std::cout << duration.count() << "s" <<  std::endl;
        return 0;
         */



/* Code I found on internet for reading fasta file properly
#include <iostream>
#include <fstream>
#include <string>
int main( int argc, char **argv ){
    if( argc <= 1 ){
        std::cerr << "Usage: "<<argv[0]<<" [infile]" << std::endl;
        return -1;
    }
    std::ifstream input(argv[1]);
    if(!input.good()){
        std::cerr << "Error opening '"<<argv[1]<<"'. Bailing out." << std::endl;
        return -1;
    }
    std::string line, name, content;
    while( std::getline( input, line ).good() ){
        if( line.empty() || line[0] == '>' ){ // Identifier marker
            if( !name.empty() ){ // Print out what we read from the last entry
                std::cout << name << " : " << content << std::endl;
                name.clear();
            }
            if( !line.empty() ){
                name = line.substr(1);
            }
            content.clear();
        } else if( !name.empty() ){
            if( line.find(' ') != std::string::npos ){ // Invalid sequence--no spaces allowed
                name.clear();
                content.clear();
            } else {
                content += line;
            }
        }
    }
    if( !name.empty() ){ // Print out what we read from the last entry
        std::cout << name << " : " << content << std::endl;
    }
    return 0;
}
*/


/*
int main(){
    vector<std::string> fasta_files1({"../tests/files/mm9_fasta/chr1.fa", "../tests/files/mm9_fasta/chr2.fa", "../tests/files/mm9_fasta/chr3.fa", "../tests/files/mm9_fasta/chr4.fa"});
    vector<std::string> fasta_files2({"../tests/files/promoters.fasta"});
    std::string bed_file2("../tests/files/bed6cols.bed");
    std::string bed_file3("../tests/files/fasta_bed.bed");
    std::string bed_file4("../tests/files/BMAL1_sites.bed");
    PWM pwm(6,fasta_files2, bed_file1, true );
    std::deque<std::string> seqs1({"GGAAA","CCAATCCACCC","GGTTAGGTGGG"});
    std::deque<std::string> seqs2({"CCAATCCACCC","GGTTAGGTGGG","GGAAA"});
    PWM pwm8(3,fasta_files1,bed_file4, false);
    pwm8.init_pwm();

    return 0;
}*/

/*
int main(){
    Timer t;
    std::ifstream is;
    is.open("../tests/files/mm9_fasta/chr1.fa", std::ios::binary);
    is.seekg(5,is.beg);
    char c = is.get();
  std::cout << c << std::endl;
}
*/

vector<array<double, 4>> matrix;
array<char, 10> buffer;
vector<gen_region> bed_info;
vector<std::string> fasta_files1({"../tests/files/mm9_fasta/chr1.fa", "../tests/files/mm9_fasta/chr2.fa",
                                  "../tests/files/mm9_fasta/chr7.fa",
                                  "../tests/files/mm9_fasta/chr8.fa",
                                  "../tests/files/mm9_fasta/chr11.fa",});

vector<std::string> fasta_files2({"../tests/files/promoters.fasta"});
vector<std::string> fasta_files3({"../tests/files/mm9_fasta/chr7.fa", "../tests/files/mm9_fasta/chr11.fa"});
std::string bed_file2("../tests/files/bed6cols.bed");
std::string bed_file3("../tests/files/fasta_bed.bed");
std::string bed_file4("../tests/files/BMAL1_sites.bed");
std::string bed_file5("../tests/files/mm9_fasta/final_bed.bed");
std::string bed_file6("../tests/files/BMAL1_chr7.bed");
//PWM pwm(6,fasta_files2, bed_file1, true );
//PWM pwmc(5,fasta_files2, bed_file1, true );
std::deque<std::string> seqs1({"GGAAA","CCAATCCACCC","GGGTGGATTGG"});
std::deque<std::string> seqs2({"CCAATCCACCC","GGGTGGATTGG","GGAAA"});
std::deque<std::string> seqs3({"TACCTGGCA","TCCACGA","CTTGG","CTTGG", "CTAGC","AATGG","TCTGG"});
std::vector<std::array<double,4>> mat({ {0.18, 0.54, 0.03, 0.25} ,
                                        {0.22, 0.29, 0.06, 0.43} ,
                                        {0.2, 0.16, 0.08, 0.56} ,
                                        {0.07, 0.12, 0.75, 0.06} ,
                                        {0.1, 0.22, 0.6, 0.08} });


int main(){

    Timer t;
    Optimized_Matrix pwm9(8,fasta_files3,bed_file6, false,0.00001 );
    pwm9.init_pwm();

}



/*
int main(){
    Timer t;
    Reader reader(fasta_files1, mtx_file, score);
    reader.init_reading();
    TReader reader2(fasta_files1);
   reader2.init_reading();
}

*/
