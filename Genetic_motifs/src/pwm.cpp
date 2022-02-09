#include "pwm.h"


bool is_number(const string& s) {
    string::const_iterator it = s.begin();
    while (it != s.end() && isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

int to_int(char c) {
    int nuc_nb;

        switch (toupper(c)) {
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
    return nuc_nb;
}




seq_t getSequence(ifstream& is, const unsigned long int& lower, const unsigned long int& upper, unsigned long int& current_pos) {

    seq_t ret;
    char c;
    char a;
    while (current_pos < lower and is.good()){
        a = is.get();
        if(a != '\n' and a !='\t' and a!= '\r'){
            ++current_pos;
        }
    }

    if (is.good()) {
        while (is.good() and current_pos < upper) {
            c = is.get();
            if (c == '>') {
                return seq_t({5}); 
            } else if (c != '\t' and c != '\r' and c != '\n') {
                ret.push_back(to_int(c));
                ++current_pos;
            }
        }
        is >> noskipws;
    } else {
        return seq_t({5});
    }

    return ret;
}

seq_t convert_to_complement(const seq_t& seq){
    seq_t ret;
    for(size_t i(seq.size()); i > 0; --i){
        switch(seq[i-1])
        {
            case(0):
                ret.push_back(3);
                break;

            case(1):
                ret.push_back(2);
                break;

            case(2):
                ret.push_back(1);
                break;

            case(3):
                ret.push_back(0);
                break;

            default:
                ret.push_back(4);
                break;
                
        }

    }
    return ret;
}

PWM::PWM(unsigned int motif_length, vector<string> fa_files, string bed_name) :
        fasta_files (fa_files),
        bed_file (bed_name),
        motif_counter(0),
        matrix(motif_length, {0.25, 0.25, 0.25, 0.25})
{
    read_bed();
}

void PWM::read_bed() {
    ifstream bed(bed_file);
    string line;


    if (bed.is_open()) {
        while (bed.good()) {
            getline(bed, line);
            if (!line.empty()) {
                gen_region help(read_bed_line(line));
                gen_regions.insert({{help.seq_name, help.beginning}, help});
            }
        }
    } else {
        cerr << "Error opening " << bed_file << endl;
    }
    bed.close();
    if (gen_regions.empty()) {
        cerr << "No valid sequence in bed file" << endl;
    }
}

gen_region PWM::read_bed_line(string line) {
    string current;
    vector<string> bed_strings;

    for (size_t i(0); i < line.size(); i++) {
        if(line[i] != '\t' and line[i] != '\r' and line[i] != '\n' and line[i] != ' ') {
            current.push_back(line[i]);
        }else if((i+1) < line.size() and line[i+1] != '\t' and line[i+1] != '\r' and line[i+1] != '\n' and line[i+1] != ' '){
            bed_strings.push_back(current);
            current.erase();
        }
    }
    bed_strings.push_back(current);

    gen_region ret;
    bool configured(true);
    if(bed_strings.size() < 3){
        ret.seq_name = "unconfigured";
        ret.beginning = 0;
        ret.end = 0;
        ret.s = BOTH;
        return ret;
    }
    ret.seq_name = bed_strings[0];
    if(is_number(bed_strings[1]) and  is_number(bed_strings[2])) {
        ret.beginning = stoi(bed_strings[1]);
        ret.end = stoi(bed_strings[2]);
    }else{
        configured = false;
    }

    if(bed_strings.back() == "+") {
        ret.s = PLUS;
    }else if(bed_strings.back() == "-") {
        ret.s = MINUS;
    }else{
        ret.s = BOTH;
    }
    if(configured) {
        return ret;
    }else{
        ret.seq_name = "unconfigured"; 
        return ret;
    }
}



void PWM::matrix_counter(int j, size_t i, double weight) {

    if(0 <= j and j < 4) {
        matrix[i][j] += weight;
    } else {
        matrix[i][0] += (weight * 0.25);
        matrix[i][1] += (weight * 0.25);
        matrix[i][2] += (weight * 0.25);
        matrix[i][3] += (weight * 0.25);
    }
}

void PWM::compute_matrix(seq_t sequence) {
    unsigned int seq_length = sequence.size();
    unsigned int motif_lenght = matrix.size();
    unsigned int difference = seq_length -motif_lenght;
    double weight(1.0/(difference + 1)); 
    if(seq_length > motif_lenght) {
        
        for(size_t i(0); i <= difference; ++i) {
            for(size_t j(i); j < motif_lenght + i; ++j) {
                matrix_counter(sequence[j], j-i, weight);
            }
        }

    } else if(seq_length == motif_lenght) {
        for(size_t i(0); i < motif_lenght; i++) {
            matrix_counter(sequence[i], i, weight);
        }
    }
}

void PWM::divide_matrix(long double nb_of_motifs) {
    for (size_t i=0; i<matrix.size(); ++i) {
        for (size_t j=0; j<=3; ++j) {
            matrix[i][j] = matrix[i][j]/(nb_of_motifs);
        }
    }
}

void PWM::find_sequences_on_chromosome(string chr, deque<seq_t>* sequences, const multimap<pair<string, unsigned int>, gen_region>* gen_regions,
                                       unsigned int length){
    string title;
    string line;
    seq_t help;
    unsigned long int current_pos(0);

    char Buffer[100000];
    ifstream is;
    is.rdbuf()->pubsetbuf(Buffer, 100000);
    is.open ( chr );
    if(!is) cerr<< "Error opening file" << chr << endl;
    string ret;
    while(is.good()) {
        if (getline(is, line).good() and line[0] == '>') {
            current_pos = 0;
            title = line.erase(0, 1); 
            for (auto it = gen_regions->lower_bound({title,0}); it->first.first == title; it++) {
                help = getSequence(is, it->second.beginning, it->second.end, current_pos);
                if (help != seq_t{5} and help.size() >= length) {
                    lock_guard<mutex> lock(m1);
                    if (it->second.s == PLUS) {
                        sequences->push_back(help);
                    } else if (it->second.s == MINUS) {
                        sequences->push_back(convert_to_complement(help));
                    } else {
                        sequences->push_back(help);
                    }
                }
            }
        }
    }
   
    is.close();
}



void PWM::addSequence(const seq_t& sequence){
    sequences.push_back(sequence);
}



void PWM::start_searching() {

#define ASYNC 1
#if ASYNC

        for (auto &file : fasta_files) {
            fut_vec1.push_back(async(launch::async, find_sequences_on_chromosome, file, &sequences, &gen_regions,getMatrix().size()));
        }
        for (size_t i(0); i < fut_vec1.size(); ++i) {
            fut_vec1[i].wait();
        }

#else

        for (auto &file : fasta_files) {
            find_sequences_on_chromosome(file, &sequences, &gen_regions, getMatrix().size());
        }

#endif

}


matrix_t PWM::getMatrix(){
    return matrix;
}

void PWM::matrix_output(){
    ofstream outfile;
    outfile.open ("output.pwm");
    outfile<<"#INCLUSive Motif Model"<<endl;
    for (size_t i=0; i<matrix.size(); ++i) {
        for (size_t j(0); j<4; ++j) {
            outfile<<matrix[i][j]<<'\t';
        }
        outfile<< '\n';
    }
    outfile.close();
}

multimap<pair<string, unsigned int>, gen_region> PWM::getGenRegions(){
    return gen_regions;
}

void PWM::init_pwm(){
    start_searching();
    for(auto& sequence : sequences){
        compute_matrix(sequence);
    }
}

void PWM::start_pwm(){
    init_pwm();
    divide_matrix(sequences.size()+1);
    matrix_output();

}

deque<seq_t> PWM::getSeqeunces(){return sequences;}
