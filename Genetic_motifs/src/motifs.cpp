#include "motifs.h"
#include <math.h> 
#include <sstream>   
	
Motifs::Motifs(double threshold, std::string fa_filename, std::string mtx_filename)
		:score_threshold (threshold), fasta_filename(fa_filename), bufferseq({}), cbufferseq({}), matrix({}), position(0),  chromosome(" "), found_motif(), score(0), score_c(0)
{
    read_matrix(mtx_filename);
    set_length();
    outfile.open("found_motifs_"+ fasta_filename.substr(fasta_filename.rfind('/') +1, fasta_filename.rfind('.') - fasta_filename.rfind('/') -1 ) +".csv");
    if(!outfile.good())
        std::cerr << "error creating file" << std::endl;
    outfile<<"#Motifs found in" << '\t' << fasta_filename << std::endl;
}

void Motifs::read_matrix(std::string filename)
{
	std::ifstream Matrix(filename);
	std::string line;
    try {
        if (!Matrix) {
            std::cerr << "error opening" << filename << std::endl;
            throw std::invalid_argument ("Error file could not open or was not found");
        }

        std::array<double, 4> help;
        while (std::getline(Matrix, line).good()) {
            std::istringstream iss(line);
            std::string score;
            for (size_t i(0); std::getline(iss, score, '\t'); ++i) {
                if(i > 3) throw std::invalid_argument ("Error invalid position weight matrix");
                help[i] = 2 + log2(std::stod(score));
            }
            matrix.push_back(help);
        }
    }catch (std::invalid_argument &s) {
        std::cerr << s.what() << std::endl;
    }

}	

void Motifs::update(char c, unsigned long int pos){
    update_sequences(c);
    position = pos;
    update_score();
    decision_maker();
}	

void Motifs::update_sequences(unsigned int c)
{
    if (bufferseq.size() < length) {
        bufferseq.push_back(c);
        cbufferseq.push_back(transform_to_complementary(c));
    }else{
        bufferseq.pop_front();
        bufferseq.push_back(c);
        cbufferseq.pop_front();
        cbufferseq.push_back(transform_to_complementary(c));
    }
}

void Motifs::update_score()
{
    
	double score_update(0);
	double score_c_update(0);
	int b(0);
	int a(0);
	if(bufferseq.size() == length) {
        for (size_t i(0); i < length; ++i) {
            a = bufferseq[i];
            b = cbufferseq[length-(i+1)];
            if(a != 4 and b != 4) {
                score_update += matrix[i][a];
                score_c_update += matrix[i][b];
            }
        } 
        score = score_update;
        score_c =score_c_update;
    }
}

void Motifs::decision_maker(){
    if(score >= score_threshold and score_c >= score_threshold) {
        found_motif.sequence = transform_bufferseq_to_string(bufferseq);
        found_motif.score = score;
        set_found_motif(' ');
    } else if(score >= score_threshold){
        found_motif.sequence = transform_bufferseq_to_string(bufferseq);
        found_motif.score = score;
        set_found_motif('+');
    } else if(score_c >= score_threshold){
        found_motif.sequence = transform_bufferseq_to_string_comp(cbufferseq);
        found_motif.score = score_c;
        set_found_motif('-');
    }
}



void Motifs:: set_found_motif(char c)
{
	if(c == ' ' or c == '+' or c == '-') {
		found_motif.chromosome = chromosome;
		if(position >= length + 1) {
            found_motif.position = position - length + 1;
        }else{
            found_motif.position = 0;
		}
		found_motif.coding = c;
		print_found_motifs();	
	} else { 
		std::cerr << "Error : character isn't '+', '-' or ' ' " << std::endl;
	}
}

void Motifs:: print_found_motifs()
{
	{
		outfile<<found_motif.chromosome<<'\t'
			<<found_motif.position<<'\t'
			<<found_motif.coding<<'\t'
			<<found_motif.sequence<<'\t'
			<<found_motif.score
			<<std::endl;
	}
}

std::string Motifs::transform_bufferseq_to_string_comp(std::deque <unsigned int> d)
{
    std::string ret;
    char c;
    for(size_t i(0); i < d.size(); ++i){
        switch(d[length -(i+1)]) {
            case (0):
                c = 'A';
                break;

            case (1):
                c = 'C';
                break;

            case (2):
                c = 'G';
                break;

            case (3):
                c = 'T';
                break;


            default:
                c = 'N';
                break;
        }
        ret+=c;
    }
    return ret;
}

std::string Motifs::transform_bufferseq_to_string(std::deque <unsigned int> d)
{
    std::string ret;
    char c;
    for(size_t i(0); i < d.size(); ++i){
        switch(d[i]) {
            case (0):
                c = 'A';
                break;

            case (1):
                c = 'C';
                break;

            case (2):
                c = 'G';
                break;

            case (3):
                c = 'T';
                break;


            default:
                c = 'N';
                break;
        }
        ret+=c;
    }
    return ret;
}

unsigned int Motifs:: transform_to_complementary(unsigned int d)
{
    unsigned int c;
    switch(d) {
        case (0):
            c = 3;
            break;

        case (1):
            c = 2;
            break;

        case (2):
            c = 1;
            break;

        case (3):
            c = 0;
            break;

        default:
            c = 4;
            break;
    }
	return c;
}

void Motifs::reinitialize(std::deque <unsigned int> seq, std::deque <unsigned int> cseq, double score, double cscore) 
{
	set_buffer_seq(seq);
	set_score(score);
	set_cbuffer_seq(cseq);
	set_cscore(cscore);
}

std::deque <unsigned int> Motifs::get_bufferseq()
{
	return bufferseq;
}

std::deque <unsigned int> Motifs::get_cbufferseq()
{
	return cbufferseq;
}

double Motifs::get_score() 
{	
	return score;
}

double Motifs::get_cscore()
{
	return score_c;
}

Found_Motif Motifs::get_found_motif() 
{
	return found_motif;
}

Matrix Motifs::get_matrix()
{
	return matrix;
}

void Motifs::set_length()
{
	length = matrix.size() ;
}

void Motifs::set_score(double set_score) 
{
	score = set_score;
}

void Motifs::set_cscore(double set_score)
{
	score_c = set_score;
}

void Motifs::set_buffer_seq(std::deque <unsigned int> seq)
{
	bufferseq = seq;
}

void Motifs::set_cbuffer_seq(std::deque <unsigned int> seq)
{
	cbufferseq = seq;
}

void Motifs::set_chromosome(std::string seq_name)
{
    chromosome = seq_name;
    bufferseq.clear();
    cbufferseq.clear();
}

void Motifs::reset_bufferseq()
{
    bufferseq.clear();
}

Motifs::~Motifs()
{outfile.close();}
