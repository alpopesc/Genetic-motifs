#include <cmath>
#include "extension.h"



Optimized_Matrix::Optimized_Matrix(unsigned int motif_length, vector <string> fa_files, string bed_name, double threshold)
	:PWM(motif_length, fa_files, bed_name), new_matrix(motif_length, {0.25, 0.25, 0.25, 0.25}),
	threshold_(threshold)
{
	init_pwm();
	old_matrix = getMatrix();
}

void Optimized_Matrix::set_old_matrix(matrix_t m) {
	old_matrix = m;
}

void Optimized_Matrix::set_new_matrix(matrix_t m) {
	new_matrix = m;
}

void Optimized_Matrix::iterate()
{
	for(size_t i(0); i< 100 and !equal_matrix(); ++i)

	{
	    matrix = vector<array<double,4>>(matrix.size(), {0.25, 0.25, 0.25, 0.25});
	    best_score_sequences.clear();
		for (auto& sequence : sequences )
		{
			best_score_sequences.push_back(find_best_score_sequences(sequence));
		}

        for (auto& sequence : best_score_sequences )
        {
            compute_matrix(sequence);
        }
        old_matrix = new_matrix;
        new_matrix = getMatrix();
		cout << i << endl;

	}
	
    divide_matrix(sequences.size() +1);
	matrix_output();
}


bool Optimized_Matrix::equal_matrix() {

    if(old_matrix.size() != new_matrix.size())
        return false;

    
    for(size_t row = 0; row < old_matrix.size(); ++row) {
        for(size_t col = 0; col < 4; ++col) {
            if(abs(old_matrix[row][col] - new_matrix[row][col]) > threshold_)
                return false;
        }
    }
    return true;
}




seq_t Optimized_Matrix::find_best_score_sequences(seq_t sequence)
{
	double score(0.);
	double best_score(0.);
	seq_t best_seq;
	seq_t help;
	for(size_t j(0); j < sequence.size() - matrix.size(); ++j){
		for (size_t i(j); i - j < old_matrix.size() ; ++i){
            score += old_matrix[i-j][sequence[i]];
            help.push_back(sequence[i]);
		}
        if (score > best_score)
        {
            best_score = score;
            best_seq = help;
        }
        help.clear();
		score = 0.0;
	}
	return best_seq;
}
