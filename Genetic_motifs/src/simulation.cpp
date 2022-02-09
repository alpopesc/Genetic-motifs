#include "simulation.h"

using namespace TCLAP; 

Simulation::Simulation(int argc, char** argv) {

genome_="nul";
matrix="nul";
score_=0;
usage_initial=false;
length_=0;
optimised=false;
init(argc, argv);
	
}

void Simulation::init(int argc, char** argv) {

    try {


        CmdLine cmd("Motifs Team_8");
        /*!
        Allows the user to choose with which input he wants to use the program
        */
        ValueArg<std::string> bedFile("B", "bed_file", "fichier BED à lire", true, "/dev/null", "filepath");

        ValueArg<std::string> matrixFile("m", "matrixFile", "Matrice à lire", true, "/dev/null", "filepath");

        MultiArg<std::string> fastaFiles("f", "position", "liste de positions à lire", true, "filepaths");
        cmd.add(fastaFiles);

        ValueArg<double> score("s", "score", "score minimum pour la matrice", true, 0, "double");

        ValueArg<bool> optimisedMatrix("o", "optimised_matrix", "do you want an optimised matrix", false, 0, "bool");
        cmd.add(optimisedMatrix);


        ValueArg<double> optimisation_threshold("d", "difference_threshold", "score minimum pour la matrice", false,
                                                0.01,
                                                "double");
        cmd.add(optimisation_threshold);


        ValueArg<unsigned int> length("l", "length", "length of the motif", true, 0, "unsigned int");

        cmd.xorAdd(score, length);
        cmd.xorAdd(matrixFile, bedFile);
        cmd.parse(argc, argv);


        if (matrixFile.isSet()) {
            usage_initial = true;
            if (length.isSet()) throw std::invalid_argument("The length should not be provided for this usage of the program");
            if (score.getValue() <= 0) throw std::invalid_argument("The score can not be negative or equal to 0");
            if(optimisation_threshold.isSet() or optimisedMatrix.isSet()) throw std::invalid_argument("No matrix optimization possible for this usage of the program");
        } else if (bedFile.isSet()) {
            usage_initial = false;
            if (score.isSet()) {
                throw std::invalid_argument("The score should not be provided for this usage of the program");
            }

        }

        genome_ = bedFile.getValue();
        matrix = matrixFile.getValue();
        position = fastaFiles.getValue();
        score_ = score.getValue();
        length_ = length.getValue();
        optimised = optimisedMatrix.getValue();
        optimisation_t = optimisation_threshold.getValue();
    }
    catch (std::invalid_argument &s) {
        std::cerr << s.what() << std::endl;
    }


}

void Simulation::run()
{
	if(!optimised){
	if(usage_initial) {
	try {
        if (position.empty()) throw std::invalid_argument ("No FASTA file entry"); 
        Reader r(position, matrix, score_);
        r.init_reading();
        }

    catch(std::invalid_argument &s) {
        	throw (std::invalid_argument("Error due to " + position[0] + " file"));
	  }
	}
    
    else {

			PWM p(length_, position, genome_);
			p.start_pwm();

	}
	}
	else if(optimised)
	{
		Optimized_Matrix O(length_, position, genome_, optimisation_t);
		O.iterate();

	}
}


	









