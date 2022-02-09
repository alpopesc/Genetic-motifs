#ifndef TEAM8_SIMULATION_H
#define TEAM8_SIMULATION_H

#include <tclap/CmdLine.h>
#include <iostream>
#include "motifs.h"
#include "pwm.h"
#include "reader.h"
#include "extension.h"


/*! \class Simulation
  This is the main class. 
  It manages user inputs, defines the simulation parameters.
  It then runs the simulation and prints the results to the output streams.
  These streams can be files with names based on the string \ref output with a suffix.

  Simulation parameters:
  -a bed file yith position to consider
  * or 
  -a position weight mass matrix
  -a file with the genome
  -a threshold score

 */

class Simulation {
    
public:
        
/*!
  Constructor based on user inputs, takes command-line arguments and passes them to \ref parse.
 */
   Simulation(int argc, char **argv);
/*!
   Uses [TCLAP](http://tclap.sourceforge.net/html/index.html) to parse user inputs.
*/
   void init(int argc, char **argv);
/*!
  Runs the simulation and give the output expected.
 */
    void run();

private:

    void parse(int, char**);
	std::string genome_;
	std::string matrix;
	double score_;
	double optimisation_t;
	std::vector<std::string> position;
	bool usage_initial;
	bool optimised;
    unsigned int length_;

};


#endif
