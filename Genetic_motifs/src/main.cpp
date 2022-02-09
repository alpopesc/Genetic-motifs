#include <tclap/CmdLine.h>
#include <iostream>
#include "simulation.h"
#include "Timer.h"
#include "reader.h"
using namespace TCLAP;
std::vector<std::string> a({"../tests/files/promoters.fasta"});
int main(int argc, char **argv) {
    Timer t;

	try {
		Simulation Sim(argc, argv);
		Sim.run();

	}
	catch(std::string &e){
		std::cerr << e << std::endl;
		return 1;
	}
	catch(TCLAP::ArgException &e)
	{
	    std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
	}

    return 0;

}



