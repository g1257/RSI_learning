#include <cmath>
#include <iostream>
#include <string>

// Exact polynomial one-particle "free" solution for Fermions
//
// Complte solutions.txt table for fermions with this one and
// with the LanczosPlusPlus code putting FermionSign=-1
// LanczosPlusPlus is doing the combinatorial problem
double doPeriodic(unsigned int nsites, unsigned int npart)
{
	// unused npart for now because looking for the g.s. of all the sectors

	double sum = 0;
	for (unsigned m = 0; m < nsites; ++m) {
		double k = M_PI * m * 2. / nsites;
		double one_particle_ek = -2 * cos(k);
		if (one_particle_ek < 0) {
			sum += one_particle_ek;
		}
	}

	return sum;
}

double doOpen(unsigned int nsites, unsigned int npart)
{
	// unused npart for now because looking for the g.s. of all the sectors

	double sum = 0;
	for (unsigned m = 0; m < nsites; ++m) {
		double k = M_PI * (m + 1) / (nsites + 1);
		double one_particle_ek = -2 * cos(k);
		if (one_particle_ek < 0) {
			sum += one_particle_ek;
		}
	}

	return sum;
}

int main(int argc, char* argv[])
{
	if (argc != 4) {
		std::cout << "USAGE " << argv[0] << " nsites npart periodic | open\n";
		return 1;
	}

	std::cout.precision(12);
	unsigned int nsites = std::atoi(argv[1]);
	unsigned int npart = std::atoi(argv[2]);
	std::string border(argv[3]);
	double energy = 0;
	if (border == "periodic") {
		energy = doPeriodic(nsites, npart);
	} else if (border == "open") {
		energy = doOpen(nsites, npart);
	} else {
		std::cout << "Expecting periodic or open, not " << border << "\n";
		return 2;
	}

	std::cout << "Energy for " << border << " for " << nsites << " sites is ";
	std::cout << energy << " intensive= " << (energy / nsites) << "\n";
}
