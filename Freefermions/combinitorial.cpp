#include <iostream>
#include <iomanip>
#include "LanczosPlusPlus/src/Models/HubbardOneOrbital/BasisOneSpin.h"
#include "../codes/PsimagLite/src/Matrix.h"

PsimagLite::Matrix<SizeType> LanczosPlusPlus::BasisOneSpin::comb_;
SizeType LanczosPlusPlus::BasisOneSpin::nsite_;

using BasisOneSite = LanczosPlusPlus::BasisOneSpin;

/* extern "C" {
    void dgeev_(char* jobvl, char* jobvr, int* n, double* a, int* lda,
                double* wr, double* wi, double* vl, int* ldvl,
                double* vr, int* ldvr, double* work, int* lwork, int* info);
}*/

// <gs|sigmax_0 sigmax_site |gs>
double computeCorrelation(unsigned int site, const std::vector<double>& gs, const BasisOneSite& bos)
{
	unsigned int hilbert = gs.size();
	double sum = 0;
	for (unsigned int alpha = 0; alpha < hilbert; ++alpha) {
		// This masks two sites: site site, and site 0
		unsigned mask = (1<<site) | 1;
		if (site == 0) mask = 0;
		unsigned int stateAlpha = bos[alpha];
		unsigned stateBeta = stateAlpha ^ mask;
		unsigned int beta = bos.perfectIndex(stateBeta); // the index for state stateBeta
		assert(beta < gs.size()); // and also for alpha 
		sum += gs[alpha]*gs[beta];
	}

	return sum;
}

int main(int argc, char* argv[])
{
	if (argc != 4) {
		std::cerr<<"USAGE: "<<argv[0]<<" number_of_sites J\n";
		return 1;
	}

	bool periodic = true;
	// Number of sites
	unsigned int n = atoi(argv[1]); // don't use atoi atof std:: something to double need to search on the web
	unsigned int npart = atoi(argv[2]);
	double couplingJ = atof(argv[3]);

	// Basis of one site |0> and |1>
	// n sites: Hilbert space is tensor product of n sites: 2^n states for the basis
	// Example: n = 3, example basis state |alice> = |0> | 1> |1>
	// Norm squared of the basis states is 1: <0|0> = 1, <1|1> = 1
	//
	// H = J*\sum_{i, i + 1} \sigmax_i \otimes \sigmax_{i + 1} 
	// Example n = 3        ---> H = sigmax_0 sigmax_1 + sigmax_1 sigmax_2 + sigmax_2 sigmax_0
	// There are a identities implied!

	// sigmax acts only on one site: sigmax |0> = |1> and sigmax|1> = |0>
	// 
	// Example H |alice> = |1> | 0> |1> + | 0> |0> |0> + |1> | 1> |0>
	//
	// Fill the matrix H for an Ising model
	//
	LanczosPlusPlus::BasisOneSpin bos(n, npart);
	unsigned int hilbert = bos.size();
	std::cout<< "hilbert " << hilbert<<"\n";	
 	PsimagLite::Matrix<double> hamiltonian(hilbert, hilbert);	

	// More intuitive way maybe:
	// loop over row and nested loop over col <col|H|row> = complex number = H(row, col)
	//
	
	for (unsigned int row = 0; row < hilbert; ++row) {
		// Every integer from 0 to 2^n - 1 represents a state of the (computational) basis
		// Example Alice looks like 0 1 1 --> interpret it as a binary number --> 3
		// H|row> = sum of states of the basis with coefficients
		// For example H|3> = H|alice> = 1*|5> + 1*|0> + 1*|6>
		// H|3> = [1, 0, 0, 0, 0, 1, 1, 0]
		//
		unsigned int statei = bos[row];
		std::cout<<"statei"<<statei<<"\n";

		for (unsigned int site = 0; site < n; ++site) {
			
			// Starting with state row 
			// flip bit number site and bit number site + 1
			// We need to worry about border, like if site + 1 >= n
			// row =      bit(n-1) bit(n-2) ... bit(site+1) bit(site) ... bit(1) bit(0)
			// mask    =    0        0      ...  1            1       ... 0       0
			// xor op. -----------------------------------------------------------------------
			// result =    bit(n-1) bit(n-2)   bar bit(site+1) bar bit(site) ... bit(1) bit(0)
			// The operation left shift "1 << n " shifts the number 1 n places to the left
			/* unsigned int next_site = site + 1;
			if (site + 1 == n) {
				next_site = 0;
			}*/
			// "ternary operator in C and C++"
			
			bool isborder = ((site + 1) == n);
			
			if (isborder and !periodic) {
				continue; // or break
			}

			unsigned int next_site = (isborder) ? 0 : site + 1;
			unsigned int masksite = (1<<site);
			unsigned int masknext = (1<<next_site);
			unsigned int resultsite = statei & masksite;
			unsigned int resultnext = statei & masknext;
			if (resultsite == 0 && resultnext == 0) continue;
			if (resultsite > 0 && resultnext > 0) continue;

			unsigned int mask = (1<<site) | (1<<next_site); // The "|" or binary op
			unsigned int statej = statei ^ mask; // ^ is the xor binary op
			unsigned int col = bos.perfectIndex(statej);
			//std::cout<<"if I apply cdagger_"<<site<<" c_"<<next_site<<" to state ";
			//std::cout<<statei<<" ("<<row<<") I get state "<<statej<<" which has index "<<col<<"\n";
			hamiltonian(row, col) = couplingJ;
		}
	}

	std::cout<<hamiltonian;
	std::cout.precision(12);
	// Think about the ground state, what energy does it have?
	// Change the sign of the coupling and see if something qualitatively diffrent happens
	std::vector<double> eigs(hilbert);

	diag(hamiltonian, eigs, 'V');
	std::cout<<"gound state energy "<<eigs[0]<<"\n";
	std::vector<double> gs(hilbert);
	for (unsigned int k = 0; k < hilbert; ++k) {
		gs[k] = hamiltonian(k, 0);
	}

	for (unsigned int site = 0; site < n; ++site) {
		double correlation = computeCorrelation(site, gs, bos);
		std::cout<<site<<" "<<correlation<<"\n";
	}
	
	std::cout<<"\n";



}

