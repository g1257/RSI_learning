#include <iostream>
#include <iomanip>
#include <cassert>
#include "LanczosPlusPlus/src/Models/HubbardOneOrbital/BasisOneSpin.h"
#include "PsimagLite/src/Matrix.h"
#include "PsimagLite/src/BitManip.h"

PsimagLite::Matrix<SizeType> LanczosPlusPlus::BasisOneSpin::comb_;
SizeType LanczosPlusPlus::BasisOneSpin::nsite_;

using BasisOneSite = LanczosPlusPlus::BasisOneSpin;


// Symmetrize a matrix (unused for now)
void symmetrize(PsimagLite::Matrix<double>& matrix)
{
	unsigned int rank = matrix.rows();
	assert(rank == matrix.cols()); // matrix must be square
	for (unsigned int row = 0; row < rank; ++row) {
		for (unsigned int col = row + 1; col < rank; ++col) {
			matrix(col, row) = matrix(row, col);
		}
	}
}

// Sign of c^\dagger_ind c_jnd applied to state |state>
int fermionSign(unsigned int state, unsigned int ind, unsigned int jnd, int fermion_sign)
{
	if (ind == jnd) {
		err(std::to_string(ind) + " must be different than " + std::to_string(jnd) + "\n");
	}


	unsigned int maskindless = (1 << ind) - 1;
	unsigned int maskjndless = (1 << jnd) - 1;

	unsigned int tmpindless = (maskindless & state);
	unsigned int tmpjndless = (maskjndless & state);
	unsigned int count = PsimagLite::BitManip::countKernighan(tmpindless) + PsimagLite::BitManip::countKernighan(tmpjndless);

	if (ind > jnd) {
		++count;
	}

	int val = (count & 1) ? fermion_sign : 1;

	return val;
}

double computeCorrelation(unsigned int ind, unsigned int jnd, const std::vector<double>& gs, const BasisOneSite& bos, int fermion_sign)
{
	unsigned int hilbert = gs.size();
	double sum = 0;
	bool sitesEqual = (ind == jnd);
	for (unsigned int alpha = 0; alpha < hilbert; ++alpha) {
		unsigned int maski = (1<<ind);
		unsigned int maskj = (1<<jnd);
		// This masks two sites: site ind, site jnd
		unsigned mask = (sitesEqual) ? 0 : maski | maskj;
		unsigned int stateAlpha = bos[alpha];
		if (!sitesEqual && (stateAlpha & maski)) continue; // can't apply cdagger at ind
		unsigned stateBeta = stateAlpha ^ mask;
		if (!sitesEqual && (stateBeta & maskj)) continue; // can't apply c at jnd
		unsigned int beta = bos.perfectIndex(stateBeta); // the index for state stateBeta
		assert(beta < gs.size()); // and also for alpha
		int fs = (sitesEqual) ? 1 : fermionSign(stateAlpha, ind, jnd, fermion_sign);
		sum += gs[alpha]*gs[beta]*fs;
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
	int FERMION_SIGN = -1; // set to -1 for fermions, 1 for bosons
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
			unsigned int resultsite = (statei & masksite) ? 1 : 0;
			unsigned int resultnext = (statei & masknext) ? 1 : 0;

			if (resultsite == resultnext) continue;

			unsigned int mask = (1<<site) | (1<<next_site); // The "|" or binary op
			unsigned int statej = statei ^ mask; // ^ is the xor binary op
			unsigned int col = bos.perfectIndex(statej);

			int fs = 1;
			if (resultsite == 1 && resultnext == 0) {
				// destroy at site and create at next
				fs = fermionSign(statei, next_site, site, FERMION_SIGN);
			} else {
				fs = fermionSign(statei, site, next_site, FERMION_SIGN);
			}


			hamiltonian(row, col) = couplingJ * fs;
		}
	}

	std::cout<<hamiltonian;
	std::cout.precision(12);
	// Think about the ground state, what energy does it have?
	// Change the sign of the coupling and see if something qualitatively diffrent happens
	std::vector<double> eigs(hilbert);

	diag(hamiltonian, eigs, 'V');
	std::cout<<"ground state energy "<<eigs[0]<<"\n";
	std::vector<double> gs(hilbert);
	for (unsigned int k = 0; k < hilbert; ++k) {
		//std::cout<<eigs[k]<<" <---\n";
		gs[k] = hamiltonian(k, 0);
	}

	std::cout<<"cdaggeri cj\n";
	for (unsigned int i = 0; i < n; ++i) {
		for (unsigned int j = 0; j < n; ++j) {
			double correlation = computeCorrelation(i, j, gs, bos, FERMION_SIGN);
			std::cout<<correlation<<" ";
		}

		std::cout<<"\n";
	}

	std::cout<<"\n";
}
