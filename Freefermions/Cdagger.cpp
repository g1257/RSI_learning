#include "../../codes/PsimagLite/src/Matrix.h"
#include <iomanip>
#include <iostream>

// states | 0 0 1 0> for N = 4
// c|0> = 0 and c|1> = 0; c
// c^\dagger |0> = |1>    c^\dagger|1> = 0
// c^\dagger_2| 1 0 1 0> =  -| 1 1 1 0> (fermion sign)
//
// H = \sum_ij A(i, j) c^\dagger_i c_j + H.C., but A_ii = 0 for all
// ind must be diffrent than jnd
// HW correct this for when they are equal
//
//

double density(unsigned int j, unsigned int hilbert, PsimagLite::Matrix<double>& hamiltonian)
{

	unsigned int maskj = (1 << j);

	double expect = 0;
	for (unsigned int k = 0; k < hilbert; ++k) {

		expect += (k & maskj) ? hamiltonian(k, 0) * hamiltonian(k, 0) : 0;
	}

	return expect;
}

double correlation(unsigned int ind, unsigned int jnd, const std::vector<double>& gs)
{
	double sum = 0;
	for (unsigned int a = 0; a < gs.size(); ++a) {
		unsigned int maskind = (1 << ind);
		unsigned int maskjnd = (1 << jnd);

		if (ind == jnd) {
			sum += gs[a] * gs[a];
			continue;
		} else if ((maskjnd & a) == 0 || (maskind & a)) {
			continue;
		}
		unsigned int mask = maskind | maskjnd;
		unsigned int b = a ^ mask;

		sum += gs[a] * gs[b];
	}

	return sum;
}

PsimagLite::Matrix<double> correlations(unsigned int n, const std::vector<double>& gs)
{

	PsimagLite::Matrix<double> corr(n, n);
	for (unsigned int i = 0; i < n; ++i) {
		for (unsigned int j = 0; j < n; ++j) {
			corr(i, j) = correlation(i, j, gs);
		}
	}
	return corr;
}

double square(double x)
{

	return x * x;
}

int main(int argc, char* argv[])
{

	if (argc != 2) {

		std::cerr << "Usage: " << argv[0] << "\n";
		return 1;
	}

	unsigned int N = atoi(argv[1]);
	unsigned int hilbert = (1 << N);
	PsimagLite::Matrix<double> hamiltonian(hilbert, hilbert);
	unsigned int middle1 = N / 2;
	bool periodic = true;
	// This is the electric field at a site
	double v_0 = 0;
	// A_ii = v_0 if i = middle1 or midddle2
	// zero o.w.
	for (unsigned int col = 0; col < hilbert; ++col) {
		unsigned int maskmid1 = (1 << middle1);
		hamiltonian(col, col) = (col & maskmid1) ? v_0 : 0;

		// H|row> = \sum_{col} factor* |col>
		// H(row, col) = factor
		// The first and last column will always be all 0
		for (unsigned int i = 0; i < N; i++) {
			unsigned int j = i + 1;
			if (periodic) {
				if (j == N) {

					j = 0;
				}
			} else {
				if (j == N)
					continue;
			}

			unsigned int maski = (1 << i);
			unsigned int maskj = (1 << j); // puts a 1 in location j,
			unsigned int mask = maski | maskj;
			unsigned int ri = (col & maski) ? 1 : 0; // if not zero then true
			unsigned int rj = (col & maskj) ? 1 : 0;
			// We only want one 'hit'
			// putting this snip here for effecacy
			if (ri == rj) {
				continue;
			}

			unsigned int row = col ^ mask;

			// fermion sign
			// n.particles = sites with the bit 1 ("the bit is set")
			// nci = n. particle before site i = \sum_{l < i} n_l
			// ncj = ... \sum_{l < j} n_l
			// factor = (-1)^{nci+ncj} = (-1) \sum_{l < i} {n_l} + \sum_{l < i} {n_l} + \sum_{i <= l < j} {n_l}}=
			// factor = (-1)^{2*\sum_{l < i} {n_l} + \sum_{i <= l < j} {n_l}} } = (-1)^\sum_{i <= l < j} {n_l}}
			// If i < j which it is in this loop
			// count between i and j
			// Maybe j does not matter here
			unsigned int maskbari = maski - 1; // = 2^i - 1 has 1 in all bits 00001111 for l < 1 is 1
			unsigned int maskbarj = maskj - 1;
			unsigned int counti = (maskbari & col);
			unsigned int countj = (maskbarj & col);
			int signi = (counti > 0) ? -1 : 1;
			int signj = (countj > 0) ? -1 : 1;
			int signextra = (ri > 0) ? -1 : 1;
			// unsigned int maskbar = maskbari & maskbarj;

			// ncount is the number of bits in count ==> factor = (ncount & 1) ? -1 : 1
			// int ncount = std::popcount(count);
			hamiltonian(row, col) = 1; // signi * signj *signextra;
						   // factor = 1 or -1
						   // We include the ones that should be zero
						   // \daggeri      .j
						   // 0              1
						   // 1              0  <--- this situation is the Hermitian conjugate
						   //
						   // 0             0  <-- gives 0
						   // 1             1  <-- gives 0
		}
	}

	// Print results
	std::cout << hamiltonian;
	std::vector<double> eigs(hilbert);

	diag(hamiltonian, eigs, 'V');
	double energyin = eigs[0]; /// N;

	std::cout.precision(12);
	std::cout << energyin << "\n";

	/*std::vector<double> gs(hilbert);
	for (unsigned int k = 0; k < hilbert; ++k  )
	{

		gs[k] = hamiltonian(k,0);

	}
	*/
	// <gs|c^dag_1 c_1| gs>

	// std::cout<<gs<<"\n";
	// PsimagLite::Matrix<double> corr = correlations(N,gs);
	// std::cout<<corr;

	/*
	std::vector<double> dens(N);
	for (unsigned int j = 0; j < N; ++j)
	{

	dens[j] = density(j,hilbert, hamiltonian);


	}
	*/
}
