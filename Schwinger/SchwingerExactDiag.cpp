#include <string>
#include <vector>
#include <fstream>

#include "Matrix.h"
#include "PsimagLite.h"

using MatrixType = PsimagLite::Matrix<double>;
using VectorType = std::vector<double>;
using VectorUintType = std::vector<unsigned int>;

class Params {
public:
	Params(const std::string& filename)
	{
		std::ifstream fin(filename);
		if (!fin || fin.bad()) {
			err("Params::ctor(): Cannot open file " + filename + "\n");
		}

		while (!fin.eof()) {
			std::string stype;
			std::string label;
			fin >> stype;
			if (stype == "") {
				continue;
			}

			fin >> label;
			if (stype == "real") {
				double val = 0;
				fin >> val;
				set(label, val);
			} else if (stype == "integer") {
				unsigned int ival = 0;
				fin >> ival;
				set(label, ival);
			} else {
				fin.close();
				err("Unknown type " + stype + " must be real or integer\n");
			}
		}

		std::cerr<<N_<<" "<<n_<<" "<<t_<<" "<<mass_<<"\n";
	}

	void set(const std::string& label, unsigned int value)
	{
		if (label == "N") {
			N_ = value;
		} else if (label == "n") {
			n_ = value;
		} else {
			err("Params::set(): unknown param " + label + "\n");
		}
	}

	void set(const std::string& label, double value)
	{
		if (label == "mass") {
			mass_ = value;
		} else if (label == "t") {
			t_ = value;
		} else {
			err("Params::set(): unknown param " + label + "\n");
		}
	}

	double get(const std::string& label) const
	{
		if (label == "N") {
			return N_;
		} else if (label == "n") {
			return n_;
		} else if (label == "mass" || label == "m") {
			return mass_;
		} else if (label == "t") {
			return t_;
		} else {
			err("Params::get(): unknown param " + label + "\n");
		}

		return 0;
	}

private:

	double mass_ = 0;
	double t_ = 0;
	unsigned int N_ = 0;
	unsigned int n_ = 0;
};

using ParamsType = Params;

double integerToField(int ind, unsigned int n)
{
	int max = (n & 1) ? (n - 1)/2 : n/2;
	return ind - max;
}

unsigned int computeLanySite(unsigned int prev, unsigned int charge, unsigned int n)
{
	unsigned int l = prev + 1;
	if (charge && l >= n) {
		l = 0;
	}

	if (!charge) {
		l = (prev > 0) ? prev - 1 : n - 1;
	}

	return l;
}

unsigned int computeL(unsigned int prev, unsigned int charge, unsigned int site, unsigned int n)
{
	unsigned int charge_modif = (site & 1) ? charge - 1 : charge;
	return computeLanySite(prev, charge_modif, n);
}

double measureSigma(const VectorType& gs, const ParamsType& params)
{
	unsigned int N = params.get("N");
	unsigned int n = params.get("n");
	unsigned int two_to_the_N = (1 << N);
	double sum1 = 0;
	for (unsigned int l = 0; l < n; ++l) {
		for (unsigned int i = 0; i < two_to_the_N; ++i) {
			double gs_val = gs[i + l*two_to_the_N];
			unsigned int prev = l;
			double sum2 = 0;
			for (unsigned int site = 0; site < N; ++site) {
				unsigned int mask1 = (1 << site);
				unsigned int e1 = (i & mask1) ? 1 : 0;
				unsigned int l1 = computeL(prev, e1, site, n);
				sum2 += l1;
				prev = l1;
			}

			sum1 += sum2*gs_val*gs_val; // times complex conjugate in reality
		}
	}

	return sum1/N;
}

void addNonDiagonalOneL(MatrixType& m, const VectorType& d, unsigned int ind, const ParamsType& params)
{
	unsigned int N = params.get("N");
	unsigned int n = params.get("n");
	unsigned int t = params.get("t");
	unsigned int two_to_the_N = (1 << N);
	unsigned int hilbert = two_to_the_N*n;
	double factor = -t*n*0.5/M_PI;
	m.resize(hilbert, hilbert);
	for (unsigned i = 0; i < two_to_the_N; ++i) {
		unsigned int row = ind*two_to_the_N + i;
		m(row, row) = d[row];

		unsigned int prev = ind;
		// Open BC
		for (unsigned int site = 0; site < N - 1; ++site) {
			unsigned int mask1 = (1 << site);
			unsigned int mask2 = (1 << (site + 1));
			unsigned int e1 = (i & mask1) ? 1 : 0;
			unsigned int e2 = (i & mask2) ? 1 : 0;
			unsigned int l1 = computeL(prev, e1, site, n);
			if (e1 == e2) {
				prev = l1;
				continue;
			}

			unsigned int j = (i ^ (mask1 | mask2));
			bool add_l = (e1 < e2);

			unsigned int l2 = l1 + 1;
			if (add_l && l2 >= n) {
					l2 = 0;
			}

			if (!add_l) {
				l2 = (l1 > 0) ? l1 - 1 : n - 1;
			}

			unsigned int col = j + two_to_the_N*l2;
			assert(col < m.cols());
			// No fermion sign because it's one dimensional with OBC
			m(row, col) = factor;
			prev = l1;
		}
	}
}

void addNonDiagonal(MatrixType& m, const VectorType& diagonal, const ParamsType& params)
{
	unsigned int n = params.get("n");
	for (unsigned i = 0; i < n; ++i) {
		addNonDiagonalOneL(m, diagonal, i, params);
	}
}

VectorType buildDiagonalSectorZero(const ParamsType& params)
{
	unsigned int N = params.get("N");
	unsigned int n = params.get("n");
	unsigned int m = params.get("m");
	unsigned int two_to_the_N = (1 << N);
	double factor = m*n*0.5/M_PI;
	VectorType dvector(two_to_the_N);

	for (unsigned i = 0; i < two_to_the_N; ++i) {
		double sum = 0;
		for (unsigned int site = 0; site < N; ++site) {
			double sign = (site & 1) ? -1 : 1;
			unsigned int mask = (1 << site);
			if (mask & i) {
				sum += sign;
			}
		}

		dvector[i] = sum*factor;
	}

	return dvector;
}

void setDiagonal(VectorType& dvector, const VectorType& dvector0, unsigned int ind, unsigned int n)
{
	double etilde = integerToField(ind, n);
	const double e_squared = etilde * etilde;
	unsigned int two_to_the_N = dvector0.size();
	unsigned int offset = ind*two_to_the_N;
	for (unsigned i = 0; i < two_to_the_N; ++i) {
		dvector[i + offset] = dvector0[i] + e_squared;
	}
}

VectorType computeDiagonal(const ParamsType& params)
{
	unsigned int N = params.get("N");
	unsigned int n = params.get("n");
	unsigned int two_to_the_N = (1 << N);
	unsigned int hilbert = n*two_to_the_N;

	VectorType dvector0 = buildDiagonalSectorZero(params);

	VectorType dvector(hilbert);
	for (unsigned int i = 1; i < n; ++i) {
		setDiagonal(dvector, dvector0, i, n);
	}

	return dvector;
}

void buildHamiltonian(MatrixType& m, const ParamsType& params)
{
	VectorType diagonal = computeDiagonal(params);
	addNonDiagonal(m, diagonal, params);
}

void buildGroundState(VectorType& gs, const ParamsType& params)
{
	MatrixType hamiltonian;
	buildHamiltonian(hamiltonian, params);

	unsigned int hilbert = hamiltonian.rows();
	assert(hilbert == hamiltonian.cols());

	std::vector<double> eigs(hilbert);
	diag(hamiltonian, eigs, 'V');

	gs.resize(hilbert);
	for (unsigned int i = 0; i < hilbert; ++i) {
		gs[i] = hamiltonian(i, 0);
	}
}

double computeOneSigma(const ParamsType& params)
{
	std::vector<double> gs;
	buildGroundState(gs, params);
	double sigma = measureSigma(gs, params);
	return sigma;
}

std::vector<double> computeSigmas(ParamsType& params, double m_initial, unsigned int m_total, double m_step)
{
	std::vector<double> sigma(m_total);
	for (unsigned int i = 0; i < m_total; ++i) {
		double mass = m_initial + i*m_step;
		params.set("mass", mass);
		sigma[i] = computeOneSigma(params);
		std::cout<<mass<<" "<<sigma[i]<<"\n";
	}

	return sigma;
}

int main(int argc, char* argv[])
{
	if (argc != 5) {
		std::cerr<<"USAGE: "<<argv[0]<<" input_file m_initial m_total m_step\n";
		return 1;
	}

	ParamsType params(argv[1]);
	double m_initial = std::atof(argv[2]);
	double m_total = std::atoi(argv[3]);
	double m_step = std::atof(argv[4]);
	std::vector<double> sigmas = computeSigmas(params, m_initial, m_total, m_step);
}
