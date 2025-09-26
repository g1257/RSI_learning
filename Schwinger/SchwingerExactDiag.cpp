#include <string>
#include <vector>
#include <fstream>

#include "CrsMatrix.h"
#include "SparseRow.h"
#include "LanczosSolver.h"
#include "Matrix.h"
#include "PsimagLite.h"
#include "Models/HubbardOneOrbital/BasisOneSpin.h"
#include "BitManip.h"

using MatrixType = PsimagLite::Matrix<double>;
using SparseMatrixType = PsimagLite::CrsMatrix<double>;
using SparseRow = PsimagLite::SparseRow<SparseMatrixType>;
using VectorType = std::vector<double>;
using VectorUintType = std::vector<unsigned int>;
using BasisFermions = LanczosPlusPlus::BasisOneSpin;

PsimagLite::Matrix<SizeType> LanczosPlusPlus::BasisOneSpin::comb_;
SizeType LanczosPlusPlus::BasisOneSpin::nsite_;

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
			} else if (stype == "special") {
				setSpecial(label);
			} else {
				fin.close();
				err("Unknown type " + stype + " must be real or integer\n");
			}
		}

		basis_fermions = new BasisFermions(N_, this->get("npart"));
		std::cerr<<N_<<" "<<n_<<" "<<t_<<" "<<mass_<<"\n";
	}

	void set(const std::string& label, unsigned int value)
	{
		if (label == "N") {
			N_ = value;
		} else if (label == "n") {
			n_ = value;
		} else if (label == "use_lanczos") {
			use_lanczos = value;
		} else if (label == "npart") {
			npart_ = value;
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

	void setSpecial(const std::string& label)
	{
		if (label == "t") {
			if (n_ == 0) {
				err("special t requires n to be set before\n");
			}

			t_ = 2.*M_PI/n_;
		} else if (label == "npart") {
			if (N_ == 0) {
				err("special npart needs N to be set before\n");
			}

			npart_ = N_/2;
		} else {
			 err("Params: unknown param or special " + label + "\n");
		}
	}

	double get(const std::string& label) const
	{
		if (label == "N") {
			return N_;
		} else if (label == "n") {
			return n_;
		} else if (label == "mass") {
			return mass_;
		} else if (label == "t") {
			return t_;
		} else if (label == "use_lanczos") {
			return use_lanczos;
		} else if (label == "npart") {
			return npart_;
		} else {
			err("Params::get(): unknown param " + label + "\n");
		}

		return 0;
	}

	const BasisFermions& basisFermions() const
	{
		assert(basis_fermions);
		return *basis_fermions;
	}

private:

	double mass_ = 0;
	double t_ = 0;
	unsigned int N_ = 0;
	unsigned int n_ = 0;
	unsigned int use_lanczos = 0;
	unsigned int npart_ = 0;
	BasisFermions* basis_fermions = nullptr;
};

using ParamsType = Params;

int maxField(unsigned int n)
{
	return (n & 1) ? (n - 1)/2 : n/2;
}

int indexToField(int ind, unsigned int n)
{
	return ind - maxField(n);
}

unsigned int fieldToIndex(int l, unsigned int n)
{
	int index = l + maxField(n);
	assert(index >= 0 && index < static_cast<int>(n));
	return index;
}

bool isInRange(int l, unsigned int n)
{
	int max = maxField(n);
	return (l >= -max && l <= max);
}

int computeLanySite(int prev, int val, unsigned int n)
{
	int l = prev + val;

	int max = maxField(n);

	if (l > max) {
		l = -max;
	}

	if (l < -max) {
		l = max;
	}

	assert(isInRange(l, n));
	return l;
}

int computeL(int prev, int charge, unsigned int site, unsigned int n)
{
	int val = (site & 1) ? charge - 1 : charge;
	return computeLanySite(prev, val, n);
}

double measureSigma(const VectorType& gs, const ParamsType& params)
{
	unsigned int N = params.get("N");
	unsigned int n = params.get("n");
	const BasisFermions& basis_fermions = params.basisFermions();
	unsigned int hilbert_fermions = basis_fermions.size();
	double sum1 = 0;
	for (unsigned int lindex = 0; lindex < n; ++lindex) {
		int field = indexToField(lindex, n);
		for (unsigned int i = 0; i < hilbert_fermions; ++i) {
			double gs_val = gs[i + lindex*hilbert_fermions];
			int prev = field;
			double sum2 = 0;
			unsigned int state_fermions = basis_fermions[i];
			for (unsigned int site = 0; site < N; ++site) {
				unsigned int mask1 = (1 << site);
				unsigned int c1 = (state_fermions & mask1) ? 1 : 0; // charge at site
				int l1 = computeL(prev, c1, site, n);
				sum2 += l1;
				prev = l1;
			}

			sum1 += sum2*gs_val*gs_val; // times complex conjugate in reality
		}
	}

	return sum1/N;
}

void addNonDiagonalOneL(SparseMatrixType& m, unsigned int& counter, const VectorType& d, unsigned int ind, const ParamsType& params)
{
	unsigned int N = params.get("N");
	unsigned int n = params.get("n");
	double t = params.get("t");
	unsigned int hilbert_fermions = params.basisFermions().size();
	double factor = -t*n*0.5/M_PI;
	int field = indexToField(ind, n);

	assert(hilbert_fermions*n == m.rows());
	assert(hilbert == d.size());

	for (unsigned i = 0; i < hilbert_fermions; ++i) {
		unsigned int row = ind*hilbert_fermions + i;

		m.setRow(row, counter);
		SparseRow sparse_row;
		// Set the diagonal element computed elsewhere
		sparse_row.add(row, d[row]);

		unsigned int state_f = params.basisFermions()[i];
		int prev = field;
		// Open BC
		for (unsigned int site = 0; site < N - 1; ++site) {
			unsigned int mask1 = (1 << site);
			unsigned int mask2 = (1 << (site + 1));
			unsigned int c1 = (state_f & mask1) ? 1 : 0; // charge at site
			unsigned int c2 = (state_f & mask2) ? 1 : 0; // charge at site + 1
			int l1 = computeL(prev, c1, site, n);
			if (c1 == c2) {
				prev = l1;
				continue;
			}

			unsigned int j = (state_f ^ (mask1 | mask2));

			int l2 = computeLanySite(field, c2 - c1, n);

			unsigned int l2index = fieldToIndex(l2, n);
			unsigned int jindex = params.basisFermions().perfectIndex(j);
			assert(jindex < hilbert_fermions);
			unsigned int col = jindex + hilbert_fermions*l2index;
			assert(col < m.cols());
			assert(col != row);
			// No fermion sign because it's one dimensional with OBC
			sparse_row.add(col, factor);

			prev = l1;
		}

		// finish the row
		counter += sparse_row.finalize(m);
	}
}

void addNonDiagonal(SparseMatrixType& msparse, const VectorType& diagonal, const ParamsType& params)
{
	unsigned int n = params.get("n");
	unsigned int hilbert = diagonal.size();
	msparse.resize(hilbert, hilbert);
	unsigned int counter = 0; // for sparse matrix rowptr
	for (unsigned i = 0; i < n; ++i) {
		addNonDiagonalOneL(msparse, counter, diagonal, i, params);
	}

	msparse.setRow(hilbert, counter);
	msparse.checkValidity();
}

VectorType buildDiagonalNoField(const ParamsType& params)
{
	unsigned int N = params.get("N");
	unsigned int n = params.get("n");
	double m = params.get("mass");
	double factor = m*n*0.5/M_PI;
	unsigned int hilbert_f = params.basisFermions().size();
	VectorType dvector(hilbert_f);

	for (unsigned i = 0; i < hilbert_f; ++i) {
		unsigned int state_f = params.basisFermions()[i];
		double sum = 0;
		for (unsigned int site = 0; site < N; ++site) {
			double sign = (site & 1) ? -1 : 1;
			unsigned int mask = (1 << site);
			if (mask & state_f) {
				sum += sign;
			}
		}

		dvector[i] = sum*factor;
	}

	return dvector;
}

void setDiagonal(VectorType& dvector, const VectorType& dvector0, unsigned int ind, unsigned int n)
{
	double etilde = indexToField(ind, n);
	const double e_squared = etilde * etilde;
	unsigned int hilbert_f = dvector0.size();
	unsigned int offset = ind*hilbert_f;
	for (unsigned i = 0; i < hilbert_f; ++i) {
		dvector[i + offset] = dvector0[i] + e_squared;
	}
}

VectorType computeDiagonal(const ParamsType& params)
{
	unsigned int n = params.get("n");
	unsigned int hilbert_f = params.basisFermions().size();
	unsigned int hilbert = n*hilbert_f;

	VectorType dvector0 = buildDiagonalNoField(params);

	VectorType dvector(hilbert);
	for (unsigned int i = 1; i < n; ++i) {
		setDiagonal(dvector, dvector0, i, n);
	}

	return dvector;
}

void buildHamiltonian(SparseMatrixType& msparse, const ParamsType& params)
{
	VectorType diagonal = computeDiagonal(params);
	addNonDiagonal(msparse, diagonal, params);
}

void buildGroundStateLanczos(const SparseMatrixType& msparse, VectorType& gs, const ParamsType& params)
{
	unsigned int hilbert = msparse.rows();
	assert(msparse.cols() == hilbert);

	using SolverParametersType = PsimagLite::ParametersForSolver<double>;

    SolverParametersType params_lanczos;
    params_lanczos.lotaMemory = true;
	params_lanczos.tolerance = 1e-40;
	params_lanczos.options = "reortho";
	params_lanczos.minSteps = 20;

    PsimagLite::LanczosSolver<SolverParametersType,
                              PsimagLite::CrsMatrix<double>,
                              VectorType>
        lanczosSolver(msparse, params_lanczos);

    double e = 0;
	gs.resize(hilbert);
	std::fill(gs.begin(), gs.end(), 0);
    VectorType initial(hilbert);
    PsimagLite::fillRandom(initial);
    lanczosSolver.computeOneState(e, gs, initial, 0);
}

void buildGroundStateEd(VectorType& eigs, MatrixType& m, VectorType& gs)
{
	unsigned int hilbert = m.rows();
	eigs.resize(hilbert);
	diag(m, eigs, 'V');


	gs.resize(hilbert);
	for (unsigned int i = 0; i < hilbert; ++i) {
		gs[i] = m(i, 0);
	}
}

void buildGroundState(VectorType& energies, VectorType& gs, const ParamsType& params)
{
	SparseMatrixType hamiltonian;
	buildHamiltonian(hamiltonian, params);

	unsigned int hilbert = hamiltonian.rows();
	assert(hilbert == hamiltonian.cols());

	if (hilbert < 4097 && !isHermitian(hamiltonian.toDense())) {
		if (hamiltonian.rows() < 40) {
			std::cout<<hamiltonian.toDense();
		}

		err("Ham not hermitian\n");
	}

	bool use_lanczos = (params.get("use_lanczos") > 0);
	if (use_lanczos) {
		buildGroundStateLanczos(hamiltonian, gs, params);
	} else {
		if (hilbert > 4096) {
			err("Hilbert size to big for exact diag; set use_lanczos in your input file instead\n");
		}

		MatrixType denseH = hamiltonian.toDense();
		// denseH will become eigenvectors
		buildGroundStateEd(energies, denseH, gs);
	}
}

std::vector<double> computeThings(ParamsType& params, double m_initial, unsigned int m_total, double m_step)
{
	std::vector<double> sigma(m_total);
	std::vector<double> gs;
	std::vector<double> energies;
	for (unsigned int i = 0; i < m_total; ++i) {
		double mass = m_initial + i*m_step;
		params.set("mass", mass);
		buildGroundState(energies, gs, params);
		sigma[i] = measureSigma(gs, params);
		double gap = (energies.size() > 1) ? energies[1] - energies[0] : 0;
		std::cout<<mass<<" "<<sigma[i]<<" "<<gap<<"\n";
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
	unsigned int m_total = std::atoi(argv[3]);
	double m_step = std::atof(argv[4]);
	std::vector<double> sigmas = computeThings(params, m_initial, m_total, m_step);
}
