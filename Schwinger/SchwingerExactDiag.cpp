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

// Aliases
//! Dense matrix from PsimagLite
using MatrixType = PsimagLite::Matrix<double>;

//! Sparse matrix (CRS) from PsimagLite
using SparseMatrixType = PsimagLite::CrsMatrix<double>;

//! Sparse row class
using SparseRow = PsimagLite::SparseRow<SparseMatrixType>;

//! Vector aliases
using VectorType = std::vector<double>;
using VectorUintType = std::vector<unsigned int>;

/*! Basis of a spinless fermion with fixed number of electrons */
using BasisFermions = LanczosPlusPlus::BasisOneSpin;

/*! Table of combinatorials; static data */
PsimagLite::Matrix<SizeType> LanczosPlusPlus::BasisOneSpin::comb_;

/*! Number of sites for static data */
SizeType LanczosPlusPlus::BasisOneSpin::nsite_;

/*!
 * \brief This class contains the parameters used by a single
 *        diagonalization
 */
class Params {
public:
	/*!
	 * \brief Constructor:
	 *
	 *  This constructor reads an input file and stores
	 *  the parameters found in it.
	 *  The format of the file is
	 *  type label value
	 *  where type can be real, integer, or special.
	 *  The special type sets parameters to a special
	 *  predetermined value and does not take value from the input.
	 *
	 * This ctor also allocates the basis of fermions with
	 * a fixed nsites and number of electrons called npart,
	 * where 0<=npart<=nsites
	 *
	 * \param[in] filename The filename for the input file
	 */
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

	/*!
	 * \brief Destructor
	 *
	 * This dtor deallocates the basis and sets it to the nullptr
	 * to help debugging in case of memory bugs.
	 */
	~Params()
	{
		delete basis_fermions;
		basis_fermions = nullptr;
	}

	/*!
	 * \brief Sets a value for an uint parameter
	 *
	 * \param[in] label The name of the parameter
	 * \param[in] value The value of the parameter
	 */
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

	/*! \brief Sets a value for a parameter of type double
	 *
	 * \param[in] label The name of the parameter
	 * \param[in] value The value of the parameter
	 */
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

	/*! \brief Sets a parameter to a special value
	 *
	 * This helps the user so that values for predetermined parameters
	 * do not need to be specified, but still need to be listed in the
	 * input file.
	 *
	 * \param[in] label The name of the parameter
	 */
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

	/*! \brief Returns or gets the value of a parameter
	 *
	 * All parameters are returned as double, even if they
	 * are integers. The caller should convert appropriately.
	 *
	 * \param[in] label The name of the wanted parameter
	 *
	 * \returns The value of the parameter as a double
	 */
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

	/*! \brief Returns the fermionic basis object
	 *
	 * \returns A constant reference to the fermionic basis
	 */
	const BasisFermions& basisFermions() const
	{
		assert(basis_fermions);
		return *basis_fermions;
	}

private:

	// The mass parameter
	double mass_ = 0;

	// The hopping parameter
	double t_ = 0;

	// The number of sites
	unsigned int N_ = 0;

	// The number of states for one link, or n in Z(n)
	unsigned int n_ = 0;

	// Whether to use lanczos or exact diag
	unsigned int use_lanczos = 0;

	// The number of electrons
	unsigned int npart_ = 0;

	// A pointer to the fermionic basis object
	BasisFermions* basis_fermions = nullptr;
};

//! Alias for the Params class
using ParamsType = Params;

/*!
 * \brief Returns the maximum field for a single link in Z(n)
 *
 * Example, n = 3, then maximum is 1, because the set is {-1, 0, 1}
 * If n = 4, then the maximum is 2, because the set is {-2, -1, 1, 2}
 *
 * \param[in] n The number of states for a single link
 *
 * \returns The maximum field value
 */
int maxField(unsigned int n)
{
	return (n & 1) ? (n - 1)/2 : n/2;
}

/*!
 * \brief Returns the field value for a given index
 *
 * Example: n = 3, then the field for 0 is -1, the field for 1 is 0,
 * and the field for 2 is 1, given that the set of fields is
 * {-1, 0, 1}
 *
 * \param[in] ind The index of the field wanted
 * \param[in] n The number of states of a single link
 * \returns The value of the field
 */
int indexToField(int ind, unsigned int n)
{
	return ind - maxField(n);
}

/*!
 * \brief Returns the index for a given field
 *
 * Example: n = 3, then the index for field -1 is 0,
 * the index for field 0 is 1, and the index for field 1 is 2.
 *
 * \param[in] l The value of the field
 * \param[in] n The number of states of a single link
 * \returns The index of the field
 */
unsigned int fieldToIndex(int l, unsigned int n)
{
	int index = l + maxField(n);
	assert(index >= 0 && index < static_cast<int>(n));
	return index;
}

/*!
 * \brief Determines whether a field value is valid
 *
 * \param[in] l The value of the field
 * \param[in] n The number of states of a single link
 * \returns true if the field is in range, false otherwise
 */
bool isInRange(int l, unsigned int n)
{
	int max = maxField(n);
	return (l >= -max && l <= max);
}

/*!
 * \brief Adds to field values with wrapping
 *
 * \param[in] prev The first field value
 * \param[in] val The second field value
 * \param[in] n The number of states in a single link
 *
 * \returns The sum of the fields with wrapping
 */
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

/*!
 * \brief Computes a field value by Gauss law
 *
 * This function uses Gauss law to compute the field
 * at the link following site site, given the charge charge
 * at site site, and the field at the previous link
 *
 * \param[in] prev The field at the previous link
 * \param[in] charge The charge
 * \param[in] site The site
 *
 * \return The value of the field
 */
int computeL(int prev, int charge, unsigned int site, unsigned int n)
{
	int val = (site & 1) ? charge - 1 : charge;
	return computeLanySite(prev, val, n);
}

/*!
 * \brief Computes Sigma in a given state
 *
 * This function computes and returns
 * \[
 * \Sigma = \frac 1N \langle gs | \sum_i E_i | gs \rangle
 * \]
 *
 * \param[in] gs The state vector in the computational basis
 * \param[in] params The params object
 *
 * \return The value of \Sigma
 */
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

/*!
 * \brief Computes the Hamiltonian for one value of the border field
 *
 * \param[in/out] m The Hamiltonian matrix
 * \param[in/out] counter A counter for the sparse matrix format
 * \param[in] d The vector of precomputed diagonal values
 * \param[in] ind The index of the border field
 * \param[in] params The parameters object
 *
 */
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

/*!
 * \brief Adds together diagonal and non-diagonal elements of the Hamiltonian
 *
 * This functions runs a loop for all possible values of the border
 * field. The diagonal elements were computed before.
 *
 * \param[in/out] sparse The hamiltonian in CRS form
 * \param[in] diagonal The precomputed diagonal elements
 * \param[in] params The parameters
 */
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

/*!
 * \brief Computes the Hamiltonian values of the staggered mass term
 *
 * \param[in] params The parameters
 *
 * \returns The vector of diagonal elements for the mass term
 */
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

/*!
 * \brief Compute diagonal elements for one border field value
 *
 * \param[out] dvector The vector to store the diagonal values
 * \param[in] dvector0 The vector containing the mass term diagonals
 * \param[in] ind The index of the border field
 * \param[in] n The number of states of a single link
 */
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

/*!
 * \brief Compute all diagonal elements of the Hamiltonian
 *
 * \param[in] params The parameters
 *
 * \returns The vector of diagonals
 */
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

/*!
 * \brief Build the Hamiltonian in sparse matrix form
 *
 * This function first pre-computes the diagonals, and then
 * it computes the off-diagonals and adds the diagonals to the full
 * Hamiltonian matrix.
 *
 * \param[in] msparse[out] The Hamiltonian in CRS form
 * \param[in] params[in] The parameters
 */
void buildHamiltonian(SparseMatrixType& msparse, const ParamsType& params)
{
	VectorType diagonal = computeDiagonal(params);
	addNonDiagonal(msparse, diagonal, params);
}

/*!
 * \brief Compute the ground state vector using Lanczos
 *
 * \param[in] msparse The CRS Hamiltonian matrix
 * \param[out] gs The ground state vector
 * \param[in] params The parameters
 */
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

/*!
 * \brief Compute the ground state vector with dense diagonalization
 *
 * \param[out] eigs The eigenvalues of the Hamiltonian
 * \param[in/out] m On input, the Hamiltonian matrix in dense form.
 *                  On output, all the eigenvectors of the Hamiltonian
 * \param[out] gs The lowest eigenvector of the Hamiltonian
 */
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

/*!
 * \brief Compute the ground state and its energy
 *
 * \param[out] energies [Only with exact diag.] All ground state energies
 * \param[out] gs The ground state vector
 * \param[in] params The parameters
 */
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

/*!
 * \brief Compute Hamiltonian properties for a range of masses
 *
 * Compute $\Sigma$ and the ground state vector. If using exact diag.,
 * compute also all eigenvalues, and the first energy gap.
 * This function prints to the standard output three columns:
 * the mass, the value of $\Sigma$, the first energy gap.
 * If using Lanzos the first energy gap appears as zero in all cases.
 *
 * \param[in] params The parameters
 * \param[in] m_initial The initial mass
 * \param[in] m_total The total number of mass values
 * \param[in] m_step The step in mass
 *
 */
void computeThings(ParamsType& params, double m_initial, unsigned int m_total, double m_step)
{
	std::vector<double> gs;
	std::vector<double> energies;
	for (unsigned int i = 0; i < m_total; ++i) {
		double mass = m_initial + i*m_step;
		params.set("mass", mass);
		buildGroundState(energies, gs, params);
		double sigma= measureSigma(gs, params);
		double gap = (energies.size() > 1) ? energies[1] - energies[0] : 0;
		std::cout<<mass<<" "<<sigma<<" "<<gap<<"\n";
	}
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
	computeThings(params, m_initial, m_total, m_step);
}
