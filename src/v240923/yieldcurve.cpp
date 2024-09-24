
#include "Statisticality.hpp"

#include "Statisticality/Mathematics/Optimize/OptimizeKernel.hpp"

using namespace Statisticality;

// Dirty pricing - bond, contunously
template<typename _RealType = double_t, typename _IntegerType = int64_t, typename _BoolType = boolean>
_RealType bondpricing(_RealType n, _RealType cr, _RealType par, _RealType r, _RealType f = 2,
	_BoolType continuous = 1, _RealType zero_threshold = 1e-6) noexcept
{
	// dirty price: pv - continuously
	_RealType fut = ((n * f) - floor(n * f)) / f;
	_RealType x = fut / (1.0 / f);
	if (x < zero_threshold) // threshold
	{
		x = 1.0;
	}
	_RealType pv = 0.0;

	// new whole years
	_IntegerType nyr = ceil(n * f);
	if (continuous)
	{
		for (_IntegerType i = 1; i <= nyr; ++i)
		{
			pv += (par * cr / f) * exp(-r / f * (i - 1.0 + x));
		}
		pv += par * exp(-r / f * (nyr - 1.0 + x));
	}
	else
	{
		for (_IntegerType i = 1; i <= nyr; ++i)
		{
			pv += (par * cr / f) / pow(1.0 + r / f, (i - 1.0 + x));
		}
		pv += par / pow(1.0 + r / f, (nyr - 1.0 + x));
	}
	
	return pv;
}

// Use my proud Optimize(R) to calculate the ytm
int calculate_ytm(const std::string& filepath,
				  const std::string& save_to,
				  int ncol, int crcol, int parcol, int pvcol, double f = 2) noexcept
{
	dataframe file;
	dataframe header;

	// Read from file
	if (filepath.ends_with(".xlsx"))
	{
		file.readFromExcel(filepath, 0);
	}
	else if (filepath.ends_with(".csv"))
	{
		file.readFromCsv(filepath);
	}
	// Get rid of the headers
	header = file.getdataframe(0, 0); header.resize(1, header.col() + 1); header[header.size() - 1] = "r";
	file = file.getdataframe(1);

	// Get the vector of n
	auto ns = file.getcol<double>(ncol);
	// Get the vector of cr
	auto crs = file.getcol<double>(crcol);
	// Get the vector of par
	auto pars = file.getcol<double>(parcol);
	// Get the vector of pv dirty prices
	auto pvs = file.getcol<double>(pvcol);

	// Create a vector of assumed rs (0.05 here)
	auto rs = std::vector<double>(ns.size(), 0.05);

	// Create a vector of f
	auto fs = std::vector<double>(ns.size(), f);

	// Merge them into a matrix (n, cr, par, ?, f)
	matrix<double> data(ns, ns.size(), 1);
	data.cbind(matrix<double>(crs, crs.size(), 1));
	data.cbind(matrix<double>(pars, pars.size(), 1));
	data.cbind(matrix<double>(rs, rs.size(), 1));
	data.cbind(matrix<double>(fs, fs.size(), 1));

	// For each row, try to use Optimize() to get the rs
	for (size_t i = 0; i < data.row(); ++i)
	{
		// Create a new wrapper to test the rs
		auto _wrapper = [&](const std::vector<double>& x)
		{
			return pow(bondpricing(data(i, 0ULL), data(i, 1ULL), data(i, 2ULL), x[0], data(i, 4ULL)) - pvs[i], 2);
		};
		auto _wrappera = [&](const double& x)
		{
			return pow(bondpricing(data(i, 0ULL), data(i, 1ULL), data(i, 2ULL), x, data(i, 4ULL)) - pvs[i], 2);
		};

		// Optimize it
		Optimizer<double> opt(0.001, 10000, "adam");
		auto optm = opt.autostart(_wrapper, { 0.01 }, { 0.01 }, { 20 }, 2, 1000);
		auto lr = opt.autolrate(_wrapper, optm, 0.001, 1.4, 10, 2, 1000);
		opt.set_learningRate(lr);
		auto optm2 = opt.optimize(_wrapper, optm, 1000);

		// Get and assign the params
		double optmparam = optm2[0];
		data(i, 3ULL) = optmparam;
	}

	// Pool the result and try to see it
	Pool(data.cbind(data, matrix<double>(pvs, pvs.size(), 1)), "data", "c");

	// Create a new column of dataframe and cbind it
	dataframe optm_rs(data.row(), 1);
	optm_rs.setmap(data.getcol(3ULL));
	dataframe nfile = file;
	nfile.cbind(optm_rs);
	nfile = header.rbind(header, nfile);

	// Pool the nfile object into runtime
	Pool(nfile, "nfile", "c");

	// Save the new dataframe
	nfile.writeToExcel(save_to);

	// Start Runtime

	return 0;
}
//
// calculate_ytm("C:\\Users\\xxx\\Desktop\\TreasuryData.csv",
//	"C:\\Users\\xxx\\Desktop\\TreasuryData.csv.42.xlsx",
//	13, 14, 15, 12);
// Runtime();

// Drop outliers
template<typename _RealType, typename _IterableType>
std::pair<_IterableType, _IterableType> drop_outliers(const _IterableType& detect, const _IterableType& followed,
	int previous = 5, _RealType k = 3.24) noexcept {

	// Length
	const size_t _len = detect.size();

	// Index selected
	std::vector<size_t> _selected;
	_selected.reserve(_len);

	// Deque, for iterable computation
	std::deque<_RealType> _deque;

	// Deque to vector
	auto to_vector = [](const std::deque<_RealType>&deque) {
		return std::vector<_RealType>(deque.begin(), deque.end());
	};

	// First fill
	for (int i = 0; i < previous; ++i){
		// Must been selected
		_selected.push_back(i);

		// Pool into deque
		_deque.push_back(detect[i]);
	}

	// Rolling detecting
	for (size_t i = previous; i < _len; ++i){

		// Compute mean and stdev
		_RealType _mean = mean(to_vector(_deque));
		_RealType _stdev = sd(to_vector(_deque));

		// If its > mean + k*stdev or < mean - k*stdev
		const _RealType& _node = detect[i];
		if ((_node > _mean + k * _stdev) ||
			(_node < _mean - k * _stdev)) {
			;
		}
		else {
			// Good, get it
			_selected.push_back(i);
		}

		// Update deque
		_deque.pop_front();
		_deque.push_back(detect[i]);
	}
	
	// Slicing
	std::pair<_IterableType, _IterableType> _pair;
	for (size_t i = 0; i < _selected.size(); ++i) {
		_pair.first.push_back(detect[_selected[i]]);
		_pair.second.push_back(followed[_selected[i]]);
	}
	return _pair;
}

// Use my proud Optimize(R) to fit the curve - linear model
int fit_curve_lm(const std::string& filepath,
				 const std::string& save_to) noexcept
{
	dataframe x;
	x.readFromExcel(filepath, 0);
	
	// Get rid of the headers
	x = x.getdataframe(1);

	// y: ytm, x: maturity
	auto y = x.getdataframe(-1, -1, 16, 16);
	auto pred = x.getdataframe(-1, -1, 13, 13);
	const matrix<double> yo = y.toMatrix<double>();
	const matrix<double> xo = pred.toMatrix<double>();
	matrix<double> ya = y.toMatrix<double>();
	matrix<double> xa = pred.toMatrix<double>();

	// Drop outliers
	auto xymm = drop_outliers<double, matrix<double>>(ya, xa, 5, 6.28);
	ya = xymm.first;
	xa = xymm.second;

	// Create a formula
	auto formula_x = [](const matrix<double>& single_col_x) {
		// y = a + bx + cx^2 + dx^3 + f*exp(-x) + g*log(x) + h*x*exp(-x)

		const matrix<double>& x = single_col_x;
		matrix<double> xx = x;
		xx.cbind(x.pow(2));
		xx.cbind(x.pow(3));
		xx.cbind(x.exp(-1));
		xx.cbind(x.log());
		xx.cbind(x % x.exp(-1));
		return xx;
	};
	xa = formula_x(xa);

	// Estimated params
	/*
	beta
	0.02892527 -0.00404009 0.00020215 -0.00000359 0.04664002 0.01366871 -0.00970796
	t
	63.3575 -24.1961 20.2281 -18.7217 57.1363 50.0266 -10.6833
	t-p
	-0.0000 -0.0000 -0.0000 -0.0000 -0.0000 -0.0000 -0.0000
	f-p
	0
	r2
	0.95679
	*/

	// Estimate the model
	OLS<double> ols(xa, ya, true);
	ols.estimate();

	// Print summary
	print("beta");
	ols.beta().copy().t().print(" ", 1000, 8);

	print("t");
	ols.t().copy().t().print(" ", 1000, 4);

	print("t-p");
	ols.t_p().copy().t().print(" ", 1000, 4);

	print("f-p");
	print(ols.f_p());

	print("r2");
	print(ols.r2());

	// Fetch the restored data (full, without dropping)
	matrix<double> rxxx = formula_x(xo);
	matrix<double> rest_estm_y = ols.predict(rxxx);
	matrix<double> rest_real_y = yo;
	dataframe restored = dataframe(rest_estm_y.cbind(rest_estm_y, rest_real_y));
	restored = restored.rbind(dataframe({ atom("estm"),atom("real") }, 1, 2), restored);

	// Save the fitted data
	restored.writeToExcel(save_to);

	// Create a continuous curve
	std::vector<double> _seq_x = seq<double>(0.1, 30, 0.0001);
	matrix<double> cont_estm_x = matrix<double>(_seq_x, _seq_x.size(), 1);
	matrix<double> cont_estm_y = ols.predict(formula_x(cont_estm_x));
	dataframe curve = dataframe(cont_estm_y.cbind(cont_estm_x, cont_estm_y));
	curve = curve.rbind(dataframe({ atom("t"),atom("yield") }, 1, 2), curve);

	// Save the fitted curve
	curve.writeToExcel(strreplace(save_to, ".xlsx", ".curve.xlsx"));

	// Pool and return
	Pool(restored, "x", "c");
	Runtime();
	return 0;
}
// fit_curve_lm("C:\\Users\\xxx\\Desktop\\Miniproject 2 Treasury Data Fall 2024 (09-21-2024).xls",
	//	"C:\\Users\\xxx\\Desktop\\Miniproject 2 Treasury Data Fall 2024 (09-21-2024).fit.xlsx");

// Use my proud Optimize(R) to fit the curve - non linear model
int fit_curve_nlm(const std::string& filepath,
	const std::string& save_to) noexcept
{
	dataframe x;
	x.readFromExcel(filepath, 0);
	
	// Get rid of the headers
	x = x.getdataframe(1);

	// y: ytm, x: maturity
	auto y = x.getdataframe(-1, -1, 16, 16);
	auto pred = x.getdataframe(-1, -1, 13, 13);
	const matrix<double> yo = y.toMatrix<double>();
	const matrix<double> xo = pred.toMatrix<double>();
	matrix<double> ya = y.toMatrix<double>();
	matrix<double> xa = pred.toMatrix<double>();

	// Drop outliers
	auto xymm = drop_outliers<double, matrix<double>>(ya, xa, 5, 6.28);
	ya = xymm.first;
	xa = xymm.second;

	// Expand x before throwing into regression
	auto cbind_x = [](const matrix<double>& single_col_x, int to_col = 8) {
		const matrix<double>& x = single_col_x;
		matrix<double> xx = x;
		for (int i = 1; i < to_col; ++i) {
			xx.cbind(x);
		}
		return xx;
	};
	
	// Create a formula
	auto formula_bx = [](const matrix<double>& b, const matrix<double>& expanded_x) {
		// y = a + bx + cx^2 + dx^3 + f*exp(g*x) + h*log(x+1) + j*x*exp(k*x)
		// Totally 9 beta parameters including intercept

		// Note here, expanded x should have 7 cols used (with ones)

		const matrix<double>& x = expanded_x;
		matrix<double> xx = x;
		matrix<double> sum = xx.getmatrix(-1, -1, 0, 0); sum.zero();
		xx.setmatrix(x.getmatrix(-1, -1, 0, 0) * b[0],
			-1, -1, 0, 0); // a
		sum += xx.getmatrix(-1, -1, 0, 0);
		xx.setmatrix(x.getmatrix(-1, -1, 1, 1) * b[1],
			-1, -1, 1, 1); // bx
		sum += xx.getmatrix(-1, -1, 1, 1);
		xx.setmatrix(x.getmatrix(-1, -1, 2, 2).pow(2) * b[2],
			-1, -1, 2, 2); // cx^2
		sum += xx.getmatrix(-1, -1, 2, 2);
		xx.setmatrix(x.getmatrix(-1, -1, 3, 3).pow(3) * b[3],
			-1, -1, 3, 3); // dx^3
		sum += xx.getmatrix(-1, -1, 3, 3);
		xx.setmatrix(x.getmatrix(-1, -1, 4, 4).exp(b[5]) * b[4],
			-1, -1, 4, 4); // f*exp(g*x)
		sum += xx.getmatrix(-1, -1, 4, 4);
		xx.setmatrix(x.getmatrix(-1, -1, 5, 5).add(1).log() * b[6],
			-1, -1, 5, 5); // h*log(x+1)
		sum += xx.getmatrix(-1, -1, 5, 5);
		xx.setmatrix(x.getmatrix(-1, -1, 6, 6) % (x.getmatrix(-1, -1, 6, 6).exp(b[8])) * b[7],
			-1, -1, 6, 6); // j*x*exp(k*x)
		sum += xx.getmatrix(-1, -1, 6, 6);

		return sum;
	};

	// Multi-thread version
	auto formula_bx_mt = [](const matrix<double>& b, const matrix<double>& expanded_x) {
		// y = a + bx + cx^2 + dx^3 + f*exp(g*x) + h*log(x+1) + j*x*exp(k*x)
		// Totally 9 beta parameters including intercept

		// Note here, expanded x should have 7 cols used (with ones)

		const matrix<double>& x = expanded_x;
		matrix<double> xx = x;
		matrix<double> sum = xx.getmatrix(-1, -1, 0, 0); sum.zero();
		#pragma omp parallel for num_threads(20)
		for (long long i = 0; i < x.row(); ++i) {
			xx((size_t)i, 0ULL) = x((size_t)i, 0ULL) * b[0]; // a
			xx((size_t)i, 1ULL) = x((size_t)i, 1ULL) * b[1]; // bx
			xx((size_t)i, 2ULL) = ::pow(x((size_t)i, 2ULL), 2) * b[2]; // cx^2
			xx((size_t)i, 3ULL) = ::pow(x((size_t)i, 3ULL), 3) * b[3]; // dx^3
			xx((size_t)i, 4ULL) = ::exp(x((size_t)i, 4ULL) * b[5]) * b[4]; // f*exp(g*x)
			xx((size_t)i, 5ULL) = ::log(x((size_t)i, 5ULL) + 1) * b[6]; // h*log(x+1)
			xx((size_t)i, 6ULL) = ::exp(x((size_t)i, 6ULL) * b[8]) * x((size_t)i, 6ULL) * b[7]; // j*x*exp(k*x)
		}
		sum += xx.getmatrix(-1, -1, 0, 0);
		sum += xx.getmatrix(-1, -1, 1, 1);
		sum += xx.getmatrix(-1, -1, 2, 2);
		sum += xx.getmatrix(-1, -1, 3, 3);
		sum += xx.getmatrix(-1, -1, 4, 4);
		sum += xx.getmatrix(-1, -1, 5, 5);
		sum += xx.getmatrix(-1, -1, 6, 6);

		return sum;
	};

	// Estimated params
	/*
	beta
	0.66072348 3.34515148 0.70848577 -1.20041321 -0.52506380 -3166.46074780 -3.90684226 81.76467740 -102.38612116
	t
	90.6094 0.5840 0.3040 -1.8525 -39.9181 -1.4044 -0.6696 34.2304 -59.6664
	t-p
	-0.0000 0.5596 0.7614 0.0649 -0.0000 0.1612 0.5036 -0.0000 -0.0000
	f-p
	0
	r2
	0.975726
	*/

	// 1) Expand xa
	xa = cbind_x(xa, 8); // 8 means leaving 1 to ones

	// 2) Scalling, especially for GGDS gradient descending model
	Scaler<double> scl(1, 1e-8);
	auto xypair = scl.scalling(xa, ya);
	xa = xypair.first;
	ya = xypair.second;

	// Estimate the model
	GGDS<double> ggds(xa, ya, formula_bx_mt, {}, "mse", "adam", 0.85, 1e5);
	ggds.estimate();

	// Print summary
	print("beta");
	ggds.beta().copy().t().print(" ", 1000, 8);

	print("t");
	ggds.t().copy().t().print(" ", 1000, 4);

	print("t-p");
	ggds.t_p().copy().t().print(" ", 1000, 4);

	print("f-p");
	print(ggds.f_p());

	print("r2");
	print(ggds.r2());

	// Fetch the restored data (full, without dropping)
	matrix<double> rest_estm_y = scl.restoring(matrix<double>(), 
		ggds.predict(scl.scalling(cbind_x(xo, 8)).first)).second;
	matrix<double> rest_real_y = yo;
	dataframe restored = dataframe(rest_estm_y.cbind(rest_estm_y, rest_real_y));
	restored = restored.rbind(dataframe({ atom("estm"),atom("real") }, 1, 2), restored);

	// Save the fitted data
	restored.writeToExcel(save_to);

	// Create a continuous curve
	std::vector<double> _seq_x = seq<double>(0.1, 30, 0.0001);
	matrix<double> cont_estm_x = matrix<double>(_seq_x, _seq_x.size(), 1);
	matrix<double> cont_estm_y = scl.restoring(matrix<double>(),
		ggds.predict(scl.scalling(cbind_x(cont_estm_x, 8)).first)).second;
	dataframe curve = dataframe(cont_estm_y.cbind(cont_estm_x, cont_estm_y));
	curve = curve.rbind(dataframe({ atom("t"),atom("yield") }, 1, 2), curve);

	// Save the fitted curve
	curve.writeToExcel(strreplace(save_to, ".xlsx", ".curve.xlsx"));

	// Pool and return
	Pool(restored, "x", "c");
	Runtime();
	return 0;
}
// fit_curve_nlm("C:\\Users\\xxx\\Desktop\\Miniproject 2 Treasury Data Fall 2024 (09-21-2024).xls",
// 	"C:\\Users\\xxx\\Desktop\\Miniproject 2 Treasury Data Fall 2024 (09-21-2024).non-linear.fit.xlsx");

// Use my proud Optimize(R) to fit the curve - psudo neural network model
int fit_curve_pnn(const std::string& filepath,
	const std::string& save_to) noexcept
{
	dataframe x;
	x.readFromExcel(filepath, 0);

	// Get rid of the headers
	x = x.getdataframe(1);

	// y: ytm, x: maturity
	auto y = x.getdataframe(-1, -1, 16, 16);
	auto pred = x.getdataframe(-1, -1, 13, 13);
	const matrix<double> yo = y.toMatrix<double>();
	const matrix<double> xo = pred.toMatrix<double>();
	matrix<double> ya = y.toMatrix<double>();
	matrix<double> xa = pred.toMatrix<double>();

	// Drop outliers
	auto xymm = drop_outliers<double, matrix<double>>(ya, xa, 5, 6.28);
	ya = xymm.first;
	xa = xymm.second;

	// Kernel function - Hyperbolic Tangent (tanh) activation function
	auto tanh = [](double x) -> double {
		return std::tanh(x);
	};

	// Kernel function - Sigmoid activation function
	auto sigmoid = [](double x) -> double {
		return 1.0 / (1.0 + std::exp(-x));
	};

	// Kernel function - Rectified Linear Unit (ReLU) activation function
	auto relu = [](double x) -> double {
		return x > 0.0 ? x : 0.0;
	};

	// Kernel function - Leaky Rectified Linear Unit (Leaky ReLU) activation function
	auto leaky_relu = [](double x, double alpha = 0.01) -> double {
		return x > 0 ? x : alpha * x;
	};

	// Expand x before throwing into regression
	auto cbind_x = [](const matrix<double>& single_col_x, int to_col = 8) {
		const matrix<double>& x = single_col_x;
		matrix<double> xx = x;
		for (int i = 1; i < to_col; ++i) {
			xx.cbind(x);
		}
		return xx;
	};

	// Create a formula - Multi-thread version - Psudo Neural Network
	auto formula_bxnn_mt = [tanh, sigmoid](const matrix<double>& b, const matrix<double>& expanded_x) {
		
		// Layer 1
		// y = a + bx + cx^2 + dx^3 + f*exp(g*x) + h*log(x+1) + j*x*exp(k*x)
		// Totally 9 beta parameters including intercept
		//
		// Note here, expanded x should have 7 cols used (with ones)
		// 
		// Parameters number: 9

		// Activation tanh

		// Layer 2
		// Dense, with 7 nodes and 1 another constant parameter (8 total) connect with each other
		//        to reduce to 3 nodes
		// Parameters number: 24

		// Activation tanh

		// Layer 3
		// Dense, with 3 nodes and 1 another constant parameter (4 total) connect with each other
		//        to reduce to 1 nodes
		// Parameters number: 4

		// Activation sigmoid

		// Total parameters 37

		const matrix<double>& x = expanded_x;

		// Place to store nn values
		matrix<double> xx = x;
		matrix<double> xx2(x.row(), 3, 0);
		matrix<double> xx3(x.row(), 1, 0);

		// Do things in parallel
		#pragma omp parallel for num_threads(20)
		for (long long i = 0; i < x.row(); ++i) {

			// Layer 1 and activation
			xx((size_t)i, 0ULL) = tanh(x((size_t)i, 0ULL) * b[0]); // a
			xx((size_t)i, 1ULL) = tanh(x((size_t)i, 1ULL) * b[1]); // bx
			xx((size_t)i, 2ULL) = tanh(::pow(x((size_t)i, 2ULL), 2) * b[2]); // cx^2
			xx((size_t)i, 3ULL) = tanh(::pow(x((size_t)i, 3ULL), 3) * b[3]); // dx^3
			xx((size_t)i, 4ULL) = tanh(::exp(x((size_t)i, 4ULL) * b[5]) * b[4]); // f*exp(g*x)
			xx((size_t)i, 5ULL) = tanh(::log(x((size_t)i, 5ULL) + 1) * b[6]); // h*log(x+1)
			xx((size_t)i, 6ULL) = tanh(::exp(x((size_t)i, 6ULL) * b[8]) * x((size_t)i, 6ULL) * b[7]); // j*x*exp(k*x)
			
			// Layer 2 and activation
			const size_t param_index_l2 = 9;
			size_t l2_iter = 0;
			for (long long i2 = 0; i2 < 3; ++i2) {
				// Constant
				xx2((size_t)i, (size_t)i2) += b[param_index_l2 + l2_iter];
				l2_iter++;

				// Previous result, dense connection
				for (long long i2a = 0; i2a < 7; ++i2a,++l2_iter) {
					xx2((size_t)i, (size_t)i2) += b[param_index_l2 + l2_iter] * xx((size_t)i, (size_t)i2a);
				}
				
				// Activation
				xx2((size_t)i, (size_t)i2) = tanh(xx2((size_t)i, (size_t)i2));
			}

			// Layer 3 and activation
			const size_t param_index_l3 = 33;
			size_t l3_iter = 0;
			for (long long i3 = 0; i3 < 1; ++i3) {
				// Constant
				xx3((size_t)i, (size_t)i3) += b[param_index_l3 + l3_iter];
				l3_iter++;

				// Previous result, dense connection
				for (long long i3a = 0; i3a < 3; ++i3a, ++l3_iter) {
					xx3((size_t)i, (size_t)i3) += b[param_index_l3 + l3_iter] * xx2((size_t)i, (size_t)i3a);
				}

				// Activation
				xx3((size_t)i, (size_t)i3) = sigmoid(xx3((size_t)i, (size_t)i3));
			}
		}

		return xx3;
	};

	// Params:
	/*
	beta
	4.23869007 -10.43947147 9.95131868 1.98609632 4.37777135 0.51907529 1.26771631 -129.22943469 -7.18800342 1.46517368 -1.19966309 -0.61071880 2.32607774 1.76528509 -0.85463308 -7.05265164 5.89780774 -1.43926621 1.81180054 -1.33979122 -0.22757587 -2.33636955 1.04044856 1.07958929 1.25557265 -1.76374109 -2.65189505 1.33079007 -1.60269809 0.60683786 -2.48069084 3.81285243 -7.15219061 -3.32957726 -6.94216200 -4.64102974 -1.92914311
	t
	0.0000 -0.0451 0.0166 0.0501 0.0000 0.0000 0.0005 -0.0923 -0.2829 0.0000 -0.0000 -0.0000 0.0001 0.0009 -0.0000 -0.0030 0.0032 -0.0000 0.0000 -0.0356 -0.0015 -0.0066 0.0000 0.0006 0.0585 -0.0000 -0.0000 0.0163 -0.0097 0.0003 -0.0001 0.0005 -0.0563 -0.0045 -0.0094 -0.0456 -0.1381
	t-p
	1.0000 0.9640 0.9868 0.9601 1.0000 1.0000 0.9996 0.9265 0.7774 1.0000 1.0000 1.0000 0.9999 0.9993 1.0000 0.9976 0.9974 1.0000 1.0000 0.9716 0.9988 0.9948 1.0000 0.9995 0.9534 1.0000 1.0000 0.9870 0.9923 0.9997 1.0000 0.9996 0.9551 0.9964 0.9925 0.9637 0.8903
	f-p
	0
	r2
	0.9844
	*/

	// 1) Expand xa
	xa = cbind_x(xa, 8); // 8 means leaving 1 to ones

	// 2) Scalling, especially for GGDS gradient descending model
	Scaler<double> scl(1, 1e-8);
	auto xypair = scl.scalling(xa, ya);
	xa = xypair.first;
	ya = xypair.second;

	// Estimate the model
	GGDS<double> ggds(xa, ya, formula_bxnn_mt, {37, 1}, "mse", "adam", 0.05, 1e5);
	ggds.estimate();

	// Print summary
	print("beta");
	ggds.beta().copy().t().print(" ", 1000, 8);

	print("t");
	ggds.t().copy().t().print(" ", 1000, 4);

	print("t-p");
	ggds.t_p().copy().t().print(" ", 1000, 4);

	print("f-p");
	print(ggds.f_p());

	print("r2");
	print(ggds.r2());

	// Fetch the restored data (full, without dropping)
	matrix<double> rest_estm_y = scl.restoring(matrix<double>(),
		ggds.predict(scl.scalling(cbind_x(xo, 8)).first)).second;
	matrix<double> rest_real_y = yo;
	dataframe restored = dataframe(rest_estm_y.cbind(rest_estm_y, rest_real_y));
	restored = restored.rbind(dataframe({ atom("estm"),atom("real") }, 1, 2), restored);

	// Save the fitted data
	restored.writeToExcel(save_to);

	// Create a continuous curve
	std::vector<double> _seq_x = seq<double>(0.1, 30, 0.0001);
	matrix<double> cont_estm_x = matrix<double>(_seq_x, _seq_x.size(), 1);
	matrix<double> cont_estm_y = scl.restoring(matrix<double>(),
		ggds.predict(scl.scalling(cbind_x(cont_estm_x, 8)).first)).second;
	dataframe curve = dataframe(cont_estm_y.cbind(cont_estm_x, cont_estm_y));
	curve = curve.rbind(dataframe({ atom("t"),atom("yield") }, 1, 2), curve);

	// Save the fitted curve
	curve.writeToExcel(strreplace(save_to, ".xlsx", ".curve.xlsx"));

	// Pool and return
	Pool(restored, "x", "c");
	Runtime();
	return 0;
}
// fit_curve_pnn("C:\\Users\\xxx\Desktop\\Miniproject 2 Treasury Data Fall 2024 (09-21-2024).xls",
// 	"C:\\Users\\xxx\\Desktop\\Miniproject 2 Treasury Data Fall 2024 (09-21-2024).pnn-network.fit.xlsx");


int main()
{
	Init();

  // Make your own path
	fit_curve_pnn("C:\\Users\\xxx\\Desktop\\Miniproject 2 Treasury Data Fall 2024 (09-21-2024).xls",
	 	"C:\\Users\\xxx\\Desktop\\Miniproject 2 Treasury Data Fall 2024 (09-21-2024).pnn-network.fit.xlsx");
	

	Return(0);
}
