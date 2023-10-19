#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
using namespace std;

static default_random_engine generator;
#define M_PI 3.14159265358979323846


double W(int N, double T, double r, double q, double sigma, double S0) {
	normal_distribution<double> dist(0, 1);
	double dt = T / N;
	vector<double> S;
	S.push_back(S0);
	double drift = exp(dt * ((r - q) - 0.5 * pow(sigma, 2)));
	double vol = sqrt(pow(sigma, 2) * dt);
	for (int i = 1; i < N; i++) {
		double Z = dist(generator);
		S.push_back(S[i - 1] * drift * exp(vol * Z));
	}
	
	return S.back();

}



/*
* @brief Monte Carlo Simulation for Call Option Pricing with Standard Error Measure
* @param underlying_price: price of stock/underlying asset
* @param strike_price: strike price of contract
* @param risk_free_interest_rate: risk-free interest rate
* @param time_to_maturity: time until expiration date of contract (assumed contract can only be exercised then)
* @param volatility: volatility of underlying asset
* @param dividend_yield: dividend yield for instruments that pay in dividends
* @param n: number of time steps
* @param m: number of simulations/random variables
* @param stderrorp: optional pointer to store standard error into
* @returns call option price based on Monte Carlo Simulation
*/
double callOptionValue(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility, double dividend_yield, int N, int M, double* stderrorp)
{

	// standard error placeholders
	double sum_CT = 0.0;

	vector<double> simulations;

	// run simulations
	for (int i = 0; i < M; i++) {
		double ST = W(N, time_to_maturity, risk_free_interest_rate, dividend_yield, volatility, underlying_price);
		simulations.push_back(ST);
		double CT = max(0.0, ST - strike_price);
		sum_CT += CT;
	}

	double C0 = exp(-risk_free_interest_rate * time_to_maturity) * sum_CT / M;
	if (stderrorp) {
		double sum_CT2 = 0.0;
		// calculating sample variance
		for (double x : simulations) {
			sum_CT2 = pow(x - C0, 2);
		}
		double SE = sum_CT2 / (M-1);
		*stderrorp = SE;
	}
	return C0;
}

/*
* @brief Monte Carlo Simulation for Put Option Pricing with Standard Error Measure
* @param underlying_price: price of stock/underlying asset
* @param strike_price: strike price of contract
* @param risk_free_interest_rate: risk-free interest rate
* @param time_to_maturity: time until expiration date of contract (assumed contract can only be exercised then)
* @param volatility: volatility of underlying asset
* @param dividend_yield: dividend yield for instruments that pay in dividends
* @param n: number of time steps
* @param m: number of simulations/random variables
* @param stderrorp: optional pointer to store standard error into
* @returns call option price based on Monte Carlo Simulation
*/
double putOptionValue(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility, double dividend_yield, int N, int M, double* stderrorp) 
{
	// standard error placeholders
	double sum_CT = 0.0;

	vector<double> simulations;

	// run simulations
	for (int i = 0; i < M; i++) {
		double ST = W(N, time_to_maturity, risk_free_interest_rate, dividend_yield, volatility, underlying_price);
		simulations.push_back(ST);
		double CT = max(0.0, strike_price - ST);
		sum_CT += CT;
	}

	double C0 = exp(-risk_free_interest_rate * time_to_maturity) * sum_CT / M;
	if (stderrorp) {
		double sum_CT2 = 0.0;
		// calculating sample variance
		for (double x : simulations) {
			sum_CT2 = pow(x - C0, 2);
		}
		double SE = sum_CT2 / (M - 1);
		*stderrorp = SE;
	}
	return C0;

}


// Use sieve of Eratosthenes to compute first n primes
vector<int> sieve(int n) {
	vector<bool>prime(n*n, true);
	vector<int> res;
	int p = 2;
	for (int i = 0; i <= n; i++) {
		while (!prime[p]) p++;
		res.push_back(p);
		for (int j = p; j < n * n; j += p) {
			prime[j] = false;
		}
	}
	return res;

}


double scrambledRadicalInverse(vector<int> perm, int a, int b) {
	double reversedDigits = 0;

	double invBase = 1.0 / (double) b, invBaseN = 1.0;

	while (a > 0) {
		// Integer divide to determine how many times base 'b' goes into 'a'
		int next = a / b;
		// Get integer remainder, determines next digit in sequence
		int digit = a - next * b;
		// Extend sequence to new digit
		reversedDigits = (double) reversedDigits * b + perm[digit];
		// Update the power
		invBaseN = invBase;
		a = next;
	}
	return reversedDigits * invBaseN;
}

double generate(vector<int>&primes, vector<vector<int>> &perms, int currentDimension, int maxDimension, int globalSample) {
	if (currentDimension >= maxDimension) {
		default_random_engine generator;
		uniform_real_distribution<double> distribution(0.0, 1.0);
		return distribution(generator);
	}

	int base = primes[currentDimension + 1];
	return scrambledRadicalInverse(perms[currentDimension], globalSample, base);

}

double boxMuller(double U1, double U2, double mean, double sigma) {
	double u = U1 * 2 * M_PI;
	return sigma * sqrt(U2) * cos(u) + mean;
}

double haltonPath(int N, double T, double r, double q, double sigma, double S0, vector<int> primes, vector<vector<int>> &perms) {
	// dimensions of sampling
	int currDim = 0, globalSample = 0;

	double dt = T / N;
	vector<double> S;
	S.push_back(S0);
	double drift = exp(dt * ((r - q) - 0.5 * pow(sigma, 2)));
	double vol = sqrt(pow(sigma, 2) * dt);
	for (int i = 1; i < N; i++) {
		double Z1 = generate(primes, perms, currDim, N, globalSample);
		globalSample++;
		double Z2 = generate(primes, perms, currDim, N, globalSample);
		currDim++;
		globalSample++;
		double Z = boxMuller(Z1, Z2, 0, 1);
		S.push_back(S[i - 1] * drift * exp(vol * Z));
	}

	return S.back();

}

/*
* @brief Monte Carlo Simulation Using Halton Sampling for Call Option Pricing with Standard Error Measure
* @param underlying_price: price of stock/underlying asset
* @param strike_price: strike price of contract
* @param risk_free_interest_rate: risk-free interest rate
* @param time_to_maturity: time until expiration date of contract (assumed contract can only be exercised then)
* @param volatility: volatility of underlying asset
* @param dividend_yield: dividend yield for instruments that pay in dividends
* @param N: number of time steps
* @param M: number of simulations/random variables
* @param stderrorp: optional pointer to store standard error into
* @returns call option price based on Monte Carlo Simulation
*/
double callOptionValueHalton(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility, double dividend_yield, int N, int M, double* stderrorp)
{
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	// standard error placeholders
	double sum_CT = 0.0;

	vector<double> simulations;


	vector<int> primes = sieve(N + 1);

	// Pre-compute a permutation for every prime base b
	vector<vector<int>> perms;
	for (int i = 0; i < N; i++) {
		vector<int> row;
		for (int j = 0; j < primes[i + 1]; j++) {
			row.push_back(j);
		}
		shuffle(row.begin() + 1, row.end(), default_random_engine(seed));
		perms.push_back(row);
	}

	// run simulations
	for (int i = 0; i < M; i++) {
		double ST = haltonPath(N, time_to_maturity, risk_free_interest_rate, dividend_yield, volatility, underlying_price, primes, perms);
		simulations.push_back(ST);
		double CT = max(0.0, ST- strike_price);
		sum_CT += CT;
	}

	double C0 = exp(-risk_free_interest_rate * time_to_maturity) * sum_CT / M;
	if (stderrorp) {
		double sum_CT2 = 0.0;
		// calculating sample variance
		for (double x : simulations) {
			sum_CT2 = pow(x - C0, 2);
		}
		double SE = sum_CT2 / (M - 1);
		*stderrorp = SE;
	}
	return C0;
}


/*
* @brief Monte Carlo Simulation Using Halton Sampling for Put Option Pricing with Standard Error Measure
* @param underlying_price: price of stock/underlying asset
* @param strike_price: strike price of contract
* @param risk_free_interest_rate: risk-free interest rate
* @param time_to_maturity: time until expiration date of contract (assumed contract can only be exercised then)
* @param volatility: volatility of underlying asset
* @param dividend_yield: dividend yield for instruments that pay in dividends
* @param N: number of time steps
* @param M: number of simulations/random variables
* @param stderrorp: optional pointer to store standard error into
* @returns put option price based on Monte Carlo Simulation
*/
double putOptionValueHalton(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility, double dividend_yield, int N, int M, double* stderrorp)
{
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	// standard error placeholders
	double sum_CT = 0.0;

	vector<double> simulations;


	vector<int> primes = sieve(N + 1);

	// Pre-compute a permutation for every prime base b
	vector<vector<int>> perms;
	for (int i = 0; i < N; i++) {
		vector<int> row;
		for (int j = 0; j < primes[i + 1]; j++) {
			row.push_back(j);
		}
		shuffle(row.begin() + 1, row.end(), default_random_engine(seed));
		perms.push_back(row);
	}

	// run simulations
	for (int i = 0; i < M; i++) {
		double ST = haltonPath(N, time_to_maturity, risk_free_interest_rate, dividend_yield, volatility, underlying_price, primes, perms);
		simulations.push_back(ST);
		double CT = max(0.0, strike_price - ST);
		sum_CT += CT;
	}

	double C0 = exp(-risk_free_interest_rate * time_to_maturity) * sum_CT / M;
	if (stderrorp) {
		double sum_CT2 = 0.0;
		// calculating sample variance
		for (double x : simulations) {
			sum_CT2 = pow(x - C0, 2);
		}
		double SE = sum_CT2 / (M - 1);
		*stderrorp = SE;
	}
	return C0;
}


int main() {
	double S = 50.00;
	double K = 10.00;
	double vol = 0.5;
	double r = 0.05;
	double T = 3.5;
	int N = 100;
	int M = 10000;
	double q = 0.0;

	double* errorfactorCall = (double*) malloc(sizeof(double));
	double resultCall = callOptionValue(S, K, r, T, vol, q, N, M, errorfactorCall);
	printf("Call option value: %f \n", resultCall);
	printf("Standard error: %f \n", *errorfactorCall);

	double* errorfactorPut = (double*)malloc(sizeof(double));
	double resultPut = putOptionValue(S, K, r, T, vol, q, N, M, errorfactorPut);
	printf("Put option value: %f \n", resultPut);
	printf("Standard error: %f \n", *errorfactorPut);

	resultCall = callOptionValueHalton(S, K, r, T, vol, q, N, M, errorfactorCall);
	printf("Call option value: %f \n", resultCall);
	printf("Standard error: %f \n", *errorfactorCall);


	resultCall = putOptionValueHalton(S, K, r, T, vol, q, N, M, errorfactorPut);
	printf("Put option value: %f \n", resultCall);
	printf("Standard error: %f \n", *errorfactorPut);

	free(errorfactorCall);
	free(errorfactorPut);
	
}