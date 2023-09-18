#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
using namespace std;

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
	default_random_engine generator;
	normal_distribution<double> distribution(0, 1);
	double dt = time_to_maturity / N;
	double nudt = (risk_free_interest_rate - 0.5 * pow(volatility, 2)) * dt;
	double volsdt = volatility * sqrt(dt);
	double lnS = log(underlying_price * exp((risk_free_interest_rate-dividend_yield) * time_to_maturity));

	// standard error placeholders
	double sum_CT = 0.0;
	double sum_CT2 = 0.0;

	// run simulations
	for (int i = 0; i < M; i++) {
		double lnSt = lnS;
		for (int j = 0; j < N; j++) {
			lnSt = lnSt + nudt + volsdt * distribution(generator);
		}
		double ST = exp(lnSt);
		double CT = max(0.0, ST - strike_price);
		sum_CT += CT;
		sum_CT2 += pow(CT, 2);
	}

	double C0 = exp(-risk_free_interest_rate * time_to_maturity) * sum_CT / M;
	if (stderrorp) {
		double sigma = sqrt((sum_CT2 - (pow(sum_CT, 2) / M)) * exp(-2 * risk_free_interest_rate * time_to_maturity) / (M - 1));
		double SE = sigma / sqrt(M);
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
	default_random_engine generator;
	normal_distribution<double> distribution(0, 1);
	double dt = time_to_maturity / N;
	double nudt = (risk_free_interest_rate - 0.5 * pow(volatility, 2)) * dt;
	double volsdt = volatility * sqrt(dt);
	double lnS = log(underlying_price * exp((risk_free_interest_rate - dividend_yield) * time_to_maturity));

	// standard error placeholders
	double sum_CT = 0.0;
	double sum_CT2 = 0.0;

	// run simulations
	for (int i = 0; i < M; i++) {
		double lnSt = lnS;
		for (int j = 0; j < N; j++) {
			lnSt = lnSt + nudt + volsdt * distribution(generator);
		}
		double ST = exp(lnSt);
		double CT = max(0.0, strike_price - ST);
		sum_CT += CT;
		sum_CT2 += pow(CT, 2);
	}

	double C0 = exp(-risk_free_interest_rate * time_to_maturity) * sum_CT / M;
	if (stderrorp) {
		double sigma = sqrt((sum_CT2 - (pow(sum_CT, 2) / M)) * exp(-2 * risk_free_interest_rate * time_to_maturity) / (M - 1));
		double SE = sigma / sqrt(M);
		*stderrorp = SE;
	}
	return C0;

}


// Use sieve of Eratosthenes to compute first n primesprimes
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
		invBaseN *= invBase;
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
	// dimensions of sampling
	int currDim = 0, globalSample = 0;

	vector<int> primes = sieve(N+1);

	// Pre-compute a permutation for every prime base b
	vector<vector<int>> perms;
	for (int i = 0; i < N; i++) {
		vector<int> row;
		for (int j = 0; j < primes[i+1]; j++) {
			row.push_back(j);
		}
		shuffle(row.begin() + 1, row.end(), default_random_engine(seed));
		perms.push_back(row);
	}

	double dt = time_to_maturity / N;
	double nudt = (risk_free_interest_rate - 0.5 * pow(volatility, 2)) * dt;
	double volsdt = volatility * sqrt(dt);
	double lnS = log(underlying_price * exp((risk_free_interest_rate - dividend_yield) * time_to_maturity));

	// standard error placeholders
	double sum_CT = 0.0;
	double sum_CT2 = 0.0;

	// run simulations
	for (int i = 0; i < M; i++) {
		double lnSt = lnS;
		for (int j = 0; j < N; j++) {
			lnSt = lnSt + nudt + volsdt * generate(primes, perms, currDim, N, globalSample);
			currDim++;
		}
		double ST = exp(lnSt);
		double CT = max(0.0, ST - strike_price);
		sum_CT += CT;
		sum_CT2 += pow(CT, 2);
	}

	double C0 = exp(-risk_free_interest_rate * time_to_maturity) * sum_CT / M;
	if (stderrorp) {
		double sigma = sqrt((sum_CT2 - (pow(sum_CT, 2) / M)) * exp(-2 * risk_free_interest_rate * time_to_maturity) / (M - 1));
		double SE = sigma / sqrt(M);
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
	// dimensions of sampling
	int currDim = 0, globalSample = 0;

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

	double dt = time_to_maturity / N;
	double nudt = (risk_free_interest_rate - 0.5 * pow(volatility, 2)) * dt;
	double volsdt = volatility * sqrt(dt);
	double lnS = log(underlying_price * exp((risk_free_interest_rate - dividend_yield) * time_to_maturity));

	// standard error placeholders
	double sum_CT = 0.0;
	double sum_CT2 = 0.0;

	// run simulations
	for (int i = 0; i < M; i++) {
		double lnSt = lnS;
		for (int j = 0; j < N; j++) {
			lnSt = lnSt + nudt + volsdt * generate(primes, perms, currDim, N, globalSample);
			currDim++;
		}
		double ST = exp(lnSt);
		double CT = max(0.0, strike_price-ST);
		sum_CT += CT;
		sum_CT2 += pow(CT, 2);
	}

	double C0 = exp(-risk_free_interest_rate * time_to_maturity) * sum_CT / M;
	if (stderrorp) {
		double sigma = sqrt((sum_CT2 - (pow(sum_CT, 2) / M)) * exp(-2 * risk_free_interest_rate * time_to_maturity) / (M - 1));
		double SE = sigma / sqrt(M);
		*stderrorp = SE;
	}
	return C0;
}


int main() {
	double S = 31.15;
	double K = 22.75;
	double vol = 0.5;
	double r = 0.05;
	double T = 3.5;
	int N = 10;
	int M = 10000;

	double* errorfactorCall = (double*) malloc(sizeof(double));
	double resultCall = callOptionValue(S, K, r, T, vol, 0.0, N, M, errorfactorCall);
	printf("Call option value: %f \n", resultCall);
	printf("Standard error: %f \n", *errorfactorCall);

	double* errorfactorPut = (double*)malloc(sizeof(double));
	double resultPut = putOptionValue(S, K, r, T, vol, 0.0, N, M, errorfactorPut);
	printf("Put option value: %f \n", resultPut);
	printf("Standard error: %f \n", *errorfactorPut);

	resultCall = callOptionValueHalton(S, K, r, T, vol, 0.0, N, M, errorfactorCall);
	printf("Call option value: %f \n", resultCall);
	printf("Standard error: %f \n", *errorfactorCall);


	resultCall = putOptionValueHalton(S, K, r, T, vol, 0.0, N, M, errorfactorPut);
	printf("Put option value: %f \n", resultCall);
	printf("Standard error: %f \n", *errorfactorPut);

	free(errorfactorCall);
	free(errorfactorPut);
	
}