#include <iostream>
#include <cmath>
#include <random>
using namespace std;

/*
* @brief Base implementation of Binomial Call Option Pricing Model for American Options
* @param underlying_price: price of stock/underlying asset
* @param strike_price: strike price of contract
* @param risk_free_interest_rate: risk-free interest rate
* @param time_to_maturity: time until expiration date of contract (assumed contract can only be exercised then)
* @param volatility: volatility of underlying asset
* @param dividend_yield: dividend yield for instruments that pay in dividends
* @param height: height of the Binomial Tree
* @returns call option price based on Binomial Model
*/
double callOptionValue(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility, double dividend_yield, int height)
{
	double deltaT = time_to_maturity / height;
	double up_factor = exp(volatility * sqrt(deltaT));
	double p0 = (up_factor * exp(-dividend_yield * deltaT) - exp(-risk_free_interest_rate * deltaT)) / ((up_factor * up_factor) - 1);
	double p1 = exp(-risk_free_interest_rate * deltaT) - p0;

	//initializing result price vector
	vector<double> price_vector(height + 1);

	//intialize option values at expiration date
	for (int i = 0; i < height + 1; i++) {
		price_vector[i] = max(0.0, (underlying_price*pow(up_factor, 2*i - height)) - strike_price);
	}

	// evaluate option values at earlier nodes
	for (int j = height - 1; j >= 0; j--) {
		for (int i = 0; i <= j; i++) {
			price_vector[i] = p0 * price_vector[i + 1] + p1 * price_vector[i];

			// exercise value, take the max(binomial value, exercise) for American options
			double exercise = underlying_price * pow(up_factor, (2 * i) - j)- strike_price;
			if (price_vector[i] < exercise) price_vector[i] = exercise;

		}
	}
	return price_vector[0];

}

/*
* @brief Base implementation of Binomial Put Option Pricing Model for American Options
* @param underlying_price: price of stock/underlying asset
* @param strike_price: strike price of contract
* @param risk_free_interest_rate: risk-free interest rate
* @param time_to_maturity: time until expiration date of contract (assumed contract can only be exercised then)
* @param volatility: volatility of underlying asset
* @param dividend_yield: dividend yield for instruments that pay in dividends
* @param height: height of the Binomial Tree
* @returns call option price based on Binomial Model
*/
double putOptionValue(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility, double dividend_yield, int height)
{
	double deltaT = time_to_maturity / height;
	double up_factor = exp(volatility * sqrt(deltaT));
	double p0 = (up_factor * exp(-dividend_yield * deltaT) - exp(-risk_free_interest_rate * deltaT)) / ((up_factor * up_factor) - 1);
	double p1 = exp(-risk_free_interest_rate * deltaT) - p0;

	//initializing result price vector
	vector<double> price_vector(height + 1);

	//intialize option values at expiration date
	for (int i = 0; i < height + 1; i++) {
		price_vector[i] = max(0.0, strike_price - (underlying_price * pow(up_factor, 2 * i - height)));
	}

	// evaluate option values at earlier nodes
	for (int j = height - 1; j >= 0; j--) {
		for (int i = 0; i <= j; i++) {
			price_vector[i] = p0 * price_vector[i + 1] + p1 * price_vector[i];

			// exercise value, take the max(binomial value, exercise) for American options
			double exercise = underlying_price * pow(up_factor, 2 * i - j) - strike_price;
			if (price_vector[i] < exercise) price_vector[i] = exercise;

		}
	}
	return price_vector[0];
}
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
//	                                                                                     //
//	 Optimizing Binomial Tree Pricing Model Using Parallel Computing                     //
//	 Implementation derived from algorithm proposed by Popuri et al. (link to paper)     //
//   Traditional Binomial Methods are based off of iterating over discounted option      //
//   payoffs ina recursive manner. Popuri et al. propose an optmize procedure that maps  // 
//   Binomial probabilities to Bernoulli paths                                           //
//	                                                                                     //
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
* @brief Parallelized Implementation of Binomial Call Option Model
* @param underlying_price: price of stock/underlying asset
* @param strike_price: strike price of contract
* @param risk_free_interest_rate: risk-free interest rate
* @param time_to_maturity: time until expiration date of contract (assumed contract can only be exercised then)
* @param volatility: volatility of underlying asset
* @param dividend_yield: dividend yield for instruments that pay in dividends
* @param height: height of the Binomial Tree
* @returns call option price based on Binomial Model
*/
double callOptionValue(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility, double dividend_yield, int height)
{

}

int main() {
	int N = 150;
	double S0 = 100.0;
	double K = 100.0;
	double T = 1.0;
	double r = 0.1;
	double sigma = 0.2;
	double q = 0.0;
	cout << callOptionValue(S0, K, r, T, sigma, q, N) << endl;
	cout << putOptionValue(S0, K, r, T, sigma, q, N) << endl;

}
