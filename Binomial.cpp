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
* @param height: height of the Binomial Tree
* @returns call option price based on Binomial Model
*/
double callOptionValue(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility, int height)
{
	double deltaT = time_to_maturity / height;
	double up_factor = exp(volatility * sqrt(deltaT));
	double down_factor = 1.0 / up_factor;
	
	//initializing result price vector
	vector<double> price_vector(height + 1);

	//intialize option values at expiration date
	for (int i = 0; i < height + 1; i++) {
		price_vector[i] = max(0.0, underlying_price * pow(up_factor, i) * pow(down_factor, height-i) - strike_price);
	}

	double compounded_return = exp(risk_free_interest_rate * deltaT);
	double risk_neutral_up = (compounded_return - down_factor) / (up_factor - down_factor);
	double risk_neutral_down = 1.0 - risk_neutral_up;

	// evaluate option values at earlier nodes
	for (int j = height-1; j >= 0; j--) {
		for (int i = 0; i <= j; i++) {
			price_vector[i] = exp(-risk_free_interest_rate * deltaT) * (risk_neutral_up * price_vector[i + 1] + risk_neutral_down * price_vector[i]);

			// exercise value, take the max(binomial value, exercise) for American options
			double exercise = strike_price - underlying_price * pow(up_factor, 2 * i - j);
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
* @param height: height of the Binomial Tree
* @returns call option price based on Binomial Model
*/
double putOptionValue(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility, int height)
{
	double deltaT = time_to_maturity / height;
	double up_factor = exp(volatility * sqrt(deltaT));
	double down_factor = 1.0 / up_factor;

	//initializing result price vector
	vector<double> price_vector(height + 1);

	//intialize option values at expiration date
	for (int i = 0; i < height + 1; i++) {
		price_vector[i] = max(0.0, strike_price - underlying_price * pow(up_factor, i) * pow(down_factor, height - i));
	}

	double compounded_return = exp(risk_free_interest_rate * deltaT);
	double risk_neutral_up = (compounded_return - down_factor) / (up_factor - down_factor);
	double risk_neutral_down = 1.0 - risk_neutral_up;

	// evaluate option values at earlier nodes
	for (int j = height - 1; j >= 0; j--) {
		for (int i = 0; i <= j; i++) {
			price_vector[i] = exp(-risk_free_interest_rate * deltaT) * (risk_neutral_up * price_vector[i + 1] + risk_neutral_down * price_vector[i]);

			// exercise value, take the max(binomial value, exercise) for American options
			double exercise = strike_price - underlying_price * pow(up_factor, 2 * i - j);
			if (price_vector[i] < exercise) price_vector[i] = exercise;

		}
	}
	return price_vector[0];

}


int main() {
	int N = 15000;
	double S0 = 100.0;
	double K = 100.0;
	double T = 1.0;
	double r = 0.1;
	double sigma = 0.2;
	cout << callOptionValue(S0, K, r, T, sigma, N) << endl;
	cout << putOptionValue(S0, K, r, T, sigma, N) << endl;

}