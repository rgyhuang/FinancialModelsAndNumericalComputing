#include <iostream>
#include <cmath>
#include <random>
using namespace std;

/* 
* @brief Simulates the cumulative density function of a standard normal distribution
* @param x: input variable
* @returns value of cdf at x
*/
double normalCDF(double x) 
{
	return 0.5 * (1 + erf(x / sqrt(2)));
}

/*
* @brief Base implementation of Black-Scholes for European Call Otions Pricing
* @param underlying_price: price of stock/underlying asset
* @param strike_price: strike price of contract
* @param risk_free_interest_rate: risk-free interest rate, assumed to be constant
* @param time_to_maturity: time until expiration date of contract (years)
* @param volatility: volatility of underlying asset, assumed to be constant as well
* @returns call option price based on Black-Scholes model
*/

double callOptionValue(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility)
{

	double sigma_t = volatility * sqrt(time_to_maturity);
	double d_1 = (log(underlying_price / strike_price) + (time_to_maturity * (risk_free_interest_rate + (volatility / 2)))) / sigma_t;
	double d_2 = d_1 - sigma_t;

	return (normalCDF(d_1) * underlying_price) - (normalCDF(d_2) * strike_price * exp(-risk_free_interest_rate * time_to_maturity));

}

/*
* @brief Base implementation of Black-Scholes for European Put Options Pricing
* @param underlying_price: price of stock/underlying asset
* @param strike_price: strike price of contract
* @param risk_free_interest_rate: risk-free interest rate, assumed to be constant
* @param time_to_maturity: time until expiration date of contract (years)
* @param volatility: volatility of underlying asset, assumed to be constant as well
* @returns put option price based on Black-Scholes model
*/
double putOptionValue(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility)
{
	double sigma_t = volatility * sqrt(time_to_maturity);
	double d_1 = (log(underlying_price / strike_price) + (time_to_maturity * (risk_free_interest_rate + (volatility / 2)))) / sigma_t;
	double d_2 = d_1 - sigma_t;

	return (normalCDF(-d_2) * strike_price * exp(-risk_free_interest_rate * time_to_maturity)) - (normalCDF(-d_1) * underlying_price);
}

///*
//* @brief Black-Scholes for instruments paying continuous yield dividends
//* @param underlying_price: price of stock/underlying asset
//* @param strike_price: strike price of contract
//* @param risk_free_interest_rate: risk-free interest rate, assumed to be constant
//* @param time_to_maturity: time until expiration date of contract (assumed contract can only be exercised then)
//* @param volatility: volatility of underlying asset, assumed to be constant as well
//* @returns call option price based on Black-Scholes model
//*/
//
//double callOptionValue(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility)
//{
//
//	double sigma_t = volatility * sqrt(time_to_maturity);
//	double d_1 = (log(underlying_price / strike_price) + (time_to_maturity * (risk_free_interest_rate + (volatility / 2)))) / sigma_t;
//	double d_2 = d_1 - sigma_t;
//
//	return (normalCDF(d_1) * underlying_price) - (normalCDF(d_2) * strike_price * exp(-risk_free_interest_rate * time_to_maturity));
//
//}

// greeks





int main() {
	// testing with parameters
	double underlying_price = 138.38;
	double strike_price = 166.05;
	double risk_free_interest_rate = 0.0212;
	double time_to_maturity = 101;
	double volatility = 0.3364;
	
	// values
	double call = callOptionValue(underlying_price, strike_price, risk_free_interest_rate, time_to_maturity, volatility);
	double put = putOptionValue(underlying_price, strike_price, risk_free_interest_rate, time_to_maturity, volatility);
	cout << "Call: " << call << endl;
	cout << "Put: " << put << endl;
	//verify put-call parity
	cout << "Verifying put-call parity: " << endl;
	cout << "Call price: " << call << endl;
	cout << "Put: " << put + underlying_price - strike_price * exp(-risk_free_interest_rate * time_to_maturity) << endl;

}