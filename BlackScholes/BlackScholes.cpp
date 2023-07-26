#include <iostream>
#include <cmath>
#include <numbers>
#include <random>
#include <iomanip>
using namespace std;

#define M_PI 3.14159265358979323846

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
* @brief Simulates the probability density function of a standard normal distribution
* @param x: input variable
* @returns value of pdf at x
*/
double normalPDF(double x)
{
	return (1 / sqrt(2 * M_PI)) * exp(-0.5 * x * x);
}


/*
* @brief Base implementation of Black-Scholes for European Call Otions Pricing
* @param underlying_price: price of stock/underlying asset
* @param strike_price: strike price of contract
* @param risk_free_interest_rate: risk-free interest rate, assumed to be constant
* @param time_to_maturity: years until expiration date of contract
* @param volatility: volatility of underlying asset, assumed to be constant as well
* @param dividend_yield: for instruments paying dividends
* @returns call option price based on Black-Scholes model
*/

double callOptionValue(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility, double dividend_yield)
{

	double sigma_t = volatility * sqrt(time_to_maturity);
	double F = underlying_price * exp((risk_free_interest_rate - dividend_yield) * time_to_maturity);
	double d_1 = (log(underlying_price / strike_price) + (time_to_maturity * (risk_free_interest_rate - dividend_yield + (volatility*volatility / 2)))) / sigma_t;
	double d_2 = d_1 - sigma_t;

	return exp(-risk_free_interest_rate*time_to_maturity)*((F*(normalCDF(d_1))) - (normalCDF(d_2) * strike_price));

}


/*
* @brief Base implementation of Black-Scholes for European Put Options Pricing
* @param underlying_price: price of stock/underlying asset
* @param strike_price: strike price of contract
* @param risk_free_interest_rate: risk-free interest rate, assumed to be constant
* @param time_to_maturity: years until expiration date of contract
* @param volatility: volatility of underlying asset, assumed to be constant as well
* @param dividend_yield: for instruments paying dividends
* @returns put option price based on Black-Scholes model
*/
double putOptionValue(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility, double dividend_yield)
{
	double sigma_t = volatility * sqrt(time_to_maturity);
	double F = underlying_price * exp((risk_free_interest_rate - dividend_yield) * time_to_maturity);
	double d_1 = (log(underlying_price / strike_price) + (time_to_maturity * (risk_free_interest_rate - dividend_yield + (volatility * volatility / 2)))) / sigma_t;
	double d_2 = d_1 - sigma_t;

	return exp(-risk_free_interest_rate * time_to_maturity) * ((normalCDF(-d_2) * strike_price) - (F * normalCDF(-d_1)));
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
//	                                                                                     //
//	 Options Greeks Implementations                                                      //
//	 @brief the following Greeks functions have the same parameters:                     //
//	 @param underlying_price: price of stock/underlying asset                            //
//	 @param strike_price: strike price of contract                                       //
//	 @param risk_free_interest_rate: risk-free interest rate, assumed to be constant     //
//	 @param time_to_maturity: years until expiration date of contract                    //
//	 @param volatility: volatility of underlying asset, assumed to be constant as well   //
//	 @param dividend_yield: for instruments paying dividends                             //
//	                                                                                     //
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Call options

// Delta call option greek
// Measures change in an option's price or premium (dV) from a change in the underlying asset (dS)
double callDelta(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility, double dividend_yield) {
	double sigma_t = volatility * sqrt(time_to_maturity);
	double d_1 = (log(underlying_price / strike_price) + (time_to_maturity * (risk_free_interest_rate - dividend_yield + (volatility * volatility / 2)))) / sigma_t;
	return exp(-dividend_yield*time_to_maturity)*normalCDF(d_1);
}

// Gamma call option greek
// Measures the delta's rate of change over time, as well as the rate of change in the underlying asset (d^2V/dS^2)
// Helps forecast price moves in the underlying asset
double callGamma(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility, double dividend_yield) {
	double sigma_t = volatility * sqrt(time_to_maturity);
	double d_1 = (log(underlying_price / strike_price) + (time_to_maturity * (risk_free_interest_rate - dividend_yield + (volatility * volatility / 2)))) / sigma_t;
	return exp(-dividend_yield * time_to_maturity) * (normalPDF(d_1) / (underlying_price * sigma_t));
}


// Vega call option greek
// Measures the risk of changes in implied volatility or the forward looking expected volatility of the underlying asset price (dV/d\sigma)
// Note: result is divided by 100 to get the option price change for one percentage point change in volality
double callVega(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility, double dividend_yield) {
	double sigma_t = volatility * sqrt(time_to_maturity);
	double d_1 = (log(underlying_price / strike_price) + (time_to_maturity * (risk_free_interest_rate - dividend_yield + (volatility * volatility / 2)))) / sigma_t;
	// price change per 1% change in volatility
	return (underlying_price * exp(-dividend_yield * time_to_maturity) * normalPDF(d_1) * sqrt(time_to_maturity)) / 100;
}

// Theta call option greek
// Measures time decay inthe value of an option or its premium (dV/dt)
// Note: extra parameter days_in_year (365 if calendar and about 252 for trading days)
double callTheta(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility, double dividend_yield, double days_in_year) {
	double sigma_t = volatility * sqrt(time_to_maturity);
	double _coeff = exp(-dividend_yield * time_to_maturity);
	double d_1 = (log(underlying_price / strike_price) + (time_to_maturity * (risk_free_interest_rate - dividend_yield + (volatility * volatility / 2)))) / sigma_t;
	double d_2 = d_1 - sigma_t;


	return (-(volatility * _coeff / (2 * sqrt(time_to_maturity)) * normalPDF(d_1))
			- (risk_free_interest_rate * strike_price * exp(-risk_free_interest_rate * time_to_maturity) * normalCDF(d_2))
			+ (dividend_yield * _coeff * normalCDF(d_1))) / days_in_year;
}

// Rho call option greek
// Measures change in option price with respect to risk free interest rate (dV/dr, per 1% change in interest rate) 
double callRho(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility, double dividend_yield) {
	double sigma_t = volatility * sqrt(time_to_maturity);
	double d_1 = (log(underlying_price / strike_price) + (time_to_maturity * (risk_free_interest_rate - dividend_yield + (volatility * volatility / 2)))) / sigma_t;
	double d_2 = d_1 - sigma_t;
	return (strike_price * time_to_maturity * exp(-risk_free_interest_rate * time_to_maturity) * normalCDF(d_2)) / 100;
}

// Put options

// Delta put option greek
// Measures change in an option's price or premium (dV) from a change in the underlying asset (dS)
double putDelta(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility, double dividend_yield) {
	double sigma_t = volatility * sqrt(time_to_maturity);
	double d_1 = (log(underlying_price / strike_price) + (time_to_maturity * (risk_free_interest_rate - dividend_yield + (volatility * volatility / 2)))) / sigma_t;
	return -exp(-dividend_yield * time_to_maturity) * normalCDF(-d_1);
}

// Gamma put option greek
// Measures the delta's rate of change over time, as well as the rate of change in the underlying asset (d^2V/dS^2)
// Helps forecast price moves in the underlying asset
double putGamma(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility, double dividend_yield) {
	double sigma_t = volatility * sqrt(time_to_maturity);
	double d_1 = (log(underlying_price / strike_price) + (time_to_maturity * (risk_free_interest_rate - dividend_yield + (volatility * volatility / 2)))) / sigma_t;
	return exp(-dividend_yield * time_to_maturity) * (normalPDF(d_1) / (underlying_price * sigma_t));
}


// Vega put option greek
// Measures the risk of changes in implied volatility or the forward looking expected volatility of the underlying asset price (dV/d\sigma)
// Note: result is divided by 100 to get the option price change for one percentage point change in volality
double putVega(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility, double dividend_yield) {
	double sigma_t = volatility * sqrt(time_to_maturity);
	double d_1 = (log(underlying_price / strike_price) + (time_to_maturity * (risk_free_interest_rate - dividend_yield + (volatility * volatility / 2)))) / sigma_t;
	// price change per 1% change in volatility
	return (underlying_price * exp(-dividend_yield * time_to_maturity) * normalPDF(d_1) * sqrt(time_to_maturity)) / 100;
}

// Theta put option greek
// Measures time decay inthe value of an option or its premium (dV/dt)
// Note: extra parameter days_in_year (365 if calendar and about 252 for trading days)
double putTheta(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility, double dividend_yield, double days_in_year) {
	double sigma_t = volatility * sqrt(time_to_maturity);
	double _coeff = exp(-dividend_yield * time_to_maturity);
	double d_1 = (log(underlying_price / strike_price) + (time_to_maturity * (risk_free_interest_rate - dividend_yield + (volatility * volatility / 2)))) / sigma_t;
	double d_2 = d_1 - sigma_t;


	return (-(volatility * _coeff / (2 * sqrt(time_to_maturity)) * normalPDF(d_1))
		+ (risk_free_interest_rate * strike_price * exp(-risk_free_interest_rate * time_to_maturity) * normalCDF(-d_2))
		- (dividend_yield * _coeff * normalCDF(-d_1))) / days_in_year;
}

// Rho put option greek
// Measures change in option price with respect to risk free interest rate (dV/dr, per 1% change in interest rate) 
double putRho(double underlying_price, double strike_price, double risk_free_interest_rate, double time_to_maturity, double volatility, double dividend_yield) {
	double sigma_t = volatility * sqrt(time_to_maturity);
	double d_1 = (log(underlying_price / strike_price) + (time_to_maturity * (risk_free_interest_rate - dividend_yield + (volatility * volatility / 2)))) / sigma_t;
	double d_2 = d_1 - sigma_t;
	return (strike_price * time_to_maturity * exp(-risk_free_interest_rate * time_to_maturity) * normalCDF(-d_2)) / -100;
}


int main() {
	// testing with parameters
	double underlying_price = 31.55;
	double strike_price = 22.75;
	double risk_free_interest_rate = 0.05;
	double time_to_maturity = 3.5;
	double volatility = 0.5;
	double dividend_yield = 0.15;

	// values
	double call = callOptionValue(underlying_price, strike_price, risk_free_interest_rate, time_to_maturity, volatility, dividend_yield);
	double put = putOptionValue(underlying_price, strike_price, risk_free_interest_rate, time_to_maturity, volatility, dividend_yield);
	cout << setprecision(3) <<	"Call: $" << call << endl;
	cout << setprecision(3) << "Put: $" << put << endl;
	
	// greeks 
	double delta_c = callDelta(underlying_price, strike_price, risk_free_interest_rate, time_to_maturity, volatility, dividend_yield);
	double gamma_c = callGamma(underlying_price, strike_price, risk_free_interest_rate, time_to_maturity, volatility, dividend_yield);
	double vega_c = callVega(underlying_price, strike_price, risk_free_interest_rate, time_to_maturity, volatility, dividend_yield);
	double theta_c = callTheta(underlying_price, strike_price, risk_free_interest_rate, time_to_maturity, volatility, dividend_yield, 365);
	double rho_c = callRho(underlying_price, strike_price, risk_free_interest_rate, time_to_maturity, volatility, dividend_yield);

	double delta_p = putDelta(underlying_price, strike_price, risk_free_interest_rate, time_to_maturity, volatility, dividend_yield);
	double gamma_p = putGamma(underlying_price, strike_price, risk_free_interest_rate, time_to_maturity, volatility, dividend_yield);
	double vega_p = putVega(underlying_price, strike_price, risk_free_interest_rate, time_to_maturity, volatility, dividend_yield);
	double theta_p = putTheta(underlying_price, strike_price, risk_free_interest_rate, time_to_maturity, volatility, dividend_yield, 365);
	double rho_p = putRho(underlying_price, strike_price, risk_free_interest_rate, time_to_maturity, volatility, dividend_yield);

	cout << "Greeks" << endl;
	cout << "Call: " << endl;
	printf("delta: %f || gamma: %f || vega: %f || theta: %f || rho: %f\n", delta_c, gamma_c, vega_c, theta_c, rho_c);
	cout << "Put: " << endl;
	printf("delta: %f || gamma: %f || vega: %f || theta: %f || rho: %f\n", delta_p, gamma_p, vega_p, theta_p, rho_p);


}
