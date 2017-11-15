
// Pricing an European Option using Simulation
// Written by Prof. Sreenivas

#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <fstream>
#include <cstdlib>
#include "normdist.h" 
#include <random>
#include <time.h>
#define UNIT_STEP(x) ((x)>0?(1):(0))
using namespace std;
typedef std::chrono::high_resolution_clock Clock; // to use clock!
												  //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
												  //std::default_random_engine generator(seed);
std::default_random_engine generator;


double risk_free_rate, strike_price, initial_stock_price, expiration_time, volatility, barrier_price, no_of_barrier;

double put_option_monte_price = 0.0;
double call_option_monte_price = 0.0;
double brownian_call_price = 0.0;
double brownian_put_price = 0.0;

int no_of_trials;


double N(const double& z) {
	if (z > 6.0) { return 1.0; }; // this guards against overflow 
	if (z < -6.0) { return 0.0; };
	double b1 = 0.31938153;
	double b2 = -0.356563782;
	double b3 = 1.781477937;
	double b4 = -1.821255978;
	double b5 = 1.330274429;
	double p = 0.2316419;
	double c2 = 0.3989423;
	double a = fabs(z);
	double t = 1.0 / (1.0 + a*p);
	double b = c2*exp((-z)*(z / 2.0));
	double n = ((((b5*t + b4)*t + b3)*t + b2)*t + b1)*t;
	n = 1.0 - b*n;
	if (z < 0.0) n = 1.0 - n;
	return n;
}

double max(double a, double b) {
	return (b < a) ? a : b;
}

double get_uniform()
{
	std::uniform_real_distribution <double> distribution(0.0, 1.0);
	double number = distribution(generator);
	return (number);
}


double sigma(double i) {
	double ti = i * expiration_time * (1 / no_of_barrier);
	double sigma = sqrt(((ti - 0)*(expiration_time - ti) / (expiration_time - 0)));
	return sigma;
}

double mu(double a, double b, double i) {
	double ti = i * expiration_time * (1 / no_of_barrier);
	double mu = a + (ti - 0)*(b - a) / (expiration_time - 0);
	return mu;
}

double brownian_bridge(double a, double b) {
	double prob_d = 1;

	for (double i = 1; i < no_of_barrier; i++) {
		prob_d = prob_d * (1 - N((barrier_price - mu(a, b, i)) / sigma(i)));
	}

	return prob_d;
}

double barrier(double S) {
	if (S < barrier_price)
		return 0;
	else
		return S;
}

bool hit_barrier(double S, bool S_validity) {
	if (S_validity != true) {
		if (S < barrier_price)
			return true;
		else
			return false;
	}
	else
		return S_validity;
}

double proper_price(double S, bool option) {
	if (S == 0)
		return 0;
	else {
		if (option == true)
			return max(0.0, S - strike_price);
		else
			return max(0.0, strike_price - S);
	}
}

void run() {
	double delta_T = expiration_time / no_of_barrier;
	double delta_R = (risk_free_rate - 0.5*pow(volatility, 2))*delta_T;
	double delta_SD = volatility*sqrt(delta_T);
	bool S1_val = false;
	bool S2_val = false;
	bool S3_val = false;
	bool S4_val = false;

	for (int j = 0; j < no_of_trials; j++) {

		double S1 = initial_stock_price;
		double S2 = initial_stock_price;
		double S3 = initial_stock_price;
		double S4 = initial_stock_price;
		S1_val = false;
		S2_val = false;
		S3_val = false;
		S4_val = false;

		for (int i = 0; i < no_of_barrier; i++)
		{
			double x = get_uniform();
			double y = get_uniform();
			double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
			double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);

			S1 = S1 * exp(delta_R + delta_SD*a);
			S1_val = hit_barrier(S1, S1_val);

			S2 = S2 * exp(delta_R - delta_SD*a);
			S2_val = hit_barrier(S2, S2_val);

			S3 = S3 * exp(delta_R + delta_SD*b);
			S3_val = hit_barrier(S3, S3_val);

			S4 = S4 * exp(delta_R - delta_SD*b);
			S4_val = hit_barrier(S4, S4_val);

		}
		S1 = barrier(S1);
		S2 = barrier(S2);
		S3 = barrier(S3);
		S4 = barrier(S4);

		double S1_brownian_call = brownian_bridge(initial_stock_price, S1)* proper_price(S1, true);
		double S2_brownian_call = brownian_bridge(initial_stock_price, S2)* proper_price(S2, true);
		double S3_brownian_call = brownian_bridge(initial_stock_price, S3)* proper_price(S3, true);
		double S4_brownian_call = brownian_bridge(initial_stock_price, S4)* proper_price(S4, true);

		brownian_call_price += (S1_brownian_call + S2_brownian_call + S3_brownian_call + S4_brownian_call) / 4.0;

		double S1_brownian_put = brownian_bridge(initial_stock_price, S1)* proper_price(S1, false);
		double S2_brownian_put = brownian_bridge(initial_stock_price, S2)* proper_price(S2, false);
		double S3_brownian_put = brownian_bridge(initial_stock_price, S3)* proper_price(S3, false);
		double S4_brownian_put = brownian_bridge(initial_stock_price, S4)* proper_price(S4, false);

		brownian_put_price += (S1_brownian_put + S2_brownian_put + S3_brownian_put + S4_brownian_put) / 4.0;


		if (S1_val == true) S1 = 0.0;
		if (S2_val == true) S2 = 0.0;
		if (S3_val == true) S3 = 0.0;
		if (S4_val == true) S4 = 0.0;

		call_option_monte_price += (proper_price(S1, true) +
			proper_price(S2, true) +
			proper_price(S3, true) +
			proper_price(S4, true)) / 4.0;
		put_option_monte_price += (proper_price(S1, false) +
			proper_price(S2, false) +
			proper_price(S3, false) +
			proper_price(S4, false)) / 4.0;

	}
	call_option_monte_price = exp(-risk_free_rate*expiration_time)*(call_option_monte_price / ((double)no_of_trials));
	put_option_monte_price = exp(-risk_free_rate*expiration_time)*(put_option_monte_price / ((double)no_of_trials));
	brownian_call_price = exp(-risk_free_rate*expiration_time)*(brownian_call_price / ((double)no_of_trials));
	brownian_put_price = exp(-risk_free_rate*expiration_time)*(brownian_put_price / ((double)no_of_trials));

}


int main(int argc, char* argv[])
{
	sscanf_s(argv[1], "%lf", &expiration_time);
	sscanf_s(argv[2], "%lf", &risk_free_rate);
	sscanf_s(argv[3], "%lf", &volatility);
	sscanf_s(argv[4], "%lf", &initial_stock_price);
	sscanf_s(argv[5], "%lf", &strike_price);
	sscanf_s(argv[6], "%d", &no_of_trials);
	sscanf_s(argv[7], "%lf", &no_of_barrier);
	sscanf_s(argv[8], "%lf", &barrier_price);

	auto t_1 = Clock::now();

	run();

	auto t_2 = Clock::now();


	double time_diff_call = chrono::duration_cast<chrono::nanoseconds>(t_2 - t_1).count() / pow(10, 9);

	cout << "--------------------------------" << endl;
	cout << "European Down-and-Out Discrete Barrier Options Pricing via Monte Carlo Simulation" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "Barrier Price = " << barrier_price << endl;
	cout << "Number of Trials = " << no_of_trials << endl;
	cout << "Number of Barriers = " << no_of_barrier << endl;
	cout << "--------------------------------" << endl;

	cout << "The average Call Price by explicit simulation = " << call_option_monte_price << endl;
	cout << "The average Call Price with Brownian-Bridge correction on the final price  = " << brownian_call_price << endl;
	
	cout << "The average Put Price by explicit simulation " << put_option_monte_price << endl;
	cout << "The average Put Price with Brownian-Bridge correction on the final price  = " << brownian_put_price << endl;
	
	cout << "--------------------------------" << endl;
	cout << "Total  " << time_diff_call << " sec" << endl;

}
