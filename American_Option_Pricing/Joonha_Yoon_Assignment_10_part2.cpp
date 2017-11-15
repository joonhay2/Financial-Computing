// Black-Scholes European Option Pricing Code
// Adapted from Prof. Odegaard's Notes
// Written by Prof. Sreenivas for IE523: Financial Computing
/*
#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include "normdist.h"          // this defines the normal distribution from Odegaard's files
using namespace std;

float up_factor, uptick_prob, risk_free_rate, strike_price, downtick_prob, notick_prob;
float initial_stock_price, expiration_time, volatility, R;
float temp;
int no_of_divisions;
float **memo;

float max(float a, float b) {
	return (b < a) ? a : b;
}

/*
void make_memo() {
	memo = new float *[no_of_divisions];
	for (int i = 0; i <= no_of_divisions; i++)
		memo[i] = new float[(2*i+1)];

	cout <<endl<< "					 MAKE memo complete@@@@@@@@@@@@@@@@@@@@" << endl;
}


void clean_memo(float *memo_array[]) {			//왜 Array 크기를 0으로만들어도 Array가 생기는지 질문!
												// 왜 Array를 만들면, -1, -2, -3 에도 자리가 생기는지 질문.
	for (int i = 0; i <= no_of_divisions; i++)
		for (int j = 0; j <= 2 * i + 1; j++)  // put zeros to all.
			memo_array[i][j] = -1;

	cout <<endl<< "					 CLEAN memo complete@@@@@@@@@@@@@@@@@@@@" << endl;
}

void print_memo(float *memo_array[]) {
	for (int i = 0; i < no_of_divisions; i++){
		for (int j = 0; j < 2 * i + 1; j++) {  // put zeros to all.
			cout << memo_array[i][j] << " ";
		}
		cout << endl;				
	}
	//cout << " 메모 0,1" << memo[1][0] << endl;
}





void make_memo() {
	memo = new float *[2*no_of_divisions];
	for (int i = 0; i <= 2*no_of_divisions; i++)
		memo[i] = new float[no_of_divisions];
}

void clean_memo(float **memo_array) {
	for (int i = 0; i <= 2*no_of_divisions; i++)
		for (int j = 0; j <= no_of_divisions; j++)  // put zeros to all.
			memo_array[i][j] = -1;
}
/*
void print_memo(float **memo_array) {
	for (int i = 0; i <= 2 * no_of_divisions; i++){
		for (int j = 0; j <= no_of_divisions; j++) {  
			cout << memo_array[i][j] << " ";
		}
		cout << endl;
	}
	//cout << " 메모 0,1" << memo[1][0] << endl;
}
*/
/*
void print_memo(float *memo_array[]) {
	for (int i = 0; i < no_of_divisions; i++) {
		for (int j = 0; j <= no_of_divisions; j++) {
			cout << memo_array[i][j] << " ";
		}
		cout << endl;
	}
	cout << endl << endl << " apah " << memo[2 * no_of_divisions - 1][0];
	delete memo[1];
	cout << endl<<endl << " apah " << memo[2 * no_of_divisions - 1][0];
	for (int i = no_of_divisions+3; i <= 2 * no_of_divisions; i++) {
		for (int j = 0; j <= no_of_divisions; j++) {
			cout << memo_array[i][j] << " ";
		}
		cout << endl;
	}
	//cout << " 메모 0,1" << memo[1][0] << endl;
}


double option_price_put_black_scholes(const double& S,      // spot price
	const double& K,      // Strike (exercise) price,
	const double& r,      // interest rate
	const double& sigma,  // volatility
	const double& time) {
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r*time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
	double d2 = d1 - (sigma*time_sqrt);
	return K*exp(-r*time)*N(-d2) - S*N(-d1);
};

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
	const double& K,       // strike (exercise) price,
	const double& r,       // interest rate
	const double& sigma,   // volatility 
	const double& time) {  // time to maturity 
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r*time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
	double d2 = d1 - (sigma*time_sqrt);
	return S*N(d1) - K*exp(-r*time)*N(d2);
};

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
};


float european_call_option(int k, int i) {

	if (memo[i][k] != -1) { 
		
		
		return memo[i][k]; 
	}
	else if (k >= no_of_divisions) 
	{
		//cout << "the price is   " << max(0.0, (strike_price - initial_stock_price*pow(up_factor, ((float)(i - k))))) << endl;
	//	cout << "k is   "<< k <<"  i is   "<<i<<"     reach the end " << endl;
		memo[i][k] = max(0.0, (initial_stock_price*pow(up_factor, ((float)(i-no_of_divisions)))) - strike_price);
		return memo[i][k];
	}
	else 
	{
			memo[i][k] =((uptick_prob*european_call_option(k + 1, i + 1) + notick_prob*european_call_option(k + 1, i) + downtick_prob*european_call_option(k + 1, i-1)) / R);
			if ((i > no_of_divisions)&&(i<2*no_of_divisions-1)&&(k>1)) {
			//	delete[] memo[i + 2];
			}
			return memo[i][k];
		}
}
/*
float european_call_option(int k, int i) {

	//cout << endl << "logic start!!" << endl;

	if (memo[k][i] != -1) { 
		temp = memo[k][i];
		delete[] memo[k];
		return memo[k][i]; 
	}
	else if (k >= no_of_divisions)
	{
		//	cout << "k is   "<< k <<"  i is   "<<i<<"     reach the end " << endl;
		memo[k][i] = max(0.0, (initial_stock_price*pow(up_factor, ((float)(i - k)))) - strike_price);
		return memo[k][i];
	}
	else
	{
		memo[k][i] = ((uptick_prob*european_call_option(k + 1, i + 2) + notick_prob*european_call_option(k + 1, i + 1) + downtick_prob*european_call_option(k + 1, i )) / R);
		return memo[k][i];
	}

	//if
	//	cout << "the price of  k = "<<k<<"  i = "<<i<<"  is  " << ((uptick_prob*memo[k + 1][i + 2] + notick_prob*memo[k + 1][i + 1] + downtick_prob*memo[k + 1][i  ]) / R) << endl;
	//return ((uptick_prob*memo[k + 1][i + 2] + notick_prob*memo[k + 1][i + 1] + downtick_prob*memo[k + 1][i]) / R);

}
*/



/*
float european_call_option(int k, int i) {
	if (k == no_of_divisions)
		return max(0.0, (initial_stock_price*pow(up_factor, ((float)i))) - strike_price);
	else
		return ((uptick_prob*european_call_option(k + 1, i + 1) + (notick_prob)*european_call_option(k + 1, i)+
		(downtick_prob)*european_call_option(k + 1, i - 1)) / R);
}


float european_put_option(int k, int i) {
	if (k == no_of_divisions)
		return max(0.0, strike_price - (initial_stock_price*pow(up_factor, ((float)i))));
	else
		return ((uptick_prob*european_put_option(k + 1, i + 1) +
		(1 - uptick_prob)*european_put_option(k + 1, i - 1)) / R);
}

int main(int argc, char* argv[])
{

	sscanf(argv[1], "%f", &expiration_time);
	sscanf(argv[2], "%d", &no_of_divisions);
	sscanf(argv[3], "%f", &risk_free_rate);
	sscanf(argv[4], "%f", &volatility);
	sscanf(argv[5], "%f", &initial_stock_price);
	sscanf(argv[6], "%f", &strike_price);


	make_memo();
	clean_memo(memo);
//	print_memo(memo);


	up_factor = exp(volatility*sqrt(2 * expiration_time / ((float)no_of_divisions)));
	R = exp(risk_free_rate*expiration_time / ((float)no_of_divisions));
	uptick_prob = pow(((sqrt(R) - (1 / (sqrt(up_factor)))) / (sqrt(up_factor) - (1 / sqrt(up_factor)))), 2);
	downtick_prob = pow(((sqrt(up_factor) - sqrt(R)) / (sqrt(up_factor) - (1 / sqrt(up_factor)))), 2);
	notick_prob = 1 - (downtick_prob + uptick_prob);


	cout << "Recursive Binomial European Option Pricing" << endl;
	cout << "Expiration Time (Years) = " << expiration_time << endl;
	cout << "Number of Divisions = " << no_of_divisions << endl;
	cout << "Risk Free Interest Rate = " << risk_free_rate << endl;
	cout << "Volatility (%age of stock value) = " << volatility * 100 << endl;
	cout << "Initial Stock Price = " << initial_stock_price << endl;
	cout << "Strike Price = " << strike_price << endl;
	cout << "--------------------------------------" << endl;
	cout << "Up Factor = " << up_factor << endl;
	cout << "Uptick Probability = " << uptick_prob << endl;
	cout << "--------------------------------------" << endl;
	double call_price = european_call_option(0, no_of_divisions);
//	double call_price = european_call_option(0, 0);
	cout << "Trinomial Price of an European Call Option = " << call_price << endl;
	cout << "Call Price according to Black-Scholes = " <<
		option_price_call_black_scholes(initial_stock_price, strike_price, risk_free_rate,
			volatility, expiration_time) << endl;
	cout << "--------------------------------------" << endl;
	//double put_price = european_put_option(0, 0);
//	cout << "Binomial Price of an European Put Option = " << put_price << endl;
	cout << "Put Price according to Black-Scholes = " <<
		option_price_put_black_scholes(initial_stock_price, strike_price, risk_free_rate,
			volatility, expiration_time) << endl;
	cout << "--------------------------------------" << endl;
	cout << "Verifying Put-Call Parity: S+P-C = Kexp(-r*T)" << endl;
//	cout << initial_stock_price << " + " << put_price << " - " << call_price;
	cout << " = " << strike_price << "exp(-" << risk_free_rate << " * " << expiration_time << ")" << endl;
	//cout << initial_stock_price + put_price - call_price << " = " << strike_price*exp(-risk_free_rate*expiration_time) << endl;
	cout << "--------------------------------------" << endl;
}

*/