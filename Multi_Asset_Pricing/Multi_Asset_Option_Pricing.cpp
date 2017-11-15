/*
IE525 Numerical method in Finance
Pricing Multi-asset option

created by Joonha Yoon, Sid
*/


#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <string.h>
#include <algorithm>
#include <iomanip>
#include "newmat.h"
#include "newmatio.h"
#include <chrono>
#include <random>
#include <time.h>
typedef std::chrono::high_resolution_clock Clock; // to use clock!
												  //unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
												  //std::default_random_engine generator(seed);
using namespace std;
std::default_random_engine generator;

vector <string> Ticker_name;
vector <double> Volatility;
vector <double> Initial_Price;
vector <double> Nu, up_factor, down_factor;
vector <double> temp;
double Risk_free_rate;
double Strike_price;
double Time_to_maturity;
double No_of_steps;
//double global_memo_total;
//double global_binary_total;
double *memo;

int No_of_Assets;
int group_no;
Matrix Corr;
Matrix factor_index;
Matrix sign_assets;
Matrix Probability;
double delta_T, discount;

int counter = 0;
// Importing data --------------------------------


void Print_data() {

	cout << "Ticker	S_0	Rf	Sigma		T	K" << endl;
	for (int i = 0; i < (int)Ticker_name.size(); i++) {
		cout << Ticker_name[i] << "	" << Initial_Price[i] << "	" << Risk_free_rate << "	" << Volatility[i] << "	" << Time_to_maturity << "	" << Strike_price << endl;
	}
	cout << endl << "Correlation Coefficient Matrix :" << endl;
	cout <<endl<< setw(10) << setprecision(5) << Corr;


}

double max(double a, double b) {
	return (b < a) ? a : b;
}

double even_odd_check(int k) {
	if (k % 2 == 0) return 1;
	else return -1;
}

double summation(vector <double> k) {
	double temp = 0;
	for (int i = 0; i < (int)k.size(); i++)
		temp += k[i];
	return temp;
}

double max_new_general(vector <double> gen_asset) {
	auto biggest = max_element(begin(gen_asset), end(gen_asset));
	return *biggest;
}


void make_memo() { // NEVER try memoization if N > 4

	int memo_size = (pow(2, (No_of_steps + 1)*No_of_Assets) - 1) / (pow(2, No_of_Assets) - 1);
	memo = new double [memo_size]; // makeing memo function. using dynamic array.
	for (int i = 0; i < memo_size; i++) // put -1 to all.
			memo[i] = -1;

}

int assign_memo_group(int k) {
	if (k == 0) return 0;
	group_no = 0;
	for (int i = 1; i <= k - 1; i++) {
		group_no += pow(2, No_of_Assets*(i - 1));
	}
	return group_no;
}
/*
double memoization_option_pricing(int k, Matrix index) {
	if (k >= No_of_steps) {
		for (int i = 0; i < No_of_Assets; i++)
			temp[i] = Initial_Price[i] * pow(up_factor[i], index(1, i + 1));
		return max(max_new_general(temp) - Strike_price, 0.0);
	}
	else

	{
		double global_memo_total = 0;
		for (int i = 0; i < pow(2, No_of_Assets); i++)
			global_memo_total += discount*(Probability(i + 1, 1)*memoization_option_pricing(k + 1, (sign_assets.Row(i + 1) + index)));

		

		return global_memo_total;
	}
}


*/

void Get_Parameter_data() {
	ifstream file_1("C:/Users/com/Documents/info.csv");
	string value;
	bool read = true;
	bool pass_column_names = true;
	while (file_1.good())
	{
		if (pass_column_names) {
			for (int j = 0; j < 6; j++)
				getline(file_1, value, ',');
			pass_column_names = false;
		}
		getline(file_1, value, ',');
		getline(file_1, value, ',');
		if (string(value) == "") break; // if Ticker name is empty, stop reading
		Ticker_name.push_back(string(value));
		getline(file_1, value, ',');
		Initial_Price.push_back(atof(string(value).c_str()));
		getline(file_1, value, ',');
		Volatility.push_back(atof(string(value).c_str()));
		getline(file_1, value, ',');
		if (read) 	Risk_free_rate = atof(string(value).c_str());
		getline(file_1, value, ',');
		if (read) 	Time_to_maturity = atof(string(value).c_str());
		getline(file_1, value);
		if (read) {
			Strike_price = atof(string(value).c_str());
			read = false;
		}
	}

	No_of_Assets = Ticker_name.size();
	Matrix A(No_of_Assets, No_of_Assets);

	ifstream file_2("C:/Users/com/Documents/corr.csv");
	if (file_2.good())
	{
		for (int i = 1; i <= No_of_Assets + 1; i++)
				for (int j = 1; j <= No_of_Assets+1; j++) {
						if (j == No_of_Assets + 1) getline(file_2, value);
						else getline(file_2, value, ',');
						if ((i >= 2)&&(j>=2)) {
						A(i - 1, j - 1) = atof(string(value).c_str());
				}
			}
		Corr = A;
	}
	cout << " Data Read Complete " << endl;
}

void Initializing() {

	

	////////for testing////////////////

	No_of_steps = 4;
	Strike_price = 100;
	
	/////////////////////////////////////

	int No_of_prob = (int)pow(2, No_of_Assets);
	int No_of_corr = No_of_Assets*(No_of_Assets - 1) / 2;
	
	Matrix A(No_of_prob, No_of_Assets);
	Matrix sign_rho(No_of_prob, No_of_corr);
	Matrix Nu_devided_by_sigma(No_of_Assets, 1);
	Matrix Corr_column_vector(No_of_corr, 1);
	Matrix B(No_of_prob, 1);
	Matrix C(1, No_of_Assets);
	sign_assets = A;
	Probability = B;
	factor_index = C;

	delta_T = Time_to_maturity / No_of_steps;
	discount = exp(-Risk_free_rate*delta_T);


	for (int i = 0; i < No_of_Assets; i++) {
		Nu.push_back(Risk_free_rate - (0.5)*pow(Volatility[i], 2));
		up_factor.push_back(exp(Volatility[i] * sqrt(delta_T)));
		down_factor.push_back(pow(up_factor[i],-1));
		factor_index(1, i + 1) = 0;
		temp.push_back(0);
	}
	

	for (int i = 1; i <= No_of_Assets; i++) 
		Nu_devided_by_sigma(i, 1) = Nu[i - 1] / Volatility[i - 1];

	int count = 1;
	for (int j = 1; j <= No_of_Assets - 1; j++)
			for (int z = j + 1; z <= No_of_Assets; z++) {
				Corr_column_vector(count, 1) = Corr(j, z);
				count++;
			}

	int temp;
	for (int i = 1; i <= No_of_Assets; i++) {
		for (int j = 1; j <= No_of_prob; j++) {
			temp = (j-1) / ((int)pow(2, No_of_Assets - i));
			sign_assets(j, i) = even_odd_check(temp);
		}
	}

	
	for (int i = 1; i <= No_of_prob; i++) {
		 count = 1;
		for (int j = 1; j <= No_of_Assets-1; j++)
			for(int z = j+1; z<=No_of_Assets; z++){
				sign_rho(i, count) = sign_assets(i, j) * sign_assets(i, z);
				count++;
		}
	}

	Probability = (1 / (pow(2, No_of_Assets)) * (1 + (sign_assets*Nu_devided_by_sigma) *sqrt(delta_T) + sign_rho * Corr_column_vector));

	cout << endl << "A :" << endl;
	cout << endl << setw(10) << setprecision(2) << sign_assets;

	cout << endl << "risk_neutral_probs :" << endl;
	cout << endl << setw(10) << setprecision(7) << Probability;

}




double binomial_option_pricing(int k, Matrix index){	// memory problem !
	if (k >= No_of_steps) {	
		for (int i = 0; i < No_of_Assets; i++)
			temp[i] = Initial_Price[i] * pow(up_factor[i], index(1,i+1));
		return max(max_new_general(temp) - Strike_price, 0.0);
	}
	else 
	{ 
		double global_binary_total = 0;
		for (int i = 0; i < pow(2, No_of_Assets); i++) 
			global_binary_total += discount*(Probability(i + 1, 1)*binomial_option_pricing(k + 1, (sign_assets.Row(i + 1) + index)));
		return global_binary_total;
	}
}

//------------------------------------------------



int main() {

	cout << "---------------------------------------" << endl;

	Get_Parameter_data();
	
	Initializing();
	cout << "---------------------------------------" << endl;

	//Print_data();

	cout << "Option price using recursion : " << binomial_option_pricing(0,factor_index) << endl;
	cout << counter << endl;


	return 0;

}


