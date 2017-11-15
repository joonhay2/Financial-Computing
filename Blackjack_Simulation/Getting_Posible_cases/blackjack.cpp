#define _SECURE_SCL 0

#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <time.h>
#include <sstream>
#include <chrono>
#include <fstream>


using namespace std;

typedef std::chrono::high_resolution_clock Clock;

int card_deq[52];
int trial;
int result_data[13][13][6];
int sub_total[13][13];
int table[10];
int cal_num[10];
double player_prob_data[13][13][6];
double dealer_prob_data[13][6];

void making_deq() {

	for (int i = 0; i < 52; i++) {
		card_deq[i] = (i / 4) + 1;
	}
	random_shuffle(&card_deq[0], &card_deq[51]);
}

void printing_deq() { // dont use
	for (int i = 0; i < 52; i++) {
		cout << i + 1 << endl;
		cout << "	[ " << card_deq[i] << " ]" << endl;

	}
}

string decide_card_name(int i) {
	string card_name;
	if (i == 1) card_name = "A";
	else if (i == 11) card_name = "J";
	else if (i == 12) card_name = "Q";
	else if (i == 13) card_name = "K";
	else {
		stringstream sstr;
		sstr << i;
		card_name = sstr.str();
	}

	return card_name;
}

int cal_total_num(int k) {
	int card_trial_sum = 0;
	for (int j = 0; j < 13; j++) {
		card_trial_sum += sub_total[k - 1][j];
	}
	return card_trial_sum;
}

int sub_total_num(int a, int b) {
	int sub_trial_sum = 0;
	for (int j = 0; j < 6; j++) {
		sub_trial_sum += result_data[a - 1][b - 1][j];
	}
	sub_total[a - 1][b - 1] = sub_trial_sum;
	return sub_trial_sum;
}

int test_num(int val) {
	if (val >= 10)
	{
		return 10;
	}
	else if (val == 1) return 11;
	else return val;
}
void clear_table() {
	for (int i = 0; i < 10; i++) {
		table[i] = 0;
	}
}

int summation(int i) {
	int sum = 0;

	for (int m = 0; m < i; m++) {
		sum += table[m];
	}

	return sum;
}

int hitting(int i) {

	table[i - 1] = test_num(card_deq[i - 1]);

	int temp_sum = summation(i);

	if (temp_sum > 21) {
		for (int j = 0; j < i; j++) {
			if (table[j] == 11) {
				table[j] = 1;
				return hitting(i);
			}
			else {
				return temp_sum;
			}
		}
	}
	else if (temp_sum < 17) {
		return hitting(i + 1);
	}
	else {
		return temp_sum;
	}

}

int dealer_sum(int a, int b) {

	int sum = 0;
	for (int i = 0; i < 13; i++) {
		sum += result_data[a - 1][i][b];
	}

	return sum;
}

void process_result() {

	for (int i = 1; i <= 13; i++) {

		for (int j = 1; j <= 13; j++) {

			int sub_card_sum = sub_total_num(i, j);

			player_prob_data[i - 1][j - 1][0] = (double)result_data[i - 1][j - 1][0] / (double)sub_card_sum * 100;
			player_prob_data[i - 1][j - 1][1] = (double)result_data[i - 1][j - 1][1] / (double)sub_card_sum * 100;
			player_prob_data[i - 1][j - 1][2] = (double)result_data[i - 1][j - 1][2] / (double)sub_card_sum * 100;
			player_prob_data[i - 1][j - 1][3] = (double)result_data[i - 1][j - 1][3] / (double)sub_card_sum * 100;
			player_prob_data[i - 1][j - 1][4] = (double)result_data[i - 1][j - 1][4] / (double)sub_card_sum * 100;
			player_prob_data[i - 1][j - 1][5] = (double)result_data[i - 1][j - 1][5] / (double)sub_card_sum * 100;
		}

		int card_sum = cal_total_num(i);

		dealer_prob_data[i - 1][0] = (double)dealer_sum(i, 0) / (double)card_sum * 100;
		dealer_prob_data[i - 1][1] = (double)dealer_sum(i, 1) / (double)card_sum * 100;
		dealer_prob_data[i - 1][2] = (double)dealer_sum(i, 2) / (double)card_sum * 100;
		dealer_prob_data[i - 1][3] = (double)dealer_sum(i, 3) / (double)card_sum * 100;
		dealer_prob_data[i - 1][4] = (double)dealer_sum(i, 4) / (double)card_sum * 100;
		dealer_prob_data[i - 1][5] = (double)dealer_sum(i, 5) / (double)card_sum * 100;

	}
}

void writing_dealer_data() {
	ofstream os;
	os.open("dealer_prob");

	for (int i = 1; i <= 13; i++) {

		os << i << ends;
		os << dealer_prob_data[i - 1][0] << ends;
		os << dealer_prob_data[i - 1][1] << ends;
		os << dealer_prob_data[i - 1][2] << ends;
		os << dealer_prob_data[i - 1][3] << ends;
		os << dealer_prob_data[i - 1][4] << ends;
		os << dealer_prob_data[i - 1][5] << ends;
	
	}
}


void writing_player_data() {
	ofstream os;
	os.open("player_prob");


	for (int i = 1; i <= 13; i++) {
		for (int j = 1; j <= 13; j++) {
			os << i << ends;
			os << j << ends;
			os << player_prob_data[i - 1][j - 1][0] << ends;
			os << player_prob_data[i - 1][j - 1][1] << ends;
			os << player_prob_data[i - 1][j - 1][2] << ends;
			os << player_prob_data[i - 1][j - 1][3] << ends;
			os << player_prob_data[i - 1][j - 1][4] << ends;
			os << player_prob_data[i - 1][j - 1][5] << ends;
		}
	}
}

void print_data() {
	cout << "-----------------------------------------" << endl;
	cout << "Total trial is " << trial << endl;
	cout << "-----------------------------------------" << endl;

	string first, second;

	cout << " For the dealer side, " << endl << endl;

	for (int i = 1; i <= 13; i++) {

		first = decide_card_name(i);

		cout << "  CARD :  " << first << endl;
		cout << "		 17 : " << dealer_prob_data[i - 1][0] << " %" << endl;
		cout << "		 18 : " << dealer_prob_data[i - 1][1] << " %" << endl;
		cout << "		 19 : " << dealer_prob_data[i - 1][2] << " %" << endl;
		cout << "		 20 : " << dealer_prob_data[i - 1][3] << " %" << endl;
		cout << "		 BJ : " << dealer_prob_data[i - 1][4] << " %" << endl;
		cout << "		 BUST : " << dealer_prob_data[i - 1][5] << " %" << endl;
		cout << "-----------------------------------------" << endl;
	}

	for (int i = 1; i <= 13; i++) {
		for (int j = 1; j <= 13; j++) {

			first = decide_card_name(i);
			second = decide_card_name(j);

			cout << "  CARD :  " << first << " & " << second << endl;
			cout << "		 17 : " << player_prob_data[i - 1][j - 1][0] << " %" << endl;
			cout << "		 18 : " << player_prob_data[i - 1][j - 1][1] << " %" << endl;
			cout << "		 19 : " << player_prob_data[i - 1][j - 1][2] << " %" << endl;
			cout << "		 20 : " << player_prob_data[i - 1][j - 1][3] << " %" << endl;
			cout << "		 BJ : " << player_prob_data[i - 1][j - 1][4] << " %" << endl;
			cout << "		 BUST : " << player_prob_data[i - 1][j - 1][5] << " %" << endl;
			cout << "-----------------------------------------" << endl;
		}
	}

}


void get_data(){
	ifstream input_file("dealer_prob");

}


void dealer() {


	making_deq();
	clear_table();
	int sum;

	table[0] = test_num(card_deq[0]);

	sum = hitting(2);


	if (sum > 21) result_data[card_deq[0] - 1][card_deq[1] - 1][5]++;
	else if (sum == 21) result_data[card_deq[0] - 1][card_deq[1] - 1][4]++;
	else if (sum == 20) result_data[card_deq[0] - 1][card_deq[1] - 1][3]++;
	else if (sum == 19) result_data[card_deq[0] - 1][card_deq[1] - 1][2]++;
	else if (sum == 18) result_data[card_deq[0] - 1][card_deq[1] - 1][1]++;
	else if (sum == 17) result_data[card_deq[0] - 1][card_deq[1] - 1][0]++;

}

int main() {

	cout << "Number of Trial? >> ";
	cin >> trial;
	srand(time(NULL));

	auto time_1 = Clock::now();

	for (int j = 0; j<13; j++)
		for (int k = 0; k<13; k++)
			for (int i = 0; i < 6; i++)
				result_data[j][k][i] = 0;

	// result_data[0~3] = seventeen ~ twenty
	// result data[4] = bj
	// result data[5] = bust


	for (int i = 0; i < trial; i++) {
		dealer();
	}

	process_result();

	auto time_2 = Clock::now();

	writing_dealer_data();
	writing_player_data();

	print_data();


	double time_diff_call = chrono::duration_cast<chrono::nanoseconds>(time_2 - time_1).count() / pow(10, 9);
	cout << " Time = " << time_diff_call << " sec" << endl;


}