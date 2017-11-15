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
double winning_prob[13][13][13];

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


void writing_player_data() {
	ofstream os;
	os.open("player_prob");


	for (int i = 1; i <= 13; i++) {
		for (int j = 1; j <= 13; j++) {
			os << i << j << endl << ends;
			os << player_prob_data[i - 1][j - 1][0] << endl << ends;
			os << player_prob_data[i - 1][j - 1][1] << endl << ends;
			os << player_prob_data[i - 1][j - 1][2] << endl << ends;
			os << player_prob_data[i - 1][j - 1][3] << endl << ends;
			os << player_prob_data[i - 1][j - 1][4] << endl << ends;
			os << player_prob_data[i - 1][j - 1][5] << endl << ends;
		}
	}
}

double winning_percentage_calculation(int pl_first, int pl_second, int dealer_card) {

	double winning_prob = 0;
	
	/*
	
	17 - bust
	18 - 17 or bust
	19 - 17, 18 or bust
	20 - 17, 18, 19 or bust
	21 
	
	*/
	
	winning_prob =
		(player_prob_data[pl_first - 1][pl_second - 1][0] * dealer_prob_data[dealer_card - 1][5] +
		player_prob_data[pl_first - 1][pl_second - 1][1] * (dealer_prob_data[dealer_card - 1][5] + dealer_prob_data[dealer_card - 1][0]) +
		player_prob_data[pl_first - 1][pl_second - 1][2] * (dealer_prob_data[dealer_card - 1][5] + dealer_prob_data[dealer_card - 1][0] + dealer_prob_data[dealer_card - 1][1]) +
		player_prob_data[pl_first - 1][pl_second - 1][3] * (dealer_prob_data[dealer_card - 1][5] + dealer_prob_data[dealer_card - 1][0] + dealer_prob_data[dealer_card - 1][1] + dealer_prob_data[dealer_card - 1][2]) +
		player_prob_data[pl_first - 1][pl_second - 1][4]*100 ) / 10000;

	return winning_prob;

}

void calculation() {
	// i, j = 플레이어 1, 2 카드; 
	//k = dealer 카드 
	for (int i = 1; i <= 13; i++) {
		for (int j = 1; j <= 13; j++) {
			for (int k = 1; k <= 13; k++) {
				winning_prob[i - 1][j - 1][k - 1] = winning_percentage_calculation(i, j, k);
				if ((i == 2) && (j == 3) && (k == 9)) {
					cout << "-----------------------------------------------------" << endl;
					cout << " Player : " << decide_card_name(i) << " & " << decide_card_name(j) << endl;
					cout << " Dealer : " << decide_card_name(k) << "	 Winning %  = " << winning_prob[i - 1][j - 1][k - 1] * 100 << " % " << endl;
				}
			}
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


void get_data() {

	int card;
	int player_first_card, player_second_card;

	ifstream input_file("dealer_prob_15000");

	if (input_file.is_open()) {
		for (int i = 1; i <= 13; i++) {

			input_file >> card;
			input_file >> dealer_prob_data[card - 1][0];
			input_file >> dealer_prob_data[card - 1][1];
			input_file >> dealer_prob_data[card - 1][2];
			input_file >> dealer_prob_data[card - 1][3];
			input_file >> dealer_prob_data[card - 1][4];
			input_file >> dealer_prob_data[card - 1][5];

		}
	}

	ifstream input_file_2("player_prob_15000");

	if (input_file_2.is_open()) {
		for (int i = 1; i <= 13; i++) {
			for (int j = 1; j <= 13; j++) {
				input_file_2 >> player_first_card;
				input_file_2 >> player_second_card;

				input_file_2 >> player_prob_data[player_first_card - 1][player_second_card - 1][0];
				input_file_2 >> player_prob_data[player_first_card - 1][player_second_card - 1][1];
				input_file_2 >> player_prob_data[player_first_card - 1][player_second_card - 1][2];
				input_file_2 >> player_prob_data[player_first_card - 1][player_second_card - 1][3];
				input_file_2 >> player_prob_data[player_first_card - 1][player_second_card - 1][4];
				input_file_2 >> player_prob_data[player_first_card - 1][player_second_card - 1][5];

			}
		}
	}

}


int main() {

	//cout << "Number of Trial? >> ";
	//cin >> trial;
	srand(time(NULL));

	auto time_1 = Clock::now();

	for (int j = 0; j<13; j++)
		for (int k = 0; k<13; k++)
			for (int i = 0; i < 6; i++)
				result_data[j][k][i] = 0;
	get_data();
	calculation();

	// result_data[0~3] = seventeen ~ twenty
	// result data[4] = bj
	// result data[5] = bust

	/*
	for (int i = 0; i < trial; i++) {
	dealer();
	}
	*/

	//ocess_result();

	auto time_2 = Clock::now();
	/*
	writing_dealer_data();
	writing_player_data();
	*/
	
	//print_data();


	double time_diff_call = chrono::duration_cast<chrono::nanoseconds>(time_2 - time_1).count() / pow(10, 9);
	cout << " Time = " << time_diff_call << " sec" << endl;


}