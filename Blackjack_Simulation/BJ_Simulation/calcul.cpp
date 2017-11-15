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
int dealer_table[10];
int player_table[10];
int cal_num[10];
int initial_money;
int betting;
int win_count = 0;
int loss_count = 0;
int money;
double player_prob_data[13][13][6];
double dealer_prob_data[13][6];
double winning_prob[13][13][13];

vector<int> game_result;



void making_deq() {

	for (int i = 0; i < 52; i++) {
		card_deq[i] = (i / 4) + 1;
	}
	random_shuffle(&card_deq[0], &card_deq[51]);
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
		dealer_table[i] = 0;
		player_table[i] = 0;
	}
}

int player_summation(int i) {
	int sum = 0;

	for (int m = 0; m < i; m++) {
		sum += player_table[m];
	}

	return sum;
}

int dealer_summation(int i) {
	int sum = 0;

	for (int m = 0; m < i; m++) {
		sum += dealer_table[m];
	}

	return sum;
}

int player_hitting(int i) {

	int temp_sum = player_summation(i);

	if (temp_sum > 21) {
		for (int j = 0; j < i; j++) {
			if (player_table[j] == 11) {
				player_table[j] = 1;
				return player_hitting(i);
			}
			else {
				return temp_sum;
			}
		}
	}
	else if (temp_sum < 17) {
		player_table[i - 1] = test_num(card_deq[i + 1]);
		return player_hitting(i + 1);
	}
	else {
		return temp_sum;
	}
}


int dealer_hitting(int i) {

	dealer_table[i-1] = test_num(card_deq[i + 15]);

	int temp_sum = dealer_summation(i);

	if (temp_sum > 21) {
		for (int j = 0; j < i; j++) {
			if (dealer_table[j] == 11) {
				dealer_table[j] = 1;
				return dealer_hitting(i);
			}
			else {
				return temp_sum;
			}
		}
	}
	else if (temp_sum < 17) {
		
		return dealer_hitting(i + 1);
	}
	else {
		return temp_sum;
	}
}


int black_jack() {
	// return 0 => loss
	// return 1 => win
	// return 2 => draw
	int play_final_sum = 0;
	int dealer_final_sum = 0;
	int temp_betting = betting;

	making_deq();
	clear_table();

	dealer_table[0] = test_num(card_deq[0]); // unknown

	player_table[0] = test_num(card_deq[1]);
	player_table[1] = test_num(card_deq[2]);

	
		if (winning_prob[player_table[0]][player_table[1]][dealer_table[0]] > 0.5) {
			if(player_summation(2) > 17){
				temp_betting = temp_betting * 2; // double!!
			}
		}

		play_final_sum = player_hitting(3);

		if (play_final_sum > 21) {
			initial_money = initial_money - temp_betting;
			return 0;
		}
		// Bust -> loss

		dealer_final_sum = dealer_hitting(2);

		if (dealer_final_sum > 21) {
			initial_money = initial_money + temp_betting;
			return 1;
		} // Bust -> win

		if (play_final_sum > dealer_final_sum) {
			initial_money = initial_money + temp_betting;
			return 1;
		}

		if (play_final_sum == dealer_final_sum) {
			return 2;
		}

		if (play_final_sum < dealer_final_sum) {
			initial_money = initial_money - temp_betting;
			return 0;
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
			player_prob_data[pl_first - 1][pl_second - 1][4] * 100) / 10000;

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
				//	cout << " Player : " << decide_card_name(i) << " & " << decide_card_name(j) << endl;
				//	cout << " Dealer : " << decide_card_name(k) << "	 Winning %  = " << winning_prob[i - 1][j - 1][k - 1] * 100 << " % " << endl;
				}
			}
		}
	}
}


void writing_winning_prob_data() {
	ofstream os;
	os.open("winning_prob");


	for (int i = 1; i <= 13; i++) {
		for (int j = 1; j <= 13; j++) {
			for (int k = 1; k <= 13; k++) {
				os << i << ends;
				os << j << ends;
				os << k << ends;
				os << winning_prob[i - 1][j - 1][k - 1] << ends;
			}
		}
	}

	cout << "writing complete" << endl;
}

void writing_simul_data() {
	ofstream os;
	os.open("simulation");

	for (int i = 0; i < trial; i++)
		os << game_result[i] << endl<< ends;

	cout << "writing complete" << endl;
}

double analyze_winning_percentage() {
	return (double)win_count / (double)(win_count + loss_count);
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

	ifstream input_file("dealer_prob_15000.txt");

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

	ifstream input_file_2("player_prob_15000.txt");

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
		//print_data();
	}

}


int main() {



	srand(time(NULL));

	cout << "Number of Trial? >> ";
	cin >> trial;

	cout << "Budget? >> ";
	cin >> initial_money;

	cout << "Set initial betting >> ";
	cin >> betting;

	money = initial_money;
	int result;


	auto time_1 = Clock::now();

	for (int j = 0; j<13; j++)
		for (int k = 0; k<13; k++)
			for (int i = 0; i < 6; i++)
				result_data[j][k][i] = 0;

	get_data();
	calculation();

	for (int k = 0; k < 10000; k++) {
		for (int i = 0; i < trial; i++) {
			result = black_jack();
			//game_result.push_back(initial_money);
			if (result == 0) loss_count++;
			else if (result == 1) win_count++;
			game_result.push_back(initial_money - money);
		}
		initial_money = money;
	}
	double average_profit = 0;
	for (int i = 0; i < game_result.size(); i++)
		average_profit += (double)game_result[i];

	average_profit = average_profit / (double)game_result.size();

	


	//writing_simul_data();

	cout << " Winning percentage is " << analyze_winning_percentage() * 100 << " %" << endl;
	//cout << " Profit is " << (double)initial_money - (double)money << " $ " <<endl;
	cout << " Average Profit is " << average_profit << " $ " << endl;
	
	auto time_2 = Clock::now();
	
	//writing_winning_prob_data();
	
	double time_diff_call = chrono::duration_cast<chrono::nanoseconds>(time_2 - time_1).count() / pow(10, 9);
	cout << " Time = " << time_diff_call << " sec" << endl;
}