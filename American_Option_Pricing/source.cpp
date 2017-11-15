#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>

using namespace std;

int no_of_divisions;
float **memo;


void make_memo() {
	memo = new float *[no_of_divisions];
		for (int i = 0; i <= no_of_divisions; i++)
		memo[i] = new float[(2*i+1)];

		cout <<endl<< "					 MAKE memo complete@@@@@@@@@@@@@@@@@@@@" << endl;
}


void clean_memo(float **memo_array) {			//왜 Array 크기를 0으로만들어도 Array가 생기는지 질문!
// 왜 Array를 만들면, -1, -2, -3 에도 자리가 생기는지 질문.
		for (int i = 0; i <= no_of_divisions; i++)
			for (int j = 0; j <= 2 * i + 1; j++)  // put zeros to all.
				memo_array[i][j] = -1;

		cout <<endl<< "					 CLEAN memo complete@@@@@@@@@@@@@@@@@@@@" << endl;
}

void print_memo(float **memo_array) {
		for (int i = 0; i < no_of_divisions; i++){
			for (int j = 0; j < 2 * i + 1; j++) {  // put zeros to all.
			cout << memo_array[i][j] << " ";
			}
		cout << endl;
		}
	cout << " 메모 0,-4      " << memo[1][-4] << endl;  // 이렇게 정의가 안되어야하는 곳에도 공간이 있어!
}


int main(int argc, char* argv[])
{

	sscanf(argv[1], "%d", &no_of_divisions);
	
	make_memo();
	clean_memo(memo);
	print_memo(memo);

	return 0;
}

