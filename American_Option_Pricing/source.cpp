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


void clean_memo(float **memo_array) {			//�� Array ũ�⸦ 0���θ��� Array�� ������� ����!
// �� Array�� �����, -1, -2, -3 ���� �ڸ��� ������� ����.
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
	cout << " �޸� 0,-4      " << memo[1][-4] << endl;  // �̷��� ���ǰ� �ȵǾ���ϴ� ������ ������ �־�!
}


int main(int argc, char* argv[])
{

	sscanf(argv[1], "%d", &no_of_divisions);
	
	make_memo();
	clean_memo(memo);
	print_memo(memo);

	return 0;
}

