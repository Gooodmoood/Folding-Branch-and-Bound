// FoldingBnB.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "pch.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>

using namespace std;

struct Element  // Structure describes one monomer placed lattice
{
	int row, col;
	char value;
	
};
struct Near {
	int row, col;
	bool flag;
};
class Matrix    // lattice class 
{
	 
public:	
	vector<Element> data;  // vector of Elements
	int width, height;
	void printMatrixCoo()
		// COO matrix output
		{
			for (auto i = data.begin(); i < data.end(); i++)  
				cout << "Element [" << (*i).row << "," << (*i).col << "] = " << (*i).value << "\n";
		}
	void printMatrixStandart() {
		// Standart matrix output
		int flag = 0;
		{
			for (auto i = 0; i < width; i++)
			{
				for (auto j = 0; j < width; j++)
				{
					for (auto t = data.begin(); t < data.end(); t++) {
						if (((*t).row == i) && ((*t).col == j)) {
							cout << (*t).value << " ";  // вывод элементов
							flag = 1;
						}
					}
					if (flag == 0) {
						cout << 0 << " ";  // вывод элементов
					}
					flag = 0;

				}
				cout << endl;
			}// цикл по элементам

		}
	}
};
int Emin = 1;
int counter = 0;
vector<int> U{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
vector<int> Z{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
vector<int> Counter{ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
int print_matrix(int ** matrix, int N, int K) {
	int i, j;
	// Print matrix
	printf("Matrix: \n");
	for (i = 0; i < N; i++) {
		for (j = 0; j < K; j++) {
			printf(" %d", matrix[i][j]);
		}
		printf("\n");
	}
	printf("\n");


	return 0;
}
int generateSeqData() {
	// Generate protein HP - sequence and prints to file
	ofstream data;

	string seq1{"PHHPHP"};
	string seq2{ "HPPPHHHHHPPPPPP" };
	string seq3{ "HPPHHHPPPHPHHPHPHPHPPHHH" };
	string seq4{ "HHPHPHPHPPPHPHPHHHPPPPHPHP" };
	string seq5{ "HPHPPPPPHHHHPPPPHPHPHPH" };

	cout << "Generation starts" << std::endl;
	cout << "Filepath: ""D:\\Sequences.txt""" << std::endl;
	data.open("D:\\Sequences.txt"); 
	if (data.is_open())
	{
		data << seq1 << std::endl;
		data << seq2 << std::endl;
		data << seq3 << std::endl;
		data << seq4 << std::endl;
		data << seq5 << std::endl;
	}

	cout << "Sequences has been generated successfully" << std::endl;
	cout << "" << endl;
	return 0;

}
int calcEnergy(Matrix lattice) {

	// calculate energy over lattce
	int energy = 0;		//final enregy of lattice
	int corelator = 0;	// lattice inner energy dont need to be cunted. Result = energy + corelator;
	int minrow = lattice.data[0].row;
	int mincol = lattice.data[0].col;
	int maxcol = lattice.data[lattice.data.size() - 1].col;
	int maxrow = lattice.data[lattice.data.size() - 1].row;


	for (auto it = lattice.data.begin(); it < lattice.data.end(); ++it) {

		if (minrow > (*it).row ) {
			minrow = (*it).row;
		}
		if (mincol > (*it).col ) {
			mincol = (*it).col;
		}
		if (maxcol < (*it).col) {
			maxcol = (*it).col;
		}
		if (maxrow < (*it).row) {
			maxrow = (*it).row;
		}
	}


	int N = maxcol - mincol + 1;
	int N1 = maxrow - minrow + 1;
	int **matrix;
	int i;
	matrix = new int*[N1];

	// Matrix initialisation
	for (auto i = 0; i < N1; i++)
		matrix[i] = new int[N];

	for (auto i = 0; i < N1; i++) {
		for (auto j = 0; j < N; j++) {
			matrix[i][j] = 0;
		}
	}

	// define frames for searching energy as square [minrow .. maxcol]
	
	for (auto it = lattice.data.begin(); it < lattice.data.end(); ++it) {

		if ((*it).value == 'H') {
			matrix[(*it).row - minrow][(*it).col - mincol] = 1;
		}
		if ((*it).value == 'P') {
			matrix[(*it).row - minrow][(*it).col - mincol] = 2;
		}
		if (it > lattice.data.begin()) {
			if (((*it).value == 'P') && ((*(it - 1)).value == 'P')) {
				corelator++;
			}
		}
		
	}
	for (auto i = 0; i < N1; i++) {
		for (auto j = 0; j < N - 1; j++) {
			if ((matrix[i][j] == 2) and(matrix[i][j + 1] == 2)) {
				energy++;
			}
		}
	}
	for (auto j = 0; j < N; j++) {
		for (auto i = 0; i < N1 - 1; i++) {
			
			if ((matrix[i][j] == 2) and (matrix[i + 1][j] == 2)) {
				energy++;
			}
		}
	}

	

	return -(energy - corelator);
}
int search(Matrix lattice, string seq, int k) {
	// lattice - не пустой!
	Matrix tmp_lattice = lattice;
	Element tmp;
	Near tmp_pos;
	vector<Near> pos;
	int prev_col, prev_row;
	vector<Element> possible;
	char monomer_liter = seq[k];
	int tmp_energy;
	int avgZ;
	double r, ro1 = 0.8, ro2 = 0.5;
	std::random_device rd;  //Will be used to obtain a seed for the random number engine
	std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
	std::uniform_real_distribution<> dis(0.0, 1.0);
	if (k < 1) {
		int i, j;
		i = seq.length();
		j = seq.length();
		tmp.col = j;
		tmp.row = i;
		tmp.value = seq[0];
		lattice.data.push_back(tmp);
		tmp.col = j;
		tmp.row = i + 1;
		tmp.value =seq[1];
		lattice.data.push_back(tmp);
		k = 2;
	}
	if (k > 1) {
		tmp.value = monomer_liter;
		tmp.row = 0;
		tmp.col = 0;
		prev_col = lattice.data[lattice.data.size() - 1].col;
		prev_row = lattice.data[lattice.data.size() - 1].row;

		// This block defines possible coordinates for monomer
		//
		{
			tmp_pos.flag = true;
			tmp_pos.col = prev_col + 1;
			tmp_pos.row = prev_row;
			pos.push_back(tmp_pos);
			tmp_pos.col = prev_col - 1;
			tmp_pos.row = prev_row;
			pos.push_back(tmp_pos);
			tmp_pos.col = prev_col;
			tmp_pos.row = prev_row - 1;
			pos.push_back(tmp_pos);
			tmp_pos.col = prev_col;
			tmp_pos.row = prev_row + 1;
			pos.push_back(tmp_pos);


			// check for nearby position avaliability
			//

			for (auto it = lattice.data.begin(); it < lattice.data.end(); ++it) {
				//cout << "lattice.row: " << (*it).row << ", Col: " << (*it).col << endl;
				for (auto i = 0; i < 4; i++) {
					if ((pos[i].col == (*it).col) && (pos[i].row == (*it).row)) {
						pos[i].flag = false;
					};
				}
			}

			for (auto it = pos.begin(); it < pos.end(); it++) {
				if ((*it).flag == true) {
					tmp.col = (*it).col;
					tmp.row = (*it).row;
					possible.push_back(tmp);
				};
			}
		}

		// Calculate energy each possible place
		// 
		if (!possible.empty()) {
			// For - loop over all possible sites 
			for (auto it = possible.begin(); it < possible.end(); it++) {
				tmp_lattice.data.push_back((*it));
				
				tmp_energy = calcEnergy(tmp_lattice); // calculating energy
				Counter[k]++;
				Z[k] += tmp_energy;
				avgZ = int(Z[k] / Counter[k]);
				//cout << "tmp_energy = " << tmp_energy << endl;
				tmp_lattice.data.pop_back();
				
				// Record
				if (tmp_energy < U[k]) {
					U[k] = tmp_energy;
				}
				// if reached end - return Emin
				if (k == seq.length()) { // if end of the protein reached
					if (Emin > tmp_energy) {
						Emin = tmp_energy;
						//cout << "Emin: " << Emin << ", k: " << k << endl;
						//tmp_lattice.printMatrixCoo();
						//cout << endl;
						
						return Emin;
					}
					return Emin;

				}
				else {
					if( seq[k] == 'H') {
						if (tmp_energy < U[k]){
							tmp_lattice.data.push_back((*it));
							search(tmp_lattice, seq, k + 1);
						}
						if (tmp_energy > avgZ + 1) {
							r = dis(gen);
							
							if (ro1 > ro1) {
								tmp_lattice.data.push_back((*it));
								search(tmp_lattice, seq, k + 1);
							}
							
						}
						if ((tmp_energy >= U[k]) && (tmp_energy <= avgZ)) {
							r = dis(gen);
							
							if ( r > ro2) {
								tmp_lattice.data.push_back((*it));
								search(tmp_lattice, seq, k + 1);
							}

						}
					}
					else {
						tmp_lattice.data.push_back((*it));
						search(tmp_lattice, seq, k + 1);
						tmp_lattice.data.pop_back();
						
						
					}
				}

				
				
			}
		}
		

	}

};
int BNB(string seq)
{


	 //Testing 2
	
	Matrix lattice;
	Element tmp;

	lattice.height = seq.length() * 2;
	lattice.width = lattice.height;
	

	std::cout << "Protein : " << seq << endl;
	search(lattice, seq, 0);
	cout << "Minimal conformation energy:" << Emin << endl;





	

	

	return 0;
}
int main()
{	
	
	string seq;
	ifstream in("D:\\Sequences.txt"); // Data file path
	int len, Mlen;
	int i, j;
	vector<string> seqData;// vector of sequences
	generateSeqData();
	// read from file to vector
	if (in.is_open())
	{
		while (getline(in, seq))
		{
			seqData.push_back(seq);
		}
	}
	
	BNB(seqData[0]);

	return 0;
}