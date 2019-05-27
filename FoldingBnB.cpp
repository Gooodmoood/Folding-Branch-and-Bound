// FoldingBnB.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <random>
#include <omp.h>
#include <time.h>
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
struct Monomer {
	char Liter;
	int U = 0;
	int Z = 0;
	int Counter = 0;
	int InnerEnergy = 0;
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
	void printMatrixStandartMk2() {
		// Standart matrix output
		int flag = 0;
		int min_row = 999;
		int max_col = 0;
		int tmp_max = 0;
		int tmp_min = 0;

		
		for (auto i = 0; i < width; i++)
		{
			for (auto j = 0; j < width; j++)
			{
				for (auto t = data.begin(); t < data.end(); t++) {
					if (min_row > (*t).row) {
						min_row = (*t).row;
					}
					if (max_col < (*t).col) {
						max_col = (*t).col;
					}
				}
			}
		}
		
		tmp_max = std::max(min_row, max_col)+2;
		tmp_min = std::min(min_row, max_col)-1;

		for (auto i = tmp_min; i < tmp_max; i++)
		{
			for (auto j = tmp_min; j < tmp_max; j++)
			{
				for (auto t = data.begin(); t < data.end(); t++) {
					if (((*t).row == i) && ((*t).col == j)) {
						cout << (*t).value << " ";  // вывод элементов
						flag = 1;
					}
				}
				if (flag == 0) {
					cout << 0 << " ";
				}
				flag = 0;

			}
			cout << endl;
		}

		
	}

};
class Protein {
public:
	vector<Monomer> data;
	int a;
};
int Emin = 1;
vector<int> EnergyMin;
void omp_pc_info() {
	int np = omp_get_max_threads();
	printf("OMP get max threads: %d. \n", np);
	omp_set_num_threads(np);
}
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

	string seq1{ "HPHPPHHPHPPHPHP"};
	string seq2{ "HPPH" };
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
int calcEnergy(Matrix lattice, Protein protein) {

	// calculate energy over lattce
	int energy = 0;		//final enregy of lattice
	int corelator = 0;	// lattice inner energy dont need to be cunted. Result = energy + corelator;
	int minrow = lattice.data[0].row;
	int mincol = lattice.data[0].col;
	int maxcol = lattice.data[lattice.data.size() - 1].col;
	int maxrow = lattice.data[lattice.data.size() - 1].row;
	int k;
	
	if (lattice.data.size() - 1 >= 2) {
		k = lattice.data.size() - 1;
	}
	else {
		k = 0;
	}
	
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
			matrix[(*it).row - minrow][(*it).col - mincol] = 2;
		}
		if ((*it).value == 'P') {
			matrix[(*it).row - minrow][(*it).col - mincol] = 1;
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


//	print_matrix(matrix, N1, N);
	return protein.data[k].InnerEnergy - energy;
	//return energy;
}
int search(Matrix lattice, string seq, Protein protein, int Emin) {
	// lattice - не пустой!

	Element tmp;
	Near tmp_pos;
	vector<Near> pos;
	int prev_col, prev_row;
	int k = (lattice.data.end() - lattice.data.begin());
	vector<Element> possible;
	
	//char monomer_liter = seq[k];
	int tmp_energy;
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
		tmp.value = seq[1];
		lattice.data.push_back(tmp);
		k = 2;
	}
	if (k > 1) {
		tmp.value = seq[k];
		tmp.row = 0;
		tmp.col = 0;
		prev_col = lattice.data[lattice.data.size() - 1].col;
		prev_row = lattice.data[lattice.data.size() - 1].row;

		// This block defines possible coordinates for monomer
		//
		{
			pos.clear();
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
			possible.clear();
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

			//#pragma omp parallel for if (k = 2)
//			#pragma omp single nowait
			for (int it = 0; it < int(possible.size()); it++) {
//				#pragma omp task
				{
					lattice.data.push_back(possible[it]);
					tmp_energy = calcEnergy(lattice, protein); // calculating energy
					lattice.data.pop_back();

					// if reached end - return Emin
					if (k == seq.length() - 1) { // if end of the protein reached
						if (Emin > tmp_energy) {
							Emin = tmp_energy;
							cout << "Emin: " << Emin << endl;
							lattice.data.push_back(possible[it]);
							lattice.printMatrixCoo();
							lattice.printMatrixStandartMk2();
							lattice.data.pop_back();
							
						}
						return 0;

					}
					else {
						if (seq[k] == 'H') {

							lattice.data.push_back(possible[it]);
							search(lattice, seq, protein, Emin);
							lattice.data.pop_back();


						}
						else {
							lattice.data.push_back(possible[it]);
							search(lattice, seq, protein, Emin);
							lattice.data.pop_back();


						}

					}
				}
			}
		}

		lattice.data.pop_back();
	}

	return 0;
};

int searchBnB_ser(Matrix lattice, string seq, Protein protein) {
	// lattice - не пустой!

	Element tmp;
	Near tmp_pos;
	vector<Near> pos;
	int prev_col, prev_row;
	int k = (lattice.data.end() - lattice.data.begin());
	vector<Element> possible;
	//char monomer_liter = seq[k];
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
		tmp.value = seq[1];
		lattice.data.push_back(tmp);
		k = 2;
	}
	if (k > 1) {
		tmp.value = seq[k];
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
				//cout << omp_get_thread_num() << endl;
				lattice.data.push_back((*it));
				tmp_energy = calcEnergy(lattice, protein); // calculating energy
				lattice.data.pop_back();
				protein.data[k].Counter++;
				protein.data[k].Z += tmp_energy;
				avgZ = int(protein.data[k].Z / protein.data[k].Counter);

				// Record
				if (tmp_energy < protein.data[k].U) {
					protein.data[k].U += tmp_energy;
				}
				// if reached end - return Emin
				if (k == seq.length() - 1) { // if end of the protein reached
					if (Emin > tmp_energy) {
						Emin = tmp_energy;
						cout << "Emin: " << Emin << endl;
						EnergyMin.push_back(Emin);
						//lattice.data.push_back((*it));
						//lattice.printMatrixCoo();
						//lattice.printMatrixStandartMk2();
						//lattice.data.pop_back();
						return 0;
					}
					return 0;
				}
				else {
					if (seq[k] == 'H') {

						if (tmp_energy < protein.data[k].U){
							lattice.data.push_back((*it));
							searchBnB_ser(lattice, seq, protein);
							lattice.data.pop_back();
							
						}
						if (tmp_energy > avgZ + 1) {
							r = dis(gen);
							
							if (ro1 > ro1) {
								lattice.data.push_back((*it));
								searchBnB_ser(lattice, seq, protein);
								lattice.data.pop_back();
								
							}
							
						}
						if ((tmp_energy >= protein.data[k].U) && (tmp_energy <= avgZ)) {
							r = dis(gen);
							
							if ( r > ro2) {
								lattice.data.push_back((*it));
								searchBnB_ser(lattice, seq, protein);
								lattice.data.pop_back();
							}
						}
					}
					else {
						lattice.data.push_back((*it));
						searchBnB_ser(lattice, seq, protein);
						lattice.data.pop_back();


					}
				}
			
			}
		}
		lattice.data.pop_back();
	}
	return 0;
};

int searchBnB(Matrix lattice, string seq, Protein protein) {
	
	// lattice - не пустой!
	Element tmp;
	Near tmp_pos;
	vector<Near> pos;
	int prev_col, prev_row;
	int k = (lattice.data.end() - lattice.data.begin());
	vector<Element> possible;
	//char monomer_liter = seq[k];
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
		tmp.value = seq[1];
		lattice.data.push_back(tmp);
		k = 2;
	}
	if (k > 1) {
		tmp.value = seq[k];
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

			
			for (int it = 0; it < possible.size(); it++) {
				
				#pragma omp task
				{
					cout << "Thread # " << omp_get_thread_num() << endl;
					lattice.data.push_back(possible[it]);
					tmp_energy = calcEnergy(lattice, protein); // calculating energy
					lattice.data.pop_back();
					protein.data[k].Counter++;
					protein.data[k].Z += tmp_energy;
					avgZ = int(protein.data[k].Z / protein.data[k].Counter);

					// Record
					if (tmp_energy < protein.data[k].U) {
						protein.data[k].U += tmp_energy;
					}
					// if reached end - return Emin
					if (k == seq.length() - 1) { // if end of the protein reached
						if (Emin > tmp_energy) {
							Emin = tmp_energy;
							//cout << "Emin: " << Emin << endl;
							EnergyMin.push_back(Emin);
							//lattice.data.push_back((*it));
							//lattice.printMatrixCoo();
							//lattice.printMatrixStandartMk2();
							//lattice.data.pop_back();
							//return 0;
						}
						//return 0;
						
					}
					else {
						if (seq[k] == 'H') {

							if (tmp_energy < protein.data[k].U){
								lattice.data.push_back(possible[it]);
								searchBnB_ser(lattice, seq, protein);
								lattice.data.pop_back();
								
							}
							if (tmp_energy > avgZ + 1) {
								r = dis(gen);
								
								if (ro1 > ro1) {
									lattice.data.push_back(possible[it]);
									searchBnB_ser(lattice, seq, protein);
									lattice.data.pop_back();
									
								}
								
							}
							if ((tmp_energy >= protein.data[k].U) && (tmp_energy <= avgZ)) {
								r = dis(gen);
								
								if ( r > ro2) {
									lattice.data.push_back(possible[it]);
									searchBnB_ser(lattice, seq, protein);
									lattice.data.pop_back();
								}
							}
						}
						else {
							lattice.data.push_back(possible[it]);
							searchBnB_ser(lattice, seq, protein);
							lattice.data.pop_back();


						}
					}				
				}
				continue;
			}
		return 0;
		}
		lattice.data.pop_back();
	}
	return 0;
};

string Hello() {
	return "Hello";
};
string World() {
	return "World";
};
int omp_test() {
	string str1, str2;
//	#pragma omp single nowait
	{
//		#pragma omp task
		{
			str1 = Hello();
		}
//		#pragma omp task
		{
			str2 = World();
		}
	}
	
	cout << str1 << " " << str2 << endl;
	return 0;
};
int BNB(string seq)
{


	
	omp_set_num_threads(omp_get_num_procs());
	cout << "Thrads: " << omp_get_num_procs() << endl;
	Monomer monomer;
	Matrix lattice;
	Protein protein;
	int Ans;
	int Emin = 1;

	lattice.height = seq.length() * 2;
	lattice.width = lattice.height;
	for (auto it = seq.begin(); it < seq.end(); it++) {
		monomer.Liter = (*it);
		monomer.U = 0;
		monomer.Z = 0;
		monomer.Counter = 0;
		monomer.InnerEnergy = 0;
		for (auto it2 = seq.begin(); it2 < it+1; it2++ ) {

			if (it2 > seq.begin()) {
				if (((*it2) == 'H') && ((*(it2 - 1)) == 'H')) {
					monomer.InnerEnergy++;
				}
			}
		}
		protein.data.push_back(monomer);

	}
	std::cout << "Protein : " << seq << endl;
	
	//cout << "Straight algorithm"<< endl;
	//search(lattice, seq, protein);
	
	//cout << "Minimal conformation energy: " << Emin << " Time spent(sec): "<< t << endl
	int avgTime = 0;
	int i;
	int t;
	// serial version
	EnergyMin.clear();
	cout <<"Bnb parallel"<< endl;
	t = time(NULL);
	{	
		#pragma omp parallel
		#pragma omp single nowait
		searchBnB(lattice, seq, protein);
	}
	t = time(NULL) - t;
	cout <<" Time spent(sec): " << t << endl;
	Ans = 999;
	for (auto it = EnergyMin.begin(); it < EnergyMin.end(); it++) {
		//std::cout << (*it) << endl;
		if (Ans > (*it))
		{
			Ans = (*it);
		}
	}
	std::cout << "Energy: " << Ans << endl;

	// parallel version
	//EnergyMin.clear();
	//Emin = 1;
	//cout <<"Bnb parallel"<< endl;
	//t = time(NULL);

	//{
	//	searchBnB(lattice, seq, protein);
	//}
	//t = time(NULL) - t;
	//cout <<" Time spent(sec): " << t << endl;
	//Ans = 999;
	//for (auto it = EnergyMin.begin(); it < EnergyMin.end(); it++) {
		//std::cout << (*it) << endl;
	//	if (Ans > (*it))
	//	{
	//		Ans = (*it);
	//	}
	//}
	//std::cout << "Energy: " << Ans << endl;

	return 0;



};
int main()
{	
	
	
	string seq;
	ifstream in("D:\\Sequences.txt"); // Data file path
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
	//omp_test();


	return 0;
}
