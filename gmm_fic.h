#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <time.h>
#include "parameter.h"

using namespace std;

class gmm {
private:
	int number_of_state;
	int number_of_data;
	int dimension;
	int dc;
	double delta;
	static const double m_pi = 3.14159265358979323846;
	string input_file_name;
	string output_parameter_file_name;
	string output_responsibility_file_name;
	vector<vector<double> > dataset;
	vector<vector<double> > responsibility;
	vector<int> shrinkage_array;

public:
	gmm(int x, double y,string input_file, string output_parameter_file, string output_responsibility_file) {
		number_of_state = x;
		delta = y;
		input_file_name = input_file;
		output_parameter_file_name = output_parameter_file;
		output_responsibility_file_name = output_responsibility_file;
	}

	void run() {
		parameter old_parameter;
		old_parameter.initiallize(dataset, number_of_state, dimension, number_of_data);
		double old_LL = old_parameter.calc_likelihood(dataset);
		responsibility.reserve(number_of_data);
		for (int i = 0; i < number_of_data; i++) {
			vector<double> temp_vector; temp_vector.reserve(number_of_state);
			for (int j = 0; j < number_of_state; j++) {
				temp_vector.push_back(0.0);
			}
			responsibility.push_back(temp_vector);
		}
		
		E_step(old_parameter);
		
		parameter new_parameter;
		new_parameter = M_step();
		double new_LL = new_parameter.calc_likelihood(dataset);
		int loop_count = 0;
		while ((loop_count < 100)) {
			old_LL = new_LL;
			old_parameter = new_parameter;
			
			E_step(old_parameter);
			shrinkage();
			new_parameter = M_step();

			new_LL = new_parameter.calc_likelihood(dataset);
			loop_count++;
		}
		output(new_parameter);
	}

	void shrinkage() {
		for (int i = 0; i < number_of_state; i++) {
			double sum = 0.0;
			for (int j = 0; j < number_of_data; j++) {
				sum += responsibility[j][i];
			}
			if (sum < number_of_data*delta) {
				shrinkage_array[i] = 1;
			}
		}
		return;
	}

	void E_step(parameter old_parameter){
		vector<Eigen::MatrixXd> inversed_covariance_matrix_array; inversed_covariance_matrix_array.reserve(number_of_state);
		vector<double> determinant_array; inversed_covariance_matrix_array.reserve(number_of_state);
		
		for (int i = 0; i < number_of_state; i++) {
			if (shrinkage_array[i] == 0) {
				Eigen::MatrixXd temp_matrix = old_parameter.covariance_matrix_array[i].inverse();
				inversed_covariance_matrix_array.push_back(temp_matrix);
				determinant_array.push_back(old_parameter.covariance_matrix_array[i].determinant());
			}else {
				Eigen::MatrixXd temp_matrix = Eigen::MatrixXd::Zero(dimension, dimension);
				inversed_covariance_matrix_array.push_back(temp_matrix);
				determinant_array.push_back(0.0);
			}
		}
		
		for (int i = 0; i < number_of_data; i++) {
			double denominator = -1000000000;
			for (int j = 0; j < number_of_state; j++) {
				if (shrinkage_array[j] == 0) {
					double exp_in = (-1.0) / 2;
					Eigen::VectorXd temp_vector = Eigen::VectorXd::Zero(dimension);
					for (int k = 0; k < dimension; k++) {
						temp_vector(k) = dataset[i][k] - old_parameter.average_array[j](k);
					}
					exp_in *= temp_vector.transpose()*inversed_covariance_matrix_array[j] * temp_vector;

					if (denominator == -1000000000) {
						denominator = exp_in + log(old_parameter.pi_array[j]) + (-dc / (2 * old_parameter.pi_array[j] * number_of_data)) - (log(pow(2 * m_pi, dimension / 2)) + log(sqrt(determinant_array[j])));
					}
					else {
						denominator = logsum(denominator, exp_in + log(old_parameter.pi_array[j]) - dc / (2 * old_parameter.pi_array[j] * number_of_data) - (log(pow(2 * m_pi, dimension / 2)) + log(sqrt(determinant_array[j]))));
					}
				}
			}
			
			for (int j = 0; j < number_of_state; j++) {
				if (shrinkage_array[j] == 0) {
					double exp_in = (-1.0) / 2;
					Eigen::VectorXd temp_vector = Eigen::VectorXd::Zero(dimension);
					for (int k = 0; k < dimension; k++) {
						temp_vector(k) = dataset[i][k] - old_parameter.average_array[j](k);
					}
					exp_in *= temp_vector.transpose()*inversed_covariance_matrix_array[j] * temp_vector;

					double numerator = exp_in + log(old_parameter.pi_array[j]) + (-dc / (2 * old_parameter.pi_array[j] * number_of_data)) - (log(pow(2 * m_pi, dimension / 2)) + log(sqrt(determinant_array[j])));
					responsibility[i][j] = exp(numerator - denominator);
				}
			}
		}

		return;
	}

	double logsum(double x, double y) {
		if(x>y){
			return((x+log(exp(y-x) + 1.0)));
		}else {
			return((y + log(exp(x - y) + 1.0)));
		}
	}

	parameter M_step() {
		parameter new_parameter;
		new_parameter.dimension = dimension;
		new_parameter.number_of_data = number_of_data;
		new_parameter.number_of_state = number_of_state;
		vector<double> Nk; Nk.reserve(number_of_state);

		for (int i = 0; i<number_of_state; i++) {
			Nk.push_back(0);
			if (shrinkage_array[i] == 0) {
				for (int j = 0; j < number_of_data; j++) {
					Nk[i] += responsibility[j][i];
				}
			}
		}

		vector<double> pi; pi.reserve(number_of_state);
		for (int i = 0; i< number_of_state; i++) {
			pi.push_back(Nk[i] / number_of_data);
		}
		new_parameter.pi_array = pi;

		vector<Eigen::VectorXd> average_array;average_array.reserve(number_of_state);

		for (int i = 0; i < number_of_state; i++) {
			Eigen::VectorXd temp_vector = Eigen::VectorXd::Zero(dimension);
			if (shrinkage_array[i] == 0) {
				for (int j = 0; j < dimension; j++) {
					double beta = 0;
					for (int k = 0; k < number_of_data; k++) {
						beta += responsibility[k][i] * dataset[k][j];
					}
					beta = beta / Nk[i];
					temp_vector(j) = beta;
				}
			}
			average_array.push_back(temp_vector);
		}
		new_parameter.average_array = average_array;

		vector<Eigen::MatrixXd> covariance_matrix_array; covariance_matrix_array.reserve(number_of_state);

		vector<Eigen::VectorXd> mu; mu.reserve(number_of_state);
		for (int i = 0; i< number_of_state; i++) {
			Eigen::VectorXd temp_vector = Eigen::VectorXd::Zero(dimension);
			for (int j = 0; j < dimension; j++) {
				temp_vector(j) = average_array[i][j];
			}
			mu.push_back(temp_vector);
		}

		for (int i = 0; i < number_of_state; i++) {
			Eigen::MatrixXd temp_matrix = Eigen::MatrixXd::Zero(dimension, dimension);

			if (shrinkage_array[i] == 0) {
				for (int j = 0; j < number_of_data; j++) {
					Eigen::VectorXd temp_vector = Eigen::VectorXd::Zero(dimension);
					for (int k = 0; k < dimension; k++) {
						temp_vector(k) = dataset[j][k] - mu[i](k);
					}
					temp_matrix += responsibility[j][i] * (temp_vector * temp_vector.transpose());
				}
				temp_matrix /= Nk[i];
			}
			covariance_matrix_array.push_back(temp_matrix);
		}

		new_parameter.covariance_matrix_array = covariance_matrix_array;
		return(new_parameter);
	}

	void output(parameter temp_parameter) {
		ofstream ofs(output_parameter_file_name.c_str());
		for (int i = 0; i < number_of_state; i++) {
			if (temp_parameter.pi_array[i] != 0) {
				ofs << temp_parameter.pi_array[i] << " ";
			}
		}
		ofs << endl;
		ofs << endl;

		for (int i = 0; i < number_of_state; i++) {
			if (temp_parameter.pi_array[i] != 0) {
				for (int j = 0; j < dimension; j++) {
					ofs << temp_parameter.average_array[i][j] << " ";
				}
				ofs << endl;

				for (int j = 0; j < dimension; j++) {
					for (int k = 0; k < dimension; k++) {
						ofs << temp_parameter.covariance_matrix_array[i](j, k) << " ";
					}
					ofs << endl;
				}
				ofs << endl;
			}
		}

		ofs.close();

		ofstream ofs2(output_responsibility_file_name.c_str());

		for (int i = 0; i < number_of_data; i++) {
			for (int j = 0; j < number_of_state; j++) {
				if (temp_parameter.pi_array[j] != 0) {
					ofs2 << responsibility[i][j] << " ";
				}
			}
			ofs2 << endl;
		}

		ofs2.close();
	}

	void load_data() {
		ifstream fp;
		fp.open(input_file_name.c_str(), ios::in);
		if (!fp) {
			cout << "Cannot open " + input_file_name << endl;
			exit(1);
		}
		char buf[100000];
		fp.getline(buf, 100000);
		char* tp;
		tp = strtok(buf, " ");
		number_of_data = atoi(tp);
		tp = strtok(NULL, " ");
		dimension = atoi(tp);

		dc = dimension + (dimension * 2 + dimension) / 2;
		shrinkage_array.reserve(number_of_data);
		for (int i = 0; i < number_of_data; i++) {
			shrinkage_array.push_back(0.0);
		}

		int count = 0;
		for (int i = 0; i < number_of_data; i++) {
			fp.getline(buf, 100000);
			vector<double> temp_vector; temp_vector.reserve(dimension);

			tp = strtok(buf, " ");
			temp_vector.push_back(atof(tp));

			for (int i = 1; i < dimension; i++) {
				tp = strtok(NULL, " ");
				temp_vector.push_back(atof(tp));
			}
			dataset.push_back(temp_vector);
			count++;
		}
		fp.close();
	}

};

void run_EM(int number_of_state, double delta,string input_file, string output_parameter_file, string output_responsibility_file);