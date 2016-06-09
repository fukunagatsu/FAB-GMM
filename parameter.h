#include <Eigen/Core>
#include <Eigen/LU>
#include <math.h>
#include <vector>
using namespace std;

class parameter{
private:
	static const double m_pi = 3.14159265358979323846;

	double box_muller(double mu, double var,int dimension) {
		int k = 2;
		vector<double> random_number_array; random_number_array.reserve(k);
		for (int i = 0; i<k; i++) {
			double a = (double)rand() / ((double)RAND_MAX + 1);
			random_number_array.push_back(a);
		}

		for (int i = 0; i<k; i += 2) {
			double a = 0; double b = 0;
			a = mu + sqrt(-2 * log(random_number_array[i])) * sin(2 * m_pi*random_number_array[i + 1]) * var;
			b = mu + sqrt(-2 * log(random_number_array[i])) * cos(2 * m_pi*random_number_array[i + 1]) * var;
			random_number_array[i] = a;
			random_number_array[i + 1] = b;
		}

		return(random_number_array[0]);
	}

public:
	vector<double> pi_array;
	vector<Eigen::VectorXd> average_array;
	vector<Eigen::MatrixXd> covariance_matrix_array;
	int number_of_state;
	int dimension;
	int number_of_data;

	double calc_likelihood(vector<vector<double> >& dataset) {
		double LL = 0;
		vector<Eigen::MatrixXd> inversed_covariance_matrix_array; inversed_covariance_matrix_array.reserve(number_of_state);
		vector<double> determinant_array; inversed_covariance_matrix_array.reserve(number_of_state);
		
		for (int i = 0; i < number_of_state; i++) {
			if (pi_array[i] != 0) {
				Eigen::MatrixXd temp_matrix = covariance_matrix_array[i].inverse();
				inversed_covariance_matrix_array.push_back(temp_matrix);
				determinant_array.push_back(covariance_matrix_array[i].determinant());
			}else {
				Eigen::MatrixXd temp_matrix = Eigen::MatrixXd::Zero(dimension, dimension);
				inversed_covariance_matrix_array.push_back(temp_matrix);
				determinant_array.push_back(0.0);
			}
		}
		for (int i = 0; i < number_of_data; i++) {
			double a = -1000000000;
			for (int j = 0; j < number_of_state; j++) {
				if (pi_array[j] != 0) {
					double exp_in = (-1.0) / 2;
					Eigen::VectorXd temp_vector = Eigen::VectorXd::Zero(dimension);
					for (int k = 0; k < dimension; k++) {
						temp_vector(k) = dataset[i][k] - average_array[j](k);
					}
					exp_in *= temp_vector.transpose()*inversed_covariance_matrix_array[j] * temp_vector;

					if (a == -1000000000) {
						a = exp_in + log(pi_array[j]) - (log(pow(2 * m_pi, dimension / 2)) + log(sqrt(determinant_array[j])));
					}
					else {
						a = logsum(a, exp_in + log(pi_array[j]) - (log(pow(2 * m_pi, dimension / 2)) + log(sqrt(determinant_array[j]))));
					}

				}
			}
			LL += a;
			
		}

		return(LL);
	}

	double logsum(double x, double y) {
		if (x>y) {
			return((x + log(exp(y - x) + 1.0)));
		}
		else {
			return((y + log(exp(x - y) + 1.0)));
		}
	}

	void initiallize(vector<vector<double> >& dataset, int temp_number_of_state, int temp_dimension, int  temp_number_of_data){
		number_of_state = temp_number_of_state;
		dimension = temp_dimension;
		number_of_data = temp_number_of_data;
		pi_array.reserve(number_of_state);
		average_array.reserve(number_of_state);
		covariance_matrix_array.reserve(number_of_state);

		//calc mean and variance
		vector<double> mu; mu.reserve(dimension);
		vector<double> var; var.reserve(dimension);
		for (int i = 0; i < dimension; i++) {
			mu.push_back(0.0); var.push_back(0.0);
		}

		for (int i = 0; i < number_of_data; i++) {
			for (int j = 0; j < dimension; j++) {
				mu[j] += dataset[i][j];
			}
		}
		for (int i = 0; i < dimension; i++) {
			mu[i] /= number_of_data;
		}

		for (int i = 0; i < number_of_data; i++) {
			for (int j = 0; j < dimension; j++) {
				var[j] += (dataset[i][j] - mu[j]) * (dataset[i][j] - mu[j]);
			}
		}
		for (int i = 0; i < dimension; i++) {
			var[i] /= number_of_data;
		}

		
		//initiallize parameter
		for (int i = 0; i < number_of_state; i++) {
			Eigen::VectorXd temp_vector = Eigen::VectorXd::Zero(dimension);
			average_array.push_back(temp_vector);
		}

		for (int i = 0; i < number_of_state; i++) {
			for (int j = 0; j < dimension; j++) {
				average_array[i](j) = box_muller(mu[j],sqrt(var[j]),dimension);
			}
		}

		for (int i = 0; i < number_of_state; i++) {
			pi_array.push_back(1.0/number_of_state);
		}

		for (int i = 0; i < number_of_state; i++) {
			Eigen::MatrixXd temp_matrix = Eigen::MatrixXd::Identity(dimension,dimension);
			covariance_matrix_array.push_back(temp_matrix);
		}
	}
};
