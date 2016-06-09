#include "gmm_fic.h"

void run_EM(int number_of_state, double delta, string input_file, string output_parameter_file, string output_responsibility_file) {
	gmm temp_object(number_of_state, delta,input_file, output_parameter_file, output_responsibility_file);
	temp_object.load_data();
	temp_object.run();
}