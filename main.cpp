#include "gmm_fic.h"

int main(int argc, char* argv[]) {
  if (argc != 5) {
    cout << "The number of argument is invald." << endl;
    return(0);
  }
  
  double delta = atof(argv[1]);
  int numbet_of_state = (int)(1.0/delta);
  string input_file = argv[2];
  string output_parameter_file = argv[3];
  string output_responsibility_file = argv[4];
  
  srand((unsigned)time(NULL));
  run_EM(number_of_state, delta,input_file, output_parameter_file, output_responsibility_file);
  
}
