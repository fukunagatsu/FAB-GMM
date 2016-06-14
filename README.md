# FAB-GMM
FAB-GMM estimates Gaussian Mixture Model (GMM) parameter based on Factorized Asymptotic Bayes (FAB) algorithm. 
The FAB algorithm is similar to the conventional expecation-maximization algorithm for fitting GMM but allows automatic estimation of the numbers of mixture components based on Factorized Information Criteria (FIC).

## Requirements
* Eigen3.2.5

You have to rewrite the location of eigen3.2.5 in makefile accordingly.

##Usage
    ./FAB-GMM <epsilon> <input_file> <output_parameter_file> <output_responsibility_file>

<epsilon>: The minimum value of mixing ratio of a component. If mixing ratio of a component is smaller than <epsilon> during calculation, FAB algorithm shirnks the component.
<input_file>      : Name of input file.
<output_parameter_file>  : Name of output file in which estimated parameter is written.
<output_responsibility_file>  : Name of output file in which estimated responsibility is written.

##Format of input file
The first line describes the number of data and dimension.
After first line, each line describes each data.

Example:
6 4  
0.1 0.1 0.1 0.1  
0.2 0.1 0.3 0.3  
-0.1 0.2 0.1 -0.1  
0.1 -0.3 0.2 0.1  
0.1 0.1 -0.1 0.2  
-0.1 0.3 0.1 0.3  

##Format of output parameter file
The first line describes the mixing ratio of each component.
After first line, the vector of mean value and covariance matrix of each component are described.
Blank line is inserted between elements (vector of mean value and covariance matrix) of each component.

## Reference
Ryohei Fujimaki, and Satoshi Morinaga. "Factorized asymptotic bayesian inference for mixture modeling." International Conference on Artificial Intelligence and Statistics. 2012.

Tsukasa Fukunaga, and Wataru Iwasaki. "Importance of considering simple factors in *C. elegans*  behavioral analysis." under submission.
