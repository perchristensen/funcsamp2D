funcsamp2D.cpp program by Per Christensen, Pixar, 2019.

This directory contains the following files:

funcsamp2D.cpp: C++ source code for a program that reads a file with
sample sequences, uses the samples to sample a specified 2D function,
and writes out sampling error for increasing numbers of samples.

users_guide.pdf: A user's guide for the funcsamp2D program.

Examples of input files: 
random_1024samples_100sequences.data
bestcand_1024samples_100sequences.data
irrational_rot_1024samples_100sequences.data
halton_base23_owen_1024samples_100sequences.data
halton_base57_owen_1024samples_100sequences.data
halton_base1113_owen_1024samples_100sequences.data
halton_base1719_owen_1024samples_100sequences.data
halton_base2329_owen_1024samples_100sequences.data
sobol_dim01_owen_1024samples_100sequences.data
sobol_dim23_owen_1024samples_100sequences.data
sobol_dim45_owen_1024samples_100sequences.data
sobol_dim67_owen_1024samples_100sequences.data
sobol_dim89_owen_1024samples_100sequences.data
pmj02_1024samples_100sequences.data

Examples of output files:
errors_quarterdisk_bestcand.data
errors_quarterdisk_halton_base23_owen.data
errors_quarterdisk_irrational_rot.data
errors_quarterdisk_pmj02.data
errors_quarterdisk_random.data
errors_quarterdisk_sobol_dim01_owen.data

plotQuarterdiskError.gp: A gnuplot script that plots the output results
as curves.

quarterdiskError.eps,quarterdiskError.pdf: Examples of gnuplot output.

