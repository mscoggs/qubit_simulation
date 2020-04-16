CC=g++


compile:
	$(CC) -o main main.cpp adiabatic.cpp check.cpp export_data.cpp hamiltonian.cpp linear_algebra.cpp mcbb.cpp mcbf.cpp mcdb.cpp operations.cpp print.cpp parameters.cpp -lgsl -llapack -lblas -std=gnu++11


submit:
	$(CC) -o cluster/make_submit cluster/make_submit.cpp -std=gnu++11
	cluster/make_submit
	condor_submit submit_file.txt

compile_gap:
	$(CC) -o gap gap.cpp hamiltonian.cpp check.cpp export_data.cpp operations.cpp parameters.cpp print.cpp linear_algebra.cpp -lgsl -llapack -lblas -std=gnu++11
