CC=g++

compile:
	$(CC) -o main main.cpp adiabatic.cpp check.cpp write_data.cpp hamiltonian.cpp linear_algebra.cpp mcbb.cpp mcbf.cpp mcdb.cpp operations.cpp print.cpp parameters.cpp -lgsl -llapack -lblas -std=gnu++11

submit:
	$(CC) -o cluster/make_submit cluster/make_submit.cpp -std=gnu++11
	cluster/make_submit
	condor_submit submit_file.txt

submit_unfinished:
	python cluster/make_submit_unfinished.py
	condor_submit submit_file.txt

submit_q:
	$(CC) -o q_cluster/make_q_submits q_cluster/make_q_submits.cpp -std=gnu++11
	q_cluster/make_q_submits 10
	chmod 700 submit_q_file.txt
	./submit_q_file.txt

fixed_protocol:
	$(CC) -o fixed fixed_mcbb_protocol.cpp check.cpp write_data.cpp hamiltonian.cpp linear_algebra.cpp operations.cpp print.cpp parameters.cpp -lgsl -llapack -lblas -std=gnu++11
	./fixed

compile_gap:
	$(CC) -o gap gap.cpp hamiltonian.cpp check.cpp write_data.cpp operations.cpp parameters.cpp print.cpp linear_algebra.cpp -lgsl -llapack -lblas -std=gnu++11
