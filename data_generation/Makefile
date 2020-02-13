CC=g++
no=2
gi=0.1
fi=0.1
gt=1.0
ft=1.0


compile:
	$(CC) -o main main.cpp adiabatic.cpp check.cpp export_data.cpp hamiltonian.cpp linear_algebra.cpp mcbb.cpp mcbf.cpp mcdb.cpp operations.cpp print.cpp parameters.cpp -lgsl -llapack -lblas -std=gnu++11

run:
	./main $(no) $(gi) $(fi) $(gt) $(ft)

submit:
	cd cluster	
	$(CC) -o make_submit make_submit.cpp -std=gnu++11
	./make_submit
	cd ..
	condor_submit submit_file.txt
