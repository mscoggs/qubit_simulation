import os, glob, re, numpy as np

NX = "4",occ = "2";

submit_file = open("submit_file.txt", "w")
params = pd.read_csv("../../data/unfinished_files/"+NX+"x"+NX+"/"+occ+"_occupants/unfinished_parameters.txt")

for ji,ki,jt,kt in zip(params["ji"],params["ki"],params["jt"],params["kt"]):

	identifier = "_num_"+num+"__ji_"+ji+"__ki_"+ki+"__jt_"+jt+"__kt_"+kt
	submit_file.write("Executable        = main\n")
	submit_file.write("Arguments         = " << num << " " << ji << " " << ki << " " << jt << " " << kt << "\n")
	submit_file.write("Log               = cluster/logs/_" + identifier  + ".log\n")
	submit_file.write("Output            = cluster/outputs/_" << identifier << "\n")
	submit_file.write("request_cpus      = 1\n")
	#submit_file << "request_memory    = 20 GB\n"
	submit_file.write("queue\n")
