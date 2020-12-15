import glob, os



dir_ = "5x5/3_occupants/"

mc_average = 0
total_files = 0
tau_step_avg = 0
for file in glob.glob(dir_+"*txt"):

    best_mc, tau = 0,0
    tau_steps = 0
    for line in open(file, "r"):
        if (("tau =" in line) and ("_" not in line)): 
            tau_steps+=1
        if ("best_mc_result " in line and "arr" not in line): 
            best_mc = float(line.split("=")[-1])



    tau_step_avg += tau_steps
    mc_average += best_mc
    total_files += 1

mc_average = mc_average/total_files
average_tau_steps = tau_step_avg/total_files

print("best_mc_result average, average tau steps", mc_average, average_tau_steps)




