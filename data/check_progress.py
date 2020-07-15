import os, glob, re, numpy as np
from shutil import copyfile
dir_name = ["4x4/2_occupants/","4x4/3_occupants/"]

total_files, finished_files, unfinished_files, unfinished_taus, unfinished_dist = 0,0,[],[],[]

for dir_ in dir_name:
    for fname in glob.glob("unfinished_files/"+dir_ +"*.txt"):
        os.remove(fname)

for dir_ in dir_name:
    for fname in glob.glob(dir_ +"*.txt"):
        total_files += 1
        distance = 100
        with open(fname) as f:
            last_tau = 0
            for line in f:
                if "distance" in line:
                    distance = min(float(line.split(" ")[-1]), distance)
                if "tau" in line and "arr" not in line and "best" not in line:
                    last_tau = float(line.split(" ")[-1])

        if distance > 1 or distance < 0:
            total_files -=1
        elif distance < 0.1:
            finished_files += 1
        else:
            if not os.path.exists("unfinished_files/"+dir_): os.makedirs("unfinished_files/"+dir_)
            copyfile(fname, "unfinished_files/"+fname)
            unfinished_files.append(fname)
            unfinished_taus.append(last_tau)
            unfinished_dist.append(distance)

# for x in range(len(unfinished_files)):
#     print()
#     print(unfinished_files[x])
#     print(unfinished_taus[x])
#     print(unfinished_dist[x])

print(str(finished_files) + "/" + str(total_files) + " files finished")
average_dist = np.mean(np.array(unfinished_dist))
print("average distance for unfinished files: ", average_dist)




for dir_ in dir_name:
    unfinished_parameters = open("unfinished_files/"+dir_+"unfinished_parameters.txt", "w")
    unfinished_parameters.write("ji ki jt kt\n")

    for fname in glob.glob("unfinished_files/"+dir_ +"MC*.txt"):

        fname = fname.split("_ji=")[-1]
        ji, fname = fname.split("_ki=")[0],fname.split("_ki=")[-1]
        ki, fname = fname.split("_jt=")[0],fname.split("_jt=")[-1]
        jt, fname = fname.split("_kt=")[0],fname.split("_kt=")[-1]
        kt = fname.split(".txt")[0]
        unfinished_parameters.write(str(ji) + " " + str(ki) + " " + str(jt) + " " + str(kt) + "\n")
