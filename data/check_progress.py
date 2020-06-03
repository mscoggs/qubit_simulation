import os, glob, re, numpy as np

dir_name = ["4x4/3_occupants/"]

total_files, finished_files, unfinished_files, unfinished_taus, unfinished_dist = 0,0,[],[],[]

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
        elif distance < 0.05:
            finished_files += 1
        else: 
            unfinished_files.append(fname)
            unfinished_taus.append(last_tau)
            unfinished_dist.append(distance)

for x in range(len(unfinished_files)):
    print()
    print(unfinished_files[x])
    print(unfinished_taus[x])
    print(unfinished_dist[x])

print(str(finished_files) + "/" + str(total_files) + " files finished")
average_dist = np.mean(np.array(unfinished_dist))
print("average distance for unfinished files: ", average_dist)


