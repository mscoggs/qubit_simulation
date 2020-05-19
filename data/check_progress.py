import os, glob, re

dir_name = ["3x3/3_occupants/", "3x3/4_occupants/"]
#dir_name = ["3x3/3_occupants/"]
total_files = 0
finished_files = 0
unfinished_files = []
unfinished_taus = []
for dir_ in dir_name:
    for fname in glob.glob(dir_ + "*.txt"):
        total_files += 1
        distance = 100
        with open(fname) as f:
            last_tau = 0
            for line in f:
                if "distance" in line:
                    distance = min(float(line.split(" ")[-1]), distance)
                if "tau" in line and "arr" not in line and "best" not in line:
                    last_tau = float(line.split(" ")[-1])
            
            unfinished_taus.append(last_tau)
        print(fname)
        print("distance: " + str(distance) + " \n\n")
        if distance < 0.05:
            finished_files += 1
        else: unfinished_files.append(fname)
        if distance > 1 or distance < 0:
            total_files -=1
print("\n\nUNFINISHED FILES:\n\n")
for x in range(len(unfinished_files)):
    print()
    print(unfinished_files[x])
    print(unfinished_taus[x])
print(str(finished_files) + "/" + str(total_files) + " files finished")



