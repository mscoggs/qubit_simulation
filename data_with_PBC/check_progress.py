import os, glob, re

dir_name = ["3x3/3_occupants/", "3x3/4_occupants/"]
#dir_name = ["3x3/3_occupants/"]
total_files = 0
finished_files = 0
for dir_ in dir_name:
    for fname in glob.glob(dir_ + "*.txt"):
        total_files += 1
        distance = 100
        with open(fname) as f:
            for line in f:
                if "distance" in line:
                    distance = min(float(line.split(" ")[-1]), distance)

        print(fname)
        print("distance: " + str(distance) + " \n\n")
        if distance < 0.05:
            finished_files += 1
        if distance > 1 or distance < 0:
            total_files -=1

print(str(finished_files) + "/" + str(total_files) + " files finished")



