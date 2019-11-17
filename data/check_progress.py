import glob
def main():
    seed_num = 10

    complete = []
    for x in range(seed_num):
        complete.append(0)
    total_files = 0
    for filepath in glob.iglob('3x3/4_occupants/*.txt'):
        for line in open(filepath):
            for x in range(2,seed_num+1):
                if "seed" in line and str(x) in line and "TOTAL" not in line:
                    complete[x-1] +=1
        total_files +=1
        
    for x in range(1,seed_num):
        print(str(complete[x]) + "/"+ str(total_files) + " files have finished seed" + str(x+1))

main()
