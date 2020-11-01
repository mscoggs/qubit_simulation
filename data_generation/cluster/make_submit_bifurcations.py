import math
import numpy as np
import os

def get_j_k(r):
    r = math.exp(r)
    k = 0.5
    j = r*k
    while(j >= 1):
        k =k/2
        j =j/2
    return truncate(j,3),truncate(k,3)


def truncate(number, digits):
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper


NX, num = 3,3
num_ri, num_rt = 2,70
ri_min, ri_max, rt_min, rt_max = -1.3,-1.4,1,3
ris, rts = np.linspace(ri_min, ri_max, num_ri),np.linspace(rt_min, rt_max, num_rt)

filepath = 'submit_file_bifurcations.txt'

if os.path.exists(filepath): os.remove(filepath)
f = open(filepath, "a")

#for total_time in [0.5,0.7,0.9,1.1,1.3,1.4,1.5,1.6,1.7]:
for total_time in np.linspace(0.5,0.9,15):
    for ri in ris:
        for rt in rts:
            ji,ki = get_j_k(ri)
            jt,kt = get_j_k(rt)
            f.write("Executable       = mcdb_bifurcation_"+str(NX)+"x"+str(NX)+"\n")
            f.write("Arguments        = "+str(num)+" "+str(ji)+" "+str(ki)+" "+str(jt)+" "+str(kt)+" "+str(total_time)+"\n")
            f.write("queue\n")
f.close()
