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


NX, num = 5,2
num_ri, num_rt = 50,50
ri_min, ri_max, rt_min, rt_max = 0.5,3,-3,1
total_time = 1.0
ris, rts = np.linspace(ri_min, ri_max, num_ri),np.linspace(rt_min, rt_max, num_rt)

filepath = 'submit_file_bifurcations.txt'

if os.path.exists(filepath): os.remove(filepath)
f = open(filepath, "a")

for ri in ris:
    for rt in rts:
        ji,ki = get_j_k(ri)
        jt,kt = get_j_k(rt)
        f.write("Executable       = mcdb_bifurcation_"+str(NX)+"x"+str(NX)+"\n")
        f.write("Arguments        = "+str(num)+" "+str(ji)+" "+str(ki)+" "+str(jt)+" "+str(kt)+" "+str(total_time)+"\n")
        f.write("queue\n")
f.close()
