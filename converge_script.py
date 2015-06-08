from sys import argv
import numpy as np
import os.path
print argv
script, res, f_phi_exact, f_phi_comp = argv 

#Read in the third column from f_phi_exact and f_phi_comp
f = open(f_phi_exact)
line = f.readline()
#print line
exact_values =[]
i = 0
res=int(res)
while (i < (res*res)):
    if (line != '\n'):
        tokens = line.split('\t')
        #print tokens
        value = float(tokens[2])
        exact_values.append(value)
        line = f.readline()
        i = i + 1
    else:
        line=f.readline()
#print exact_values

g = open(f_phi_comp)
line = g.readline()
#print line
comp_values = []
i = 0
while (i < (res*res)):
    if (line != '\n'):
        tokens = line.split('\t')
        #print tokens
        comp_value = float(tokens[2])
        comp_values.append(comp_value)
        line = g.readline()
        i = i + 1
    else:
        line=g.readline()
#print comp_values

#Substract the square of one col from the square of the other
sub = []
elem = len(comp_values)
for x in xrange(elem):
    sub.append(comp_values[x]**2-exact_values[x]**2)
sum_total = 0
for x in xrange(elem):
    sum_total += sub[x]
error_norm = np.sqrt(abs(sum_total)*(1/res**2))
print error_norm 
if not(os.path.isfile("/home/rocio/Thesis/Marker-Level-Set-2D/data/convergence.dat")):
    output_file = open("/home/rocio/Thesis/Marker-Level-Set-2D/data/convergence.dat","w")
else:
    output_file = open("/home/rocio/Thesis/Marker-Level-Set-2D/data/convergence.dat","a")
output_file.write("%f \t %f\n" %(res,error_norm))
#Add up all the values in the resulting array
#Divide by the res^2

#open the output file
#append to the end of the file the resolution and the result of the error norm
#close the file
#end.

