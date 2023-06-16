import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
from astropy.table import Table
import numpy as np
import sys
from matplotlib import cm
from matplotlib.colors import ListedColormap
from mpl_toolkits import mplot3d
from operator import itemgetter

#Checking if the program was given the right arguments

if len(sys.argv) != 4 :
    sys.exit("Usage: <InputText> <ChScProbability> <AngularDependence>")
    
InputText = sys.argv[1]
ChScProbability = sys.argv[2]
AngularDependence = sys.argv[3]


###############################################    MAIN   ######################################################
 
Names = np.array(["theta","phi","closest_PMT","Start_Time","Arrival_Time","type","NEvent"])
datas = Table.read(InputText,format="ascii",names=Names)

thetas = np.array(datas["theta"])
phis = np.array(datas["phi"])
closest_PMTs = np.array(datas["closest_PMT"])
Start_Times = np.array(datas["Start_Time"])
Arrival_Times = np.array(datas["Arrival_Time"])
types = np.array(datas["type"])
NEvent = np.array(datas["NEvent"])


N = len(thetas)

def cos_alpha (theta,phi) : #assuming the sun to be in position (0,0)
    return -np.sin(phi-np.pi/2)

NEvents = max(NEvent+1)

First_Hits = np.zeros(10)

iout = 0
First_thetas = []
First_phis = []

for j in range(0,NEvents) :
    
    Sliced_Vector = []
    Ordered_Vector = []

    for i in range (iout,N) :
        if NEvent[i] == j :
            Sliced_Vector.append([Start_Times[i],types[i],thetas[i],phis[i]])
        if NEvent[i] > j :
            iout = i
            break
        
    Ordered_Vector = sorted(Sliced_Vector, key=itemgetter(0), reverse=False)
    
    #print(Ordered_Vector)

    for k in range (0,min(10,len(Ordered_Vector))) :
        First_Hits[k] += Ordered_Vector[k][1]
        if k<5 :
            First_thetas.append(Ordered_Vector[k][2])
            First_phis.append(Ordered_Vector[k][3])

NBin = np.arange(0,10,1)

plt.hist(NBin,10,weights=First_Hits/NEvents)
plt.xlabel("# hit")
plt.ylabel("Cherenkov probability")
plt.title("Cherenkov probability")
plt.savefig(ChScProbability)
print(ChScProbability,"printed")


First_phis = np.array(First_phis)
First_thetas = np.array(First_thetas)

plt.hist(cos_alpha(First_thetas,First_phis),bins=20)
plt.xlabel("cos(a)")
plt.ylabel("# counts")
plt.title("Angular dependence in the first 5 hits")
plt.xlim(-1.1,1.1)
plt.savefig(AngularDependence)
print(AngularDependence,"printed")