
# PROGRAM FOR BEAM ANALYSIS IN FEA 
# To find deflection and slope for beams of rectangular cross section.

# import numpy 
# import pandas 
import numpy as np
import pandas as pd


# Input from user like Number of elements in beam , Cross section area, Young's modulus, Area moment of inertia .
Elements=int(input("Enter number of elements in beam: "))
nodes=Elements+1
print("Total number of nodes: ",nodes)
gsm=np.zeros((2*nodes,2*nodes), dtype=float)
dands=np.zeros((2*nodes,1))
sn=[]
en=[]
le=[]
F=np.zeros((2*nodes,1),dtype=float)
Fixed=[]
Hingroll=[]
print("\nEnter \'i\' if beam has uniform rectangular cross-section throughtout length and same material.")
print("Enter \'d\' if beam has non-uniform rectangular cross-sections or different material.")
beamcond=input("Enter cross section type : ")
if beamcond in ['I','i']:
    e=float(input("Enter Young's modulus in kN/m^2 : "))
    b = float(input("Enter width of cross-section in m : "))
    h = float(input("Enter height of cross-section in m : "))
    i=(b*(h**3))/12
    for j in range(Elements):
        l=float(input("Enter length of element %d : "%(j+1)))
        le.append(l)
        C=(e*i)/np.power(l,3)
        mat=np.array([[12,6*l,-12,6*l],
                      [6*l,4*l**2,-6*l,2*l**2],
                      [-12,-6*l,12,-6*l],
                      [6*l,2*l**2,-6*l,4*l**2]])
        elesm=np.multiply(C,mat)
        gsm[j*2:j*2+4,j*2:j*2+4]+=elesm
if beamcond in ['D','d']:
    for k in range(Elements):
        b = float(input("Enter width of cross-section in m or in mm for element %d : "%(k+1)))
        h = float(input("Enter height of cross-section in m or in mm for element %d : "%(k+1)))
        e = float(input("Enter Young's modulus in kN/m^2 or in MPa for element %d : "%(k+1)))
        i = (b*(h**3))/12
        l = float(input("Enter length of element %d : " % (k+1)))
        le.append(l)
        C = (e*i)/np.power(l, 3)
        mat = np.array([[12, 6*l, -12, 6*l],
                        [6*l, 4*l**2, -6*l, 2*l**2],
                        [-12, -6*l, 12, -6*l],
                        [6*l, 2*l**2, -6*l, 4*l**2]])
        elesm = np.multiply(C, mat)
        gsm[k*2:k*2+4, k*2:k*2+4] += elesm
print("\n")


# Formation  of Global stiffness matrix.
print("\n")
print("--------------------------------GLOBAL STIFFNESS MATRIX(10^6*K)--------------------------------")
print("\n K= \n")
K=pd.DataFrame(data=gsm,index=range(1,2*nodes+1),columns=range(1,2*nodes+1))
print(K)
print("\n")


# Formation Global Load vector.
nol = int(input("Enter number of loads acting on beam : "))
for n in range(nol):
    print("For load no.%d" % (n+1))
    print("Enter \'p\' for point load.")
    print("Enter \'udl\' for uniformly distributed load.")
    print("Enter \'uvl\' for unformly varying load.")
    loadtype=input("Enter type of load acting : ")
    if loadtype in ['p','P']:
        nn=int(input("Enter node number on which point load is acting: "))
        val=float(input("Enter load value : " ))
        F[2*nn-2,0]+=val
    if loadtype in ['udl','Udl','UDL']:
        snn=int(input("Enter start node number for udl : "))
        enn=int(input("Enter end node number for udl : "))
        val=float(input("Enter load intensity : " ))
        f=(val*le[snn-1])/2
        m=(val*(le[snn-1])**2)/12
        F[2*snn-2,0]+=f
        F[(2*snn+1)-2,0]+=m
        F[2*enn-2,0]+=f 
        F[(2*enn+1)-2,0]+=(-m)*1
    if loadtype in ['Uvl','uvl','UVL']:
        snn=int(input("Enter start node number for uvl : "))
        enn=int(input("Enter end node number for uvl : "))
        val=float(input("Enter maximum load intensity for uvl : "))
        f=(3*val*le[snn-1])/20
        f1=(7*val*le[snn-1])/20 
        m=(val*(le[snn-1])**2)/30
        m1=(val*(le[snn-1]**2))/20
        print("\n")
        print("Enter 'd' for decreasing intensity. ")
        print("Enter 'i' for increasing intensity. ")
        uvltyp=input("Enter type of uvl : ")
        if uvltyp in ['d','D']:
            F[2*snn-2,0]+=f1
            F[(2*snn+1)-2,0]+=m1
            F[2*enn-2,0]+=f
            F[(2*enn+1)-2,0]+=(-m)*1
        if uvltyp in ['i','I']:
            F[2*snn-2,0]+=f
            F[(2*snn+1)-2,0]+=m
            F[2*enn-2,0]+=f1
            F[(2*enn+1)-2,0]+=(-m1)*1


# Reduction of Global Stiffness Matrix and Global Load Vector
force=pd.DataFrame(F,index=range(2*nodes),columns=range(1))
KD=pd.DataFrame(gsm,index=range(2*nodes),columns=range(2*nodes))
nnws = int(input("Enter number of nodes having supports : "))
for m in range(nnws):
    nn = int(input("Enter node number having support: "))
    print("\n")
    print("Enter \'f\' for fixed support.")
    print("Enter \'h\' for hinge support.")
    print("Enter \'r\' for roller support.")
    scond = input("Enter type of support at node %d: " % nn)
    print("\n")
    if scond in ['f','F']:
        Fixed.append(nn)
        K1=KD.drop([nn*2-2,(nn*2+1)-2],axis=0)
        K1=K1.drop([nn*2-2,(nn*2+1)-2],axis=1)
        KD=K1
        f=force.drop(2*nn-2,axis=0)
        f=f.drop((2*nn+1)-2,axis=0)
        force=f
    elif scond in['h','H']:
        Hingroll.append(nn)
        K1=KD.drop(nn*2-2,axis=0)
        K1=K1.drop(nn*2-2,axis=1)
        KD=K1
        f=force.drop(2*nn-2,axis=0)
        force=f
    elif scond in ['r','R']:
        Hingroll.append(nn)
        K1=KD.drop(nn*2-2,axis=0)
        K1=K1.drop(nn*2-2,axis=1)
        KD=K1
        f = force.drop(2*nn-2, axis=0)
        force=f
redsmat=KD.to_numpy()
redfmat=force.to_numpy()
X=np.linalg.solve(redsmat,redfmat)


# To find values of Deflections and Slopes at Nodes  
print("\n-----------------------------------Nodal Deflections and Slopes (*10^-6) in m and rad------------------------------------------\n")
x=np.ones((2*nodes))
for y in Fixed:
    x[2*y-2]=0
    x[(2*y+1)-2]=0
for z in Hingroll:
    x[2*z-2]=0
u=0
for w in range(2*nodes):
    if x[w]==1:
        x[w]=X[u]
        u+=1
disround = np.round(x*10**6, 6)
a=disround.tolist()
deflection = a[0:2*nodes:2]
slope = a[1:2*nodes:2]
data = {"Nodes": range(1, nodes+1), "deflection": deflection, "slope": slope}
u1 = pd.DataFrame(data, columns=['Nodes', 'deflection', 'slope'])
print(u1.to_string(index=False))

# To find reactions at supports.
print("\n-------------------------------------Reactions at supports (kN)-------------------------------------\n")
vr=[]
mr=[]
for o in Fixed:
    verti=np.matmul(gsm[2*o-2],x)-F[2*o-2]
    moment=np.matmul(gsm[(2*o+1)-2],x)-F[(2*o+1)-2]
    vr.append(verti[0])
    mr.append(moment[0])
for p in Hingroll:
    verti=np.matmul(gsm[2*p-2],x)-F[2*p-2]
    vr.append(verti[0])
    mr.append(0)
N=Fixed+Hingroll
reactions={'Nodes':N,'Ry': vr,'Mz':mr}
R=pd.DataFrame(reactions,columns=['Nodes','Ry','Mz'])
print(R.to_string(index=False))


  