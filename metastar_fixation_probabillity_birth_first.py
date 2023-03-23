def phi(r,t, N):
    """this function returns the fixation probability 
    of a well-mixed population
    {Parameters: 
    r: selection parameter for birth
    t: selection parameter for death
    N: population size}
    """
    if r==1 and t==1:
        return 1/N
    else:
        return (1-t/r)/(1-(t/r)**N)

def T_0_j(t1,r1,t,r,j,N,M):
    """
    This function returns \pi_{\circ}^{j-} of the manuscript
    """
    return t1/r1* phi(1/r,1/t, N)/(t1/r1*phi(1/r,1/t, N)+(M-1-j+t1*j)*phi(r,t,N))


def T_1_j(t1,r1,t,r,j, N,M):
    """
    this function returns \pi_{\bullet}^{j+} of the manuscript
    {Parameters: 
    r: selection parameter for birth at the individual level
    t: selection parameter for death at the individual level
    r1: selection parameter for birth at the patch level
    t1: selection parameter for death at the patch level
    j: number of mutant leaves
    N: number of individuals in each patc
    M: number of patches
    }
    """
    return r1*phi(r,t, N)/(r1*phi(r,t,N)+(M-1-j+t1*j)*phi(1/r,1/t,N))

def phi_1_0(t1,r1,t,r, N,M):
    """
    This function returns the fixation probability 
    of a metastar when we start with a mutant patch
    in the center
    """
    S=0
    for j in range(1,M-1):
        S = S + (T_0_j(t1,r1,t,r,j, N, M)/T_1_j(t1,r1,t,r,j, N, M))**j
    return T_1_j(t1,r1,t,r,0, N, M)/(1+(1-T_1_j(t1,r1,t,r,0, N, M))*S)

def phi_0_1(t1,r1,t,r, N, M):
    """
    This function returns the fixation probability
    of the metastar starting from a mutant leaf
    """
    return (1-T_0_j(t1,r1,t,r,1, N, M))*phi_1_0(t1,r1,t,r, N,M)/T_1_j(t1,r1,t,r,1, N, M)

def phi_patch(t1,r1,t,r, N, M):
    """
    This function returns the average fixation probability
    of the metastar starting with a mutant patch 
    """
    return ((M-1)*phi_0_1(t1,r1,t,r, N, M) + phi_1_0(t1,r1,t,r, N,M))/M



import numpy as np
import matplotlib.pyplot as plt

#parameters
M = 5
N = 5
R = np.arange(0.1,4.1,0.1)

#MBBDD : t>1, r>1, t1>1, r1>1
phi_star = []
for r in R:
    t=1/r
    t1=1/r
    r1 =r
    
    #print(phi_patch(t1,r1,t,r, N, M) * phi(r,t, N))
    phi_star.append(phi_patch(t1,r1,t,r, N, M) * phi(r,t, N))  

plt.plot(R, phi_star, label = 'MBBDD')    

#MBBdD: t>1, r>1, t1=1, r1>1
phi_star = []
for r in R:
    t1 = 1
    t = 1/r
    r1 =r
    
    phi_star.append(phi_patch(t1,r1,t,r, N, M) * phi(r,t, N)) 

plt.plot(R, phi_star, label = 'MBBdD')

#MbBDD
phi_star = []
for r in R:
    t1 = 1/r
    t = 1/r
    r1 =1
    #print(phi_patch(t1,r1,t,r, N, M) * phi(r,t, N))
    phi_star.append(phi_patch(t1,r1,t,r, N, M) * phi(r,t, N)) 

plt.plot(R, phi_star, label = 'MbBDD')

#MbBdD
phi_star = []
for r in R:
    t1 = 1
    t = 1/r
    r1 =1
    phi_star.append(phi_patch(t1,r1,t,r, N, M) * phi(r,t, N)) 

plt.plot(R, phi_star, label = 'MbBdD')

#MBBdd: t=1, r>1, t1=1, r1>1
phi_star = []
for r in R:
    t1 = 1
    t = 1
    r1 =r
    #print(phi_patch(t1,r1,t,r, N, M) * phi(r,t, N))
    phi_star.append(phi_patch(t1,r1,t,r, N, M) * phi(r,t, N))

plt.plot(R, phi_star, label = 'MBBdd')    

#MBBDd: t=1, r>1, t1>1, r1>1
phi_star = []
for r in R:
    t = 1
    t1 = 1/r
    r1 =r
    #print(phi_patch(t1,r1,t,r, N, M) * phi(r,t, N))
    phi_star.append(phi_patch(t1,r1,t,r, N, M) * phi(r,t, N)) 

plt.plot(R, phi_star, label = 'MBBDd') 

#MBbDD
phi_star = []
for r in R:
    t1 = 1/r
    t = 1/r
    r1 =r
    r=1
    #print(phi_patch(t1,r1,t,r, N, M) * phi(r,t, N))
    phi_star.append(phi_patch(t1,r1,t,r, N, M) * phi(r,t, N)) 

plt.plot(R, phi_star, label='MBbDD') 

#MBbdD
phi_star = []
for r in R:
    t1 = 1
    t = 1/r
    r1 =r
    r=1
    #print(phi_patch(t1,r1,t,r, N, M) * phi(r,t, N))
    phi_star.append(phi_patch(t1,r1,t,r, N, M) * phi(r,t, N)) 

plt.plot(R, phi_star, label='MBbdD')     

#MbBDd
phi_star = []
for r in R:
    t1 = 1/r
    t = 1
    r1 =1
    #print(phi_patch(t1,r1,t,r, N, M) * phi(r,t, N))
    phi_star.append(phi_patch(t1,r1,t,r, N, M) * phi(r,t, N)) 

plt.plot(R, phi_star, label='MbBDd') 

#MbBdd
phi_star = []
for r in R:
    t1 = 1
    t = 1
    r1 =1
    phi_star.append(phi_patch(t1,r1,t,r, N, M) * phi(r,t, N))

plt.plot(R, phi_star, label='MbBdd')

#MbbDD
phi_star = []
for r in R:
    t1 = 1/r
    t = 1/r
    r1 =1
    r=1
    phi_star.append(phi_patch(t1,r1,t,r, N, M) * phi(r,t, N))

plt.plot(R, phi_star, label='MbbDD')    

#MbbdD

phi_star = []
for r in R:
    t1 = 1
    t = 1/r
    r1 =1
    r=1
    phi_star.append(phi_patch(t1,r1,t,r, N, M) * phi(r,t, N)) 

plt.plot(R, phi_star, label='MbbdD')

#MBbDd
phi_star = []
for r in R:
    t1 = 1/r
    t = 1
    r1 =r
    r=1
    #print(phi_patch(t1,r1,t,r, N, M) * phi(r,t, N))
    phi_star.append(phi_patch(t1,r1,t,r, N, M) * phi(r,t, N)) 

plt.plot(R, phi_star, label='MBbDd')

#MBbdd
phi_star = []
for r in R:
    t1 = 1
    t = 1
    r1 =r
    r=1
    #print(phi_patch(t1,r1,t,r, N, M) * phi(r,t, N))
    phi_star.append(phi_patch(t1,r1,t,r, N, M)* phi(r,t, N)) 

plt.plot(R, phi_star, label='MBbdd')

#MbbDd
phi_star = []
for r in R:
    t1 = 1/r
    t = 1
    r1 =1
    r=1
    phi_star.append(phi_patch(t1,r1,t,r, N, M) * phi(r,t, N)) 
    
plt.plot(R, phi_star, label='MbbDd')

#Mbbdd

phi_star = []
for r in R:
    t1 = 1
    t = 1
    r1 =1
    r=1
    phi_star.append(phi_patch(t1,r1,t,r, N, M)* phi(r,t, N)) 

plt.plot(R, phi_star, label='Mbbdd')    

plt.legend()

plt.show()