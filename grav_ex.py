# Import necessary libraries/routines

import numpy as np
import sys as sys
import matplotlib
import matplotlib.pyplot as plt
plt.rc('font', family='serif')
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
from mpl_toolkits.mplot3d import Axes3D
import time

#####################################################################
# analytical solutions for exercise 1 (or 3) and 2
#####################################################################

def analytical_solution(x,y,z):
    if ex==1 or ex==3:
        gx=0.
        gy=0.
        gz=Ggrav*volth*rho0/z**2
        U=-Ggrav*volth*rho0/z
    if ex==2:
        r2=x**2+y**2+z**2
        r=np.sqrt(r2)
        if r2<R**2/4.:
            gx=0.
            gy=0.
            gz=0.
            U=2*np.pi*Ggrav*rho0*(R**2/4.-R**2)  #this
        elif r2<R**2:
            gx=4*np.pi/3.*Ggrav*rho0*(r-R**3/8./r2) #this
            gy=0.
            gz=0.
            U=4*np.pi/3.*Ggrav*rho0*(r2/2+R**3/8./r)-2*np.pi*Ggrav*rho0*R**2  #
        else:
            gx=Ggrav*(4.*np.pi/3.*rho0*(R**3-R**3/8.))/r2 #
            gy=0.
            gz=0.
            U=-Ggrav*(4.*np.pi/3.*rho0*(R**3-R**3/8.))/r  #
        #end if
    return gx,gy,gz,U

#####################################################################
# PREM density
#####################################################################

def prem_density(radius):
    x=radius/6371.0e3
    if radius>6371.0e3:
        densprem=0
    elif radius<=1221.5:
        densprem=13.0885-8.8381*x**2
    elif radius<=3480.0:
        densprem=12.5815-1.2638*x-3.6426*x**2-5.5281*x**3
    elif radius<=3630.0e3:
        densprem=7.9565-6.4761*x+5.5283*x**2-3.0807*x**3
    elif radius<=5600.0e3:
        densprem=7.9565-6.4761*x+5.5283*x**2-3.0807*x**3
    elif radius<=5701.0e3:
        densprem=7.9565-6.4761*x+5.5283*x**2-3.0807*x**3
    elif radius<=5771.0e3:
        densprem=5.3197-1.4836*x
    elif radius<=5971.0e3:
        densprem=11.2494-8.0298*x
    elif radius<=6151.0e3:
        densprem=7.1089-3.8045*x
    elif radius<=6291.0e3:
        densprem=2.6910+0.6924*x
    elif radius<=6346.0e3:
        densprem=2.6910+0.6924*x
    elif radius<=6356.0e3:
        densprem=2.9
    elif radius<=6368.0e3:
        densprem=2.6
    else:
        densprem=1.020
    return densprem*1000

#####################################################################
# plot gravity potential and vector components
#####################################################################

def plotting_figs(distance,g,gth,U,Uth,ex,N):
    plt.figure(figsize=(20,10))
    plt.rc('font', size=20)
    
    if ex==1 or ex==3:
        legend_g='g$_{z}$'
        legend_gth='g$_{zth}$'
        x_label='z (m)'
        y_label='g$_{z}$ (ms$^{-2}$)'
    elif ex==2:
        legend_g='g$_{x}$'
        legend_gth='g$_{xth}$' 
        x_label='x (m)'
        y_label='g$_{x}$ (ms$^{-2}$)'
    # plot gravity vectors
    plt.subplot(1,2,1)
    plt.plot(distance,g, label=legend_g, lw=2.0, marker='o', c='black')
    plt.plot(distance,gth, label=legend_gth, lw=2.0, ls='--', marker='o', c='cyan')
    plt.legend()
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if ex==1 or ex==3:
        plt.xscale('log')
    elif ex==2:
        pass
    plt.grid(which='both')

    # plot gravity potential
    plt.subplot(1,2,2)
    plt.plot(distance,U, label='U', lw=2.0, marker='o', c='green')
    plt.plot(distance,Uth, label='U$_{th}$', lw=2.0, ls='--', marker='o', c='orange')
    plt.legend()
    plt.xlabel(x_label)
    plt.ylabel('U (Jkg$^{-1}$)')
    if ex==1 or ex==3:
        plt.xscale('log')
    elif ex==2:
        pass
    plt.grid(which='both')
    plt.savefig('./gravity_plot_ex%s_n%s.png'%(ex,N), dpi=200) # save figure
    plt.close()

#####################################################################
# convert spherical to cartesian and vice-versa [ref.1]
#####################################################################
def Cartesian2Spherical(p):
    # p = (x,y,z)
    # theta  in (0,pi) and phi in (0,2pi)
    x,y,z = p
    r = np.sqrt(x*x+y*y+z*z)
    theta = np.arccos(z/r)    # lattitude (inclination)
    phi = np.arctan2(y,x)     # longitude (azimuth)
    q = np.array([r,theta,phi])
    return q

def Spherical2Cartesian(q):
    # q = (r,theta,phi)
    # theta  in (0,pi) and phi in (0,2pi)
    r,theta,phi = q
    SinTheta = np.sin(theta)
    CosTheta = np.cos(theta)
    SinPhi = np.sin(phi)
    CosPhi = np.cos(phi)
    rSinTheta = r*SinTheta
    x = rSinTheta*CosPhi
    y = rSinTheta*SinPhi
    z = r*CosTheta
    p  = np.array([x,y,z])
    return p

#####################################################################
# initial set-up
#####################################################################

# choose a value for N
N=10

ex=1  # exercise 1 or 2

R=6371.0e3
rho0=3000.
Ggrav=6.6738480e-11

h=2.*R/N
dV=h**3                                #compute dV
NP=N**3

if ex==1 or ex==3:
    volth=4./3.*np.pi*R**3             #analytical volume of solid sphere
    m=9

if ex==2:
    volth=4./3.*np.pi*(R**3-R**3/8.0)  #analytical volume of hollow sphere
    m=100

#####################################################################
# point layout
#####################################################################

if ex==1 or ex==2:                                  
    x = np.empty(NP,dtype=np.float64) 
    y = np.empty(NP,dtype=np.float64)
    z = np.empty(NP,dtype=np.float64)
    rho = np.empty(NP,dtype=np.float64)

    counter=0
    for i in range(0,N):
        for j in range(0,N):
            for k in range(0,N):
                #implement code to compute x-coordinate >>> x[counter]=
                #implement code to compute y-coordinate >>> y[counter]=
                #implement code to compute z-coordinate >>> z[counter]=
                counter+=1
            #end for
        #end for
    #end for
elif ex==3:                                         
    h1=R/N
    h2=np.pi/N
    h3=2*np.pi/N
    r = np.empty(NP,dtype=np.float64) 
    theta = np.empty(NP,dtype=np.float64)
    phi = np.empty(NP,dtype=np.float64)
    rho = np.empty(NP,dtype=np.float64)
    counter=0
    for i in range(0,N):
        for j in range(0,N):
            for k in range(0,N):
                #implement code to compute x-coordinate >>> r[counter]=
                #implement code to compute y-coordinate >>> theta[counter]=
                #implement code to compute z-coordinate >>> phi[counter]=           
                counter+=1
    x,y,z=Spherical2Cartesian([r,theta,phi])   #convert r,theta,phi to x,y,z

    # plot clouds of points and save the figure               
    fig = plt.figure(figsize=(13,10))
    plt.rc('font', size=20)
    ax = Axes3D(fig)
    pp = ax.scatter(x, y, z, c=x, s=40, cmap='winter', alpha=1, edgecolor='k')
    ax.tick_params(axis='X', pad=15)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.xaxis.labelpad= 25
    ax.yaxis.labelpad= 25
    ax.zaxis.labelpad= 25
    clb = fig.colorbar(pp,shrink=0.55, pad=0.05)
    clb.ax.tick_params(labelsize=20)
    clb.ax.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
    clb.formatter.set_powerlimits((0, 0))
    clb.ax.yaxis.set_offset_position('left')  
    clb.update_ticks()
    clb.ax.set_title('x', pad=30)
    ax.view_init(30,295)
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    ax.ticklabel_format(style='sci', axis='z', scilimits=(0,0))
    plt.savefig("./plot_point_clouds_ex%s_n%s.png"%(ex,N), dpi=200)
    plt.close()

#####################################################################
# setup density 
#####################################################################

if ex==1 or ex==3:
    for i in range(0,NP):
        r2=x[i]**2+y[i]**2+z[i]**2
        if r2<R**2:
            #implement density acoording to the question >>> rho[i]=
            #and remove the "pass" the following pass line
            pass
            #density according to PREM
            #rho[i]=prem_density(np.sqrt(r2))
        else:
            #implement density acoording to the question >>> rho[i]=
            #and remove the "pass" the following pass line
            pass
        #end if
    #end for

elif ex==2:
    for i in range(0,NP):
        r2=x[i]**2+y[i]**2+z[i]**2
        if r2<R**2 and r2>R**2/4.:
            #implement density acoording to the question >>> rho[i]= 
            #and remove the "pass" the following pass line
            pass
        else:
            #implement density acoording to the question >>> rho[i]=
            #and remove the "pass" the following pass line
            pass
        #end if
    #end for

#####################################################################
# save the mesh information
#####################################################################

np.savetxt('mesh_ex%s_n%s.ascii'%(ex,N),np.array([x,y,z,rho]).T,header='# x,y,z,rho')

#####################################################################
# compute numerical volume
#####################################################################

if ex==1 or ex==2:
    vol=0.
    mass=0.
    for i in range(0,NP):
        if rho[i]>0:
            vol+=dV
            mass+=rho[i]*dV

    print('N=',N,'vol=',vol,'vol. rel. err.=',(vol-volth)/volth, 'mass=',mass)

elif ex==3:
    r_e = np.linspace(0,R,N+1); r_e_center = [i+h1/2. for i in r_e]
    theta_e = np.linspace(0,np.pi,N+1); theta_e_center = [i+h2/2. for i in theta_e]
    phi_e = np.linspace(-np.pi,np.pi,N+1); phi_e_center = [i+h3/2. for i in phi_e]
    vol = 0.
    for r_ in r_e_center:
        for ang_ in theta_e_center:
            for azm_ in phi_e_center:
                if (r_<R) and (ang_<np.pi) and (azm_<np.pi):
                    vol += (r_)**2*np.sin(ang_)*h1*h2*h3  

    print('N=',N,'vol=',vol,'vol. rel. err.=',(vol-volth)/volth)

#####################################################################
# compute moment of inertia 
#####################################################################

I=0.
for i in range(0,NP):
    #implement moment of inertia acoording to the question >>> I+=
    #and remove the "pass" the following pass line
    pass
#end for
I*=8.*np.pi/3.

print('I=',I)

#####################################################################
# setup measurement points
#####################################################################

xm = np.empty(m,dtype=np.float64) 
ym = np.empty(m,dtype=np.float64)
zm = np.empty(m,dtype=np.float64)

if ex==1 or ex==3:
    for k in range(0,m):
        xm[k]=0.
        ym[k]=0.
        zm[k]=R+10**k

if ex==2:
    dx=3.0*R/m
    for k in range(0,m):
        xm[k]=k*dx+dx/2
        ym[k]=0.
        zm[k]=0.

#####################################################################
# compute gravity fields
#####################################################################

gx = np.zeros(m,dtype=np.float64) 
gy = np.zeros(m,dtype=np.float64)
gz = np.zeros(m,dtype=np.float64)
U = np.zeros(m,dtype=np.float64)
gxth = np.zeros(m,dtype=np.float64) 
gyth = np.zeros(m,dtype=np.float64)
gzth = np.zeros(m,dtype=np.float64)
Uth = np.zeros(m,dtype=np.float64)

for k in range(0,m):
    for i in range(0,NP):
        #implement code to compute numerical gravity potential and vector components i.e. gx,gy,gz,U
        #<<<code here>>>
        #and remove the "pass" the following pass line
        pass
    #end for
    gxth[k],gyth[k],gzth[k],Uth[k]=analytical_solution(xm[k],ym[k],zm[k])
#end for

#####################################################################
# save gravity potential and vector components
#####################################################################

np.savetxt('grav_ex%s_n%s.ascii'%(ex,N),np.array([xm,ym,zm,gx,gy,gz,U,gxth,gyth,gzth,Uth]).T,\
               header='# xm,ym,zm,gx,gy,gz,U,gxth,gyth,gzth,Uth')

#####################################################################
# plotting gravity potential and vector components
#####################################################################

# now call plotting function
if ex==1 or ex==3:
    plotting_figs(zm,gz,gzth,U,Uth,ex,N)
elif ex==2:
    plotting_figs(xm,gx,gxth,U,Uth,ex,N)

