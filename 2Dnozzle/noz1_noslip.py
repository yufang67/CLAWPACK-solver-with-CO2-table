import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.colors as colors
#from matplotlib import gridspec

#with open('fort.q0007') as f:
#    f=[x.strip() for x in f if x.strip()]
#    datap=[tuple(map(float,x.split())) for x in f[2:]]
#line3 = f[2]
#line4 = f[3]
#nx = int(line3.split()[0])
#ny = int(line4.split()[0])        
#dx = 1.0/nx
#dy = 1.0/ny
#data = np.loadtxt('fort.q0007', skiprows=9)
#############################################
#x_comp = np.zeros((nx))
#y_comp = np.zeros((ny))
#for i in range(0,nx):
#        x_comp[i] = (i)*dx
#for j in range(0,ny):
#    y_comp[j] = (j)*dy
#xmap,ymap = np.meshgrid(x_comp,y_comp)
#rho_map = np.zeros((nx,ny))
#rho_map = np.reshape(rho,(nx,ny))
############################################
#f = np.loadtxt('1Dshocktube_vri/TEST.C0002')
data = np.loadtxt('NOZ1.C0040')
data2 = np.loadtxt('NOZ1.C0010')
data3 = np.loadtxt('naka09_noz1_937.tsv')
data4 = np.loadtxt('naka09_noz1_937_points.tsv')

x4sat = data4[3:11,0]
p4sat = data4[3:11,1]

x4m = data4[11:15,0]
p4m = data4[11:15,1]

x3 = data3[5:31,0]
p3 = data3[5:31,1]




nx = 504
ny = 44


xx  = data[:,8] 
yy  = data[:,9]
rho = data[:,1]
press  = data[:,0]
sound  = data[:,2]
mass   = data[:,10]
temp   = data[:,4]
eint   = data[:,5]
ux     = data[:,6]
uy     = data[:,7]
#
xx2  = data2[:,8] 
yy2  = data2[:,9]
press2  = data2[:,0]
rho2 = data2[:,1]
sound2  = data2[:,2]
ux2     = data2[:,6]
uy2     = data2[:,7]
temp2   = data2[:,4]
mass2   = data2[:,10]
#
r = np.zeros((ny,nx))
x = np.zeros((ny,nx))
y = np.zeros((ny,nx))
p = np.zeros((ny,nx))
c = np.zeros((ny,nx))
Y = np.zeros((ny,nx))
T = np.zeros((ny,nx))
e = np.zeros((ny,nx))
u = np.zeros((ny,nx))
v = np.zeros((ny,nx))
x_dis = []   #np.zeros((1,100))
y_dis = []   #np.zeros((1,100))
q_dis = []   #np.zeros((1,100))

x_dis2 = []   #np.zeros((1,100))
y_dis2 = []   #np.zeros((1,100))
Q_dis2 = []   #np.zeros((1,100))
#
bc = []   #np.zeros((1,100))
#
x = np.reshape(xx,(ny,nx))
y = np.reshape(yy,(ny,nx))
r = np.reshape(rho,(ny,nx))
p = np.reshape(press,(ny,nx))
c = np.reshape(sound,(ny,nx))
yyy = np.reshape(mass,(ny,nx))
T = np.reshape(temp,(ny,nx))
e = np.reshape(eint,(ny,nx))
u = np.reshape(ux,(ny,nx))
v = np.reshape(uy,(ny,nx))
# 
x2 = np.reshape(xx2,(ny,nx))
y2 = np.reshape(yy2,(ny,nx))  
p2 = np.reshape(press2,(ny,nx)) 
r2 = np.reshape(rho2,(ny,nx))
u2 = np.reshape(ux2,(ny,nx))
c2 = np.reshape(sound2,(ny,nx))
v2 = np.reshape(uy2,(ny,nx))
T2 = np.reshape(temp2,(ny,nx))
yyy2 = np.reshape(mass2,(ny,nx))
#
for i in range(nx):
    x_dis.append(x[ny/2,i]*1000-2)
#    Q_dis.append(yyy[ny/2,i])
    q_dis.append(p[ny/2,i]/1e6)
#    Q_dis.append(c[ny/2,i])
#    Q_dis.append((p[ny/2,i] + r[ny/2,i]*(u[ny/2,i]*u[ny/2,i]+v[ny/2,i]*v[ny/2,i])))
#    Q_dis.append(r[ny/2,i])
#    Q_dis.append(np.sqrt(u[ny/2,i]*u[ny/2,i]+v[ny/2,i]*v[ny/2,i])/c[ny/2,i])
#    Q_dis.append(u[ny/2,i])
#    q_dis.append(T[ny/2,i])
for i in range(nx):
    x_dis2.append(x2[ny/2,i]*1000-2)
#    y_dis2.append(y2[ny/2,i])
    Q_dis2.append(p2[ny/2,i]/1e6)
#    Q_dis2.append(r2[ny/2,i])
#    Q_dis2.append(T2[ny/2,i])
#    Q_dis2.append(yyy2[ny/2,i])
#    Q_dis.append(np.sqrt(u2[ny/2,i]*u2[ny/2,i]+v2[ny/2,i]*v2[ny/2,i])/c2[ny/2,i])  
for j in range(ny):
    y_dis.append(y[j,nx-4]*1000)
#    bc.append(np.sqrt(u[j,nx-3]*u[j,nx-3]+v[j,nx-3]*v[j,nx-3])/c[j,nx-3])
#    bc.append(c[j,nx-4])
    bc.append(u[j,nx-4])
#    bc.append(v[j,nx-4])





####  MAP ##############
plt.figure(figsize=(8,4))
ax = plt.gca()
#
#plt.pcolormesh(x*1000, y*1000, r, cmap='jet',shading='gouraud')

#plt.pcolormesh(x, y, p, cmap='jet',shading='gouraud')
#plt.pcolormesh(x*1000, y*1000, r, cmap='jet',shading='gouraud')
#ax.text(35,5,'Mpa',fontsize=11,fontweight='bold')
#plt.pcolormesh(x*1000, y*1000, c, cmap='jet',shading='gouraud')
#plt.pcolormesh(x*1000, y*1000, Y, cmap='jet',shading='gouraud')
#ax.text(35,5,'Quality',fontsize=11,fontweight='bold')
#plt.pcolormesh(x*1000, y*1000, T, cmap='jet',shading='gouraud')
#plt.pcolormesh(x*1000, y*1000, e/1e3, cmap='jet',shading='gouraud')
plt.pcolormesh(x*1000, y*1000, u, cmap='jet',shading='gouraud')
#plt.pcolormesh(x*1000, y*1000, v, cmap='jet',shading='gouraud')
#
plt.xticks(size = 12)
plt.yticks(size = 12)
#ax.set_xlim([-2, 29])
#ax.set_ylim([-2, 11])
#ax.tick_params(labelsize=5)
#plt.grid(True)
plt.colorbar(fraction=0.046,pad=0.01)
#plt.ylabel('y $[mm]$',fontsize=15,position=(0,0.5),rotation = "vertical")
#plt.xlabel('x $[mm]$',fontsize=15,rotation = "horizontal")
#plt.savefig("map_xvelo.pdf")
plt.figure(figsize=(6,4))
plt.plot(bc,y_dis,linewidth = 2,linestyle='-',color = 'red')
#plt.axis([35,95,4.5,5.5])
##
plt.figure(figsize=(8,5))
plt.plot(x_dis,q_dis,linewidth = 3,linestyle=':',color = 'blue',label='CLAWPACK_slip')
plt.scatter(x4m,p4m,linewidth = 2, marker = 'd',color = 'red',label='strain gauge Nakagawa et al.2009')
plt.plot(x4sat,p4sat, 'ks--',markerfacecolor='white',label='$P_{sat}(T_{wall})$ Nakagawa et al.2009')
plt.plot(x3,p3,linewidth = 3, marker = 's',color = 'orange')
plt.plot(x_dis2,Q_dis2,linewidth = 2,linestyle='-',color = 'black',label='CLAWPACK_noslip')
#plt.xticks(np.arange(0, 40, 2))
plt.axvline(x=27.35,color='grey',linestyle='--') 
plt.xlabel('x $[mm]$',fontsize=12,rotation = "horizontal") 
plt.ylabel('Pressure $[Mpa]$',fontsize=12,position=(0,0.5),rotation = "vertical")
plt.axis([-5,84,2,8])
#ax.locator_params(numticks=12)
plt.legend()
plt.grid(True)
#plt.savefig("press.pdf")
plt.tight_layout()

#plt.figure(figsize=(8,4))
#plt.scatter(xx,yy,s=0.1)
#plt.grid(True)
#plt.savefig("grid.pdf")

plt.show()
