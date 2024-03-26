import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import math

G = 1.
M1 = 1.
#number of frames
sim_duration = 2000
# Computational timestep
dt = 0.01
# timesteps per frame
frame_interval = 10

TestMassColour = 'cyan'

def ivs(r): return r/(np.sum(r**2)**1.5)

class TestMassArray:
    def __init__(self, r, v):
        self.r = np.array(r, dtype=float)
        self.v = np.array(v, dtype=float)
        self.rnew = np.zeros(self.r.shape, dtype=float)
        self.vnew = np.zeros(self.r.shape, dtype=float)
    def dxv(self, r, v, pr, M2):
        rx = np.copy(v)
        vx = np.zeros(v.shape, dtype=float)
        for i in range(v.shape[0]):
            vx[i] = (-G)*(M1*ivs(r[i]) + M2*ivs(pr) - M2*ivs(pr-r[i]))
        return [rx,vx]
    def timestep(self, dt, p):
        rk1, vk1 = self.dxv(self.r, self.v, p.r, p.m)
        rk2, vk2 = self.dxv(self.r + rk1*dt/2.0, self.v + vk1*dt/2.0, p.r + p.rk1*dt/2.0, p.m)
        rk3, vk3 = self.dxv(self.r + rk2*dt/2.0, self.v + vk2*dt/2.0, p.r + p.rk2*dt/2.0, p.m)
        rk4, vk4 = self.dxv(self.r + rk3*dt,     self.v + vk3*dt,     p.r + p.rk3*dt,     p.m)
        self.rnew = dt*(rk1 + rk2*2.0 + rk3*2.0 + rk4)/6.0
        self.vnew = dt*(vk1 + vk2*2.0 + vk3*2.0 + vk4)/6.0
    def pivot(self):
        self.r += self.rnew
        self.v += self.vnew

class Perturber:
    def __init__(self, r, v, M):
        self.m = M
        self.r = np.array(r, dtype=float)
        self.v = np.array(v, dtype=float)
        self.rnew = np.zeros(3, dtype=float)
        self.vnew = np.zeros(3, dtype=float)
        self.rk1,self.vk1 = [np.zeros(3, dtype=float),np.zeros(3, dtype=float)]
        self.rk2,self.vk2 = [np.zeros(3, dtype=float),np.zeros(3, dtype=float)]
        self.rk3,self.vk3 = [np.zeros(3, dtype=float),np.zeros(3, dtype=float)]
        self.rk4,self.vk4 = [np.zeros(3, dtype=float),np.zeros(3, dtype=float)]
    def dxv(self,r,v):
        rx = v
        vx = (-G*(M1+self.m)*ivs(r))
        return [rx,vx]
    def timestep(self, dt):
        self.rk1, self.vk1 = self.dxv(self.r, self.v)
        self.rk2, self.vk2 = self.dxv(self.r + self.rk1*dt/2.0, self.v + self.vk1*dt/2.0)
        self.rk3, self.vk3 = self.dxv(self.r + self.rk2*dt/2.0, self.v + self.vk2*dt/2.0)
        self.rk4, self.vk4 = self.dxv(self.r + self.rk3*dt,     self.v + self.vk3*dt)
        self.rnew = dt*(self.rk1 + self.rk2*2.0 + self.rk3*2.0 + self.rk4)/6.0
        self.vnew = dt*(self.vk1 + self.vk2*2.0 + self.vk3*2.0 + self.vk4)/6.0
    def pivot(self):
        self.r += self.rnew
        self.v += self.vnew

class DoubleSystem:
    def __init__(self, perturber, particles):
        self.perturber = perturber
        self.particles = particles
        self.t = 0.
        self.p0 = ax.scatter(0, 0, 0, color='yellow', s=9.0, edgecolors=None, zorder=10)
        self.p1 = ax.scatter(0, 0, 0, color='red', s=9.0, edgecolors=None, zorder=10)
        self.pt = ax.scatter(0, 0, 0, color='cyan', s=2.5, edgecolors=None, zorder=10)
    def timestep(self, dt):
        self.perturber.timestep(dt)
        self.particles.timestep(dt, self.perturber)
    def pivot(self):
        self.perturber.pivot()
        self.particles.pivot()
    def evolve(self, dt):
        self.timestep(dt)
        self.pivot()
    def plot(self):
        plots=[]
        self.p0.set_offsets([0,0])
        self.p0.set_3d_properties(0, zdir='z')
        plots.append(self.p0)
        self.p1.set_offsets(self.perturber.r[:2])
        self.p1.set_3d_properties(self.perturber.r[2], zdir='z')
        plots.append(self.p1)
        self.pt.set_offsets(self.particles.r[:,:2])
        self.pt.set_3d_properties(self.particles.r[:,2], zdir='z')
        plots.append(self.pt)
        return plots

def circularOrbit(phase, dist):
    cspd = math.sqrt(G*M1/dist)
    r = [math.cos(phase*2.0*math.pi)*dist, math.sin(phase*2.0*math.pi)*dist,0.]
    v = [math.sin(phase*2.0*math.pi)*cspd,-math.cos(phase*2.0*math.pi)*cspd,0.]
    return [r,v]

# Initialise galaxy
#r1 =  2, r2 =  3, r3 =  4, r4 =  5, r5 =  6
#n1 = 12, n2 = 18, n3 = 24, n4 = 30, n5 = 36
#set angular momentum unit vector as [0,0,1]
pp = []
for pt in range(12):
    pp.append(circularOrbit(pt/12.0, 2.0))
for pt in range(18):
    pp.append(circularOrbit(pt/18.0, 3.0))
for pt in range(24):
    pp.append(circularOrbit(pt/24.0, 4.0))
for pt in range(30):
    pp.append(circularOrbit(pt/30.0, 5.0))
for pt in range(36):
    pp.append(circularOrbit(pt/36.0, 6.0))
pn = np.array(pp)

testParticles = TestMassArray(pn[:,0], pn[:,1])
perturberGalaxy = Perturber([-20,20,20],[5,-5.,-5],10)

plt.ioff()
plt.style.use('dark_background')
fig = plt.figure(figsize=[6, 6])
ax = fig.add_subplot(111,projection='3d')
ax.set_xlim3d([-18,18])
ax.set_ylim3d([-18,18])
ax.set_zlim3d([-18,18])

ds = DoubleSystem(perturberGalaxy, testParticles)

def animate(i):
    for n in range(frame_interval):
        ds.evolve(dt)
    return ds.plot()
ani = animation.FuncAnimation(fig, animate, repeat=False, frames=sim_duration, blit=False, interval=30,)
plt.show()
