import matplotlib.pyplot as plt
import numpy as np

### initial conditions ###
p0 = 10
q0 = 0
m = 1
deltat = 0.5
T = 10
n = int(T/deltat)

def energy(q,p):
  return 0.5* p**2 / m + 0.5* q**2 *m

### analytic solutions ###
def q_analytic(t):
  return p0/m * np.sin(t) + q0 * np.cos(t)
def p_analytic(t):
  return p0*np.cos(t) - m*q0*np.sin(t)
def dqdt_analytic(t):
  return p0* np.cos(t) - m*q0* np.sin(t)
times = np.arange(0,T,deltat/100)
energy_analytic = energy(q_analytic(times),dqdt_analytic(times))

### Euler ###
def euler(q0,p0,deltat,T):
  q = np.zeros([n])
  p = np.zeros([n])
  p[0] = p0
  q[0] = q0
  for i in range(0,n-1):
    p[i+1] = p[i] - q[i]*m*deltat
    q[i+1] = q[i] + deltat * p[i] /m # (more stable if one takes p[i+1] (=simplectic Euler))
  return q, p
q_euler, dqdt_euler = euler(q0,p0,deltat,T)
energy_euler = energy(q_euler,dqdt_euler)

### Runge Kutta 4th order ###
# u = [q, p] = [q dqdt] (ungefähr)
u0 = [p0,q0]

def F(u,t): # this actually does not depend on t for our case
  dqdt = u[1]/m
  dpdt = -u[0]*m
  return np.array([dqdt,dpdt])

def rk4(u0,deltat,T, F):
  u = np.zeros([n,2])
  u[0,0] = q0
  u[0,1] = p0
  t=1
  for i in range(0,n-1):
    k1 = F(u[i,:], t)
    k2 = F(u[i,:] + 0.5*deltat* k1,t + 0.5*deltat)
    k3 = F(u[i,:] + 0.5*deltat* k2,t + 0.5*deltat)
    k4 = F(u[i,:] + deltat* k3,t + deltat)

    u[i+1,:] = u[i,:] + deltat*(1/6*k1 + 1/3*k2 +1/3*k3 + 1/6*k4)
  return u

u_rk4 = rk4(u0,deltat,T,F)
q_rk4 = u_rk4[:,0]
dqdt_rk4 = u_rk4[:,1]
energy_rk4 = energy(q_rk4,dqdt_rk4)

### Plotting the solutions in a phase path diagram ###
plt.plot(q_analytic(times), p_analytic(times), label ='analytic')
plt.plot(q_rk4, dqdt_rk4, label ='RK4')
plt.plot(q_euler, dqdt_euler, label ='Euler')
plt.xlabel('q(t)')
plt.ylabel('p(t)')
plt.axis('square')
plt.xlim(-12,12)
plt.ylim(-12,12)
plt.legend()
plt.grid(color = 'gainsboro')
plt.savefig("plot-phasepath.pdf")

### Plotting the solutions in a energy evolution diagram ###
# plt.plot(times, energy_analytic, label ='analytic')
# plt.plot(times, energy_rk4, label ='RK4')
# plt.plot(times, energy_euler, label ='Euler')
# plt.xlabel('t')
# plt.ylabel('energy(t)')
# # plt.axis('square')
# # plt.xlim(-12,12)
# plt.ylim(49.94,50.008)
# plt.legend()
# plt.grid(color = 'gainsboro')
# plt.savefig("plot-energy.pdf")

### Plotting the solutions in a position-time diagram, or possibly the error ###
# plt.plot(times, q_analytic(times), label ='analytic')
# plt.plot(times, dqdt_analytic(times), label ='analytic')
# plt.plot(times, q_rk4, label ='RK4')
# plt.plot(times, q_analytic(times)-q_rk4) #, label = "error = q_analytic - q_rk4")
# # plt.plot(times, dqdt_rk4, label ='RK4 velocity')
# # plt.plot(times, q_euler, label ='Euler')
# plt.xlabel('t')
# plt.ylabel('error = q_analytic - q_rk4')
# # plt.axis('square')
# # plt.xlim(-12,12)
# # plt.ylim(49.94,50.008)
# plt.legend()
# plt.grid(color = 'gainsboro')
# plt.savefig("plot-error.pdf")