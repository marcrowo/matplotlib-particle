import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

class Particle:
	def __init__(self, pos, vel, mass, charge):
		self.pos = pos
		self.vel = vel
		self.mass = mass
		self.charge = charge
class Dipole:
	def __init__(self, mu, source, moment):
		self.mu = mu
		self.source = source
		self.moment = moment
	def flux_density(self, x):
	#http://en.wikipedia.org/wiki/Magnetic_moment
		r = np.subtract (x, self.source)
		d = np.linalg.norm(r)
		d3 = d*d*d
		d5 = d3*d*d
		p = 3.0 * np.dot(self.moment, r) / d5
		b1 = r * p
		b2 = self.moment / d3
		b1 = np.subtract(b1, b2)
		b1 *= self.mu
		return b1

part1 = Particle([0.0, 0.0, 0.0], [0.5, 0.0, 0.0], 1.0, 1.0)
#part2 = Particle([0.0, 10.0, 0.0], [1.0, 0.0, 0.0], 1.0, 1.0)
#part3 = Particle([0.0, 0.0, 2.0], [1.0, 0.0, -0.02], 1.0, 1.0)

particle_group = [part1] #, part2, part3]

dipole = Dipole(1.0, [100.0, 0.0, 0.0], [200.0, 0.0, 0.0])

#ELECTRIC FIELD
e_field = [0.0, 0.0, 0.0]
magn_uniform_field = [0.0, 0.0, 0.0]

#MATPLOTLIB PARAMS
fig = plt.figure()
#ax = fig.gca(projection='3d') // INCOMPATIBILITY ISSUE??
ax = Axes3D(fig)

#Euler method
def euler(particle, max_iter):
	q_over_m = particle.charge / particle.mass
	results = []
	i = 0
	
	for i in range(max_iter):
		i += 1
		p = particle.pos
		v = particle.vel
		b = dipole.flux_density(p) + magn_uniform_field
		a = np.cross(v, b) + e_field
		a = a * q_over_m
		v = np.add(v, a)
		p = np.add(p, v)
		particle.pos = p
		particle.vel = v
		#print p
		results.append(p)
	return results
	
#Runge-Kutta method
def rk4(particle, max_iter):
	q_over_m = particle.charge / particle.mass
	results = []
	i = 0
	
	for i in range(max_iter):
		i += 1
		p1 = np.array(particle.pos)
		v1 = np.array(particle.vel)
		b = dipole.flux_density(p1) + magn_uniform_field
		a1 = np.cross(v1, b) + e_field
		a1 = a1 * q_over_m
		
		p2 = p1 + (v1 * 0.5)
		v2 = v1 + (a1 * 0.5)
		b = dipole.flux_density(p2) + magn_uniform_field
		a2 = np.cross(v2, b) + e_field
		a2 = a2 * q_over_m
		
		p3 = p1 + (v2 * 0.5)
		v3 = v1 + (a2 * 0.5)
		b = dipole.flux_density(p3) + magn_uniform_field
		a3 = np.cross(v3, b) + e_field
		a3 = a3 * q_over_m
		
		p4 = p1 + v3
		v4 = v1 + a3
		b = dipole.flux_density(p4) + magn_uniform_field
		a4 = np.cross(v4, b) + e_field
		a4 = a4 * q_over_m
		
		dv = a1 + 2.0 * a2 + 2.0 * a3 + a4 
		v = v1 + dv / 6.0
		
		dp = v1 + 2.0 * v2 + 2.0 * v3 + v4
		p = p1 + dp / 6.0
		
		particle.pos = p
		particle.vel = v
		#print p
		results.append(p)
	return results
	

#Plotting
def part_plot(particle, max_iter, method, color):
	
	ax.scatter(dipole.source[0], dipole.source[1], dipole.source[2], color = 'red')
	
	x = []
	y = []
	z = []
	
	if method == 'euler':
		results = euler(particle, max_iter)
	elif method == 'rk4':
		results = rk4(particle, max_iter)
	
	for p in results:
		x.append(p[0])
		y.append(p[1])
		z.append(p[2])	
	ax.scatter(x, y, z, color = color)

for part in particle_group:
    colors = ['green', 'purple', 'blue', 'yellow']
    part_plot(part, 500, 'rk4', colors[particle_group.index(part)])

#part_plot(part1, 160, 'rk4', 'green')
#part_plot(part2, 500, 'euler', 'purple')
#part_plot(part2, 160, 'blue')
#part_plot(part3, 800, 'rk4', 'purple')

plt.show()