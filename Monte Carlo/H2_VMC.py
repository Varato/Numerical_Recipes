import numpy as np
from numpy import linalg as la
import sys
import random
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D


global n_walkers  # Number of walker
global step_size  # The max step in a walk
global r10        # old coordinates of 2 electrons
global r20
global E_sum, Esqrd_sum, nAccept
global E_ave, E_var
global MCSteps
global sample

R_norm = 1.398 # in unit of Bohr radias
R = np.array([0,0,R_norm],dtype=float)

E_sum = Esqrd_sum = nAccept = 0
MCSteps = 5000
n_walkers = 150
step_size = 1.0
xMin = -10.0
xMax = 10.0

def flushPrint(string,num):
	sys.stdout.write('\r')
	sys.stdout.write(string+'%s' % num)
	sys.stdout.flush()

def initialize():
	global n_walkers
	global r10, r20
	r10 = np.zeros([n_walkers,3])
	r20 = np.zeros([n_walkers,3])
	for walker in range(n_walkers):
		r10[walker]=np.random.uniform(-0.5,0.5,3)
		r20[walker]=np.random.uniform(-0.5,0.5,3)

def zeroAccumulate():
	global E_sum, Esqrd_sum
	E_sum=0
	Esqrd_sum=0

def trial_wave_func(r1, r2):
	rA1 = la.norm(r1-R)
	rA2 = la.norm(r2-R)
	rB1 = la.norm(r1)
	rB2 = la.norm(r2)
	psi = np.exp(-rA1-rB2)+np.exp(-rA2-rB1)
	return psi

def local_energy(r1, r2):
	rA1 = la.norm(r1-R)
	rA2 = la.norm(r2-R)
	rB1 = la.norm(r1)
	rB2 = la.norm(r2)
	rAB = R_norm
	r12 = la.norm(r1-r2)
	T_up = np.exp(-rB2-rA1)*(-1.+1./rA1+1./rB2) + np.exp(-rA2-rB1)*(-1.+1./rB1+1./rA2)
	T = T_up/trial_wave_func(r1, r2)
	V = -1./rA1-1./rB1-1./rA2-1./rB2+1./rAB+1./r12
	E_loc = T+V
	return E_loc

def Metropolis_step(walker):
	global r10, r20, nAccept, E_sum, Esqrd_sum, alpha, sample
	trial_r1 = r10[walker]+step_size*np.random.uniform(-1, 1, 3)
	trial_r2 = r20[walker]+step_size*np.random.uniform(-1, 1, 3)
	# if trial_x < xMin: trial_x = xMin
	# if trial_x > xMax: trial_x = xMax
	w0 = trial_wave_func(r10[walker], r20[walker])
	w1 = trial_wave_func(trial_r1, trial_r2)
	p = (w1/w0)**2
	# Is this trial acceptable?
	if random.random() < p: 
		# update
		r10[walker] = trial_r1
		r20[walker] = trial_r2
		nAccept += 1
	E_local = local_energy(r10[walker],r20[walker])
	E_sum = E_sum+E_local
	Esqrd_sum = Esqrd_sum+E_local**2

def oneMento_Carlo_step():
	global n_walkers
	for walker in range(n_walkers):
		Metropolis_step(walker)

def runMonte_Carlo():
	global MCSteps, step_size, n_walkers
	global E_ave, E_var, nAccept
	# Firstly use 20% MCSteps to adjust step_size so that 
	# the acceptance ratio is about 50%
	thermSteps = int(0.2*MCSteps)
	adjust_interval = int(0.1*thermSteps)+1
	nAccept = 0

	print "-----------------------------"
	print "Thermalizing...."
	for i in range(thermSteps):
		flushPrint("current step: ", i+1)
		oneMento_Carlo_step()
		if (i+1) % adjust_interval == 0:
			step_size *= nAccept/(0.5*n_walkers*adjust_interval)
			nAccept = 0
	print
	print "Adjusted step_size = %f" % step_size
	print "-----------------------------"

	#After adjust step_size, let's get E_ave and E_var
	zeroAccumulate()
	nAccept=0
	print "Start real computation:"
	for i in range(MCSteps):
		flushPrint("current step: ", i+1)
		oneMento_Carlo_step()
	print

	E_ave = 1.0*E_sum/n_walkers/MCSteps
	E_var = 1.0*Esqrd_sum/n_walkers/MCSteps - E_ave**2

if __name__=="__main__":
	initialize()
	
	runMonte_Carlo()
	print "Acceptance ratio: ", 1.0*nAccept/n_walkers/MCSteps
	print "E_ave = %f, E_var = %f" % (E_ave, E_var)

	# plt.subplot(121)
	# plt.plot(a, Ev, 'r*-')
	# plt.xlabel("Parameter")
	# plt.ylabel("Energy")
	# plt.yticks([0.,.2,.4,.5,.6,.8,1.,1.2,1.4])
	# plt.plot([0,1.4],[0.5,0.5],'k-.')
	# # plt.axis([0, 1.4, 0, 1.4])
	# plt.subplot(122)
	# plt.plot(a, var, 'b*-')
	# plt.xlabel("Parameter")
	# plt.ylabel("Varience")
	# plt.show()




