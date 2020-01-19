#!/usr/bin/python3
# -*- coding: utf-8 -*-

from numpy import *
from matplotlib.pyplot import *
import cmath
import scipy.interpolate
from multiprocessing import Pool, cpu_count, Manager 


# ============================================
# ==============   Define Parameters  =================
# ============================================


# layer parameters
ER = [2.0]	   					# relative permittivity of each layer
MR = [1.0]						# relative permeability of each layer
L = [500e-9/(4*sqrt(2))] 	# thickness of each layer in m

# parameters of incidence medium
er_1 = 1
mr_1 = 1

# parameters of transmission medium
er_2 = 4
mr_2 = 1

# incidence angles in ° (theta can be a vector)
theta = 45
phi   = 0   

# free-space wavelength in m (can be a vector)
lam0  = linspace(50,2000,1000)*1e-9 

# Data file provides the permittivity/permeability (NK = 0) or the refractive index n/k (NK = 1)
NK = 1








# ============================================
# ===============   Computation   ===================
# ============================================

Ncpu = cpu_count()
manager = Manager()
thetaSplit = manager.list()
R_TE = []
T_TE = []
R_TM = []
T_TM = []
I = array([[1, 0],[0, 1]])

# --------------------------------------------------------------------------------------------------------------------

def Inv2(A):
	
	return(1/(A[0,0]*A[1,1] - A[0,1]*A[1,0])*array([ [A[1,1], -A[0,1] ], [-A[1,0], A[0,0] ] ]))


# --------------------------------------------------------------------------------------------------------------------

def redheffer(SA11,SA12,SA21,SA22,SB11,SB12,SB21,SB22):

	D = SA12.dot(Inv2(I - SB11.dot(SA22)))
	F = SB21.dot(Inv2(I - SA22.dot(SB11)))

	SAB11 = SA11 + D.dot(SB11.dot(SA21))
	SAB12 = D.dot(SB12)
	SAB21 = F.dot(SA21)
	SAB22 = SB22 + F.dot(SA22.dot(SB12))

	return(SAB11,SAB12,SAB21,SAB22)

# --------------------------------------------------------------------------------------------------------------------

def Interp1(d, lam):
	
	if type(d) == str:
	
		with open(d, 'r') as f:
			data 	= loadtxt(d, comments='#')
			
			lamd 	= data[:,0]
			d1		= data[:,1]
			d2		= data[:,2]

			if NK:		# data given in terms of n and k
				rel	= (d1**2 - d2**2) + 1j*(2*d1*d2)
			else:		# data given in terms of er_r and er_i
				rel	= d1 + 1j*d2
			
		# interpolate data
		interp = scipy.interpolate.interp1d(lamd, rel)
		
		return interp(lam)
			
	else:
		return  d
	
	
# --------------------------------------------------------------------------------------------------------------------


def tmm1d(ER, MR, L, er1, mr1, er2, mr2, theta, phi, pte, ptm, lam):
	
	theta = theta*pi/180
	phi = phi*pi/180

	k0 = 2*pi/lam
	ninc = cmath.sqrt(er1*mr1)

	# normalized kinc vector (k0 is missing)
	kinc  = ninc*array([sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)])

	kx = kinc[0]
	ky = kinc[1]

	# ------------------------------------------
	# Calculate Gap Medium Parameters

	Qg = array([[kx*ky, 1+ky**2], [-(1+kx**2),  -kx*ky]])

	Vg = -1j*Qg

	# ------------------------------------------
	# Calculate source

	# vector normal to medium
	nvec = array([0, 0, -1])

	# definition of polarization vector (see Lecture 3 CEM)
	if theta == 0:
		ate = array([0, 1, 0])
	else:
		ate = cross(kinc,nvec)/linalg.norm(cross(kinc,nvec))

	atm = cross(ate,kinc)/linalg.norm(cross(ate,kinc)) 

	P = pte*ate + ptm*atm

	Ei = array([[P[0]], [P[1]]])

	# ------------------------------------------
	# Initialize Global Scattering Matrix

	S11G = array([[0, 0], [0, 0]])
	S22G = array([[0, 0], [0, 0]])
	S21G = I
	S12G = I

	# ------------------------------------------
	# Loop through layers

	for i in range(0, len(ER)):

		er = Interp1(ER[i],lam)
		mr = Interp1(MR[i],lam)

		# -------------------------------------------------------
		# Parameters for Layer i

		kzi = cmath.sqrt(mr*er - kx**2 - ky**2)
		
		if abs(kzi) == 0:
			kzi = 1e-50

		Qi = 1/mr*array([[kx*ky, mr*er - kx**2], [ky**2-mr*er, -kx*ky]])

		Omegai = 1j*kzi*I

		Vi = Qi.dot(Inv2(Omegai))

		# -------------------------------------------------------
		# Scattering Matrix for Layer i

		Ai   = I + Inv2(Vi).dot(Vg)
		Bi   = I - Inv2(Vi).dot(Vg)
		Xi   = array([[exp(Omegai[0,0]*k0*L[i]), 0], [0, exp(Omegai[1,1]*k0*L[i])]])

		D    = Ai - Xi.dot(Bi.dot(Inv2(Ai).dot(Xi.dot(Bi))))

		S11i = Inv2(D).dot(Xi.dot(Bi.dot(Inv2(Ai).dot(Xi.dot(Ai)))) - Bi)
		S12i = Inv2(D).dot(Xi.dot(Ai - Bi.dot(Inv2(Ai).dot(Bi))))
		S22i = S11i
		S21i = S12i

		S11G, S12G, S21G, S22G = redheffer(S11G,S12G,S21G,S22G,S11i,S12i,S21i,S22i)



	# ------------------------------------------
	# Connect to reflection region

	kzr = cmath.sqrt(er1*mr1 - kx**2 - ky**2)
	
	if abs(kzr) == 0:
		kzr = 1e-50

	Qr = 1/mr1*array([[kx*ky, mr1*er1 - kx**2], [ky**2-mr1*er1, -kx*ky]])

	Omegar = 1j*kzr*I

	Vr = Qr.dot(Inv2(Omegar))

	A   = I + Inv2(Vg).dot(Vr)
	B   = I - Inv2(Vg).dot(Vr)
	S11r = -Inv2(A).dot(B)
	S12r = 2*Inv2(A)
	S21r = 1/2*(A - B.dot(Inv2(A).dot(B)))
	S22r = B.dot(Inv2(A))

	S11G, S12G, S21G, S22G = redheffer(S11r,S12r,S21r,S22r,S11G,S12G,S21G,S22G)

	# ------------------------------------------
	# Connect to transmission region
	
	kzt = cmath.sqrt(er2*mr2 - kx**2 - ky**2)
		
	if abs(kzt) == 0:
		kzt = 1e-50

	Qt = 1/mr2*array([[kx*ky, mr2*er2 - kx**2], [ky**2-mr2*er2, -kx*ky]])
	Omegat = 1j*kzt*I;
	Vt = Qt.dot(Inv2(Omegat))

	A   = I + Inv2(Vg).dot(Vt)
	B   = I - Inv2(Vg).dot(Vt)
	S11t = B.dot(Inv2(A))
	S12t = 1/2*(A - B.dot(Inv2(A).dot(B)))
	S21t = 2*Inv2(A)
	S22t = -Inv2(A).dot(B)

	S11G, S12G, S21G, S22G = redheffer(S11G,S12G,S21G,S22G,S11t,S12t,S21t,S22t)


	# ------------------------------------------

	# Calculate transmitted and reflected fields
	Er = S11G.dot(Ei)
	Et = S21G.dot(Ei)

	# Calculate longitudinal components
	Ezi = -(kx*Ei[0] + ky*Ei[1])/kinc[2]
	Ezr = -(kx*Er[0] + ky*Er[1])/kzr
	Ezt = -(kx*Et[0] + ky*Et[1])/kzt

	# Calculate reflectance and transmittance
	R = (abs(Er[0])**2 + abs(Er[1])**2 + abs(Ezr)**2)/(abs(Ei[0])**2 + abs(Ei[1])**2 + abs(Ezi)**2);
	T = (abs(Et[0])**2 + abs(Et[1])**2 + abs(Ezt)**2)/(abs(Ei[0])**2 + abs(Ei[1])**2 + abs(Ezi)**2)*real(kzt/mr2/kinc[2]*mr1);

	return(R,T)



# --------------------------------------------------------------------------------------------------------------------


def InterpData(d):
		
	if type(d) == str:
	
		with open(d, 'r') as f:
			data 	= loadtxt(d, comments='#')
			
			lamd 	= data[:,0]
			d1		= data[:,1]
			d2		= data[:,2]

			if NK:		# data given in terms of n and k
				rel	= (d1**2 - d2**2) + 1j*(2*d1*d2)
			else:		# data given in terms of er_r and er_i
				rel	= d1 + 1j*d2
			
		# interpolate data
		interp = scipy.interpolate.interp1d(lamd, rel)
		
		if size(lam0) == 1:
			return [interp(lam0)]
		else:
			return interp(lam0)
			
	else:
		return  full(size(lam0), d)



# --------------------------------------------------------------------------------------------------------------------





def tmm1dLam():
	
	global R_TE, R_TM, T_TM, T_TE, er1, er2, mr1, mr2
	
	for i in range(len(lam0)):

		r, t = tmm1d(ER, MR, L, er1[i], mr1[i], er2[i], mr2[i], theta, phi, 1, 0, lam0[i])

		R_TE.append(r)
		T_TE.append(t)

		r, t = tmm1d(ER, MR, L, er1[i], mr1[i], er2[i], mr2[i], theta, phi, 0, 1, lam0[i])

		R_TM.append(r)
		T_TM.append(t)
	
	
	fig, ax = subplots()
	
	plot(lam0,T_TM,'r',label="Tp")    
	plot(lam0,R_TM,'--r',label="Rp")
	plot(lam0,T_TE,'b',label="Ts")
	plot(lam0,R_TE,'--b',label="Rs")
	legend()
	grid()
	xlabel("Wavelength (m)")
	ylabel("Reflectance & Transmittance")
	gca().set_xlim(lam0[0], lam0[-1])
	gca().set_ylim(0, 1.01)
	
	ax.tick_params(axis='both', which='major')
	ticklabel_format(style='sci', axis='x', scilimits=(0,3))
	
	show()
	
	
# --------------------------------------------------------------------------------------------------------------------
	
	
def tmm1dAng():
	
	global R_TE, R_TM, T_TM, T_TE, er1, er2, mr1, mr2

	for tet in theta:

		r, t = tmm1d(ER, MR, L, er1[0], mr1[0], er2[0], mr2[0], tet, phi, 1, 0, lam0)

		R_TE.append(r)
		T_TE.append(t)

		r, t = tmm1d(ER, MR, L, er1[0], mr1[0], er2[0], mr2[0], tet, phi, 0, 1, lam0)

		R_TM.append(r)
		T_TM.append(t)
	
	plot(theta,T_TM,'r',label="Tp")    
	plot(theta,R_TM,'--r',label="Rp")
	plot(theta,T_TE,'b',label="Ts")
	plot(theta,R_TE,'--b',label="Rs")
	legend()
	grid()
	xlabel("Incidence angle (°)")
	ylabel("Reflectance & Transmittance")
	gca().set_xlim(theta[0], theta[-1])
	gca().set_ylim(0, 1.01)
	
	show()


# --------------------------------------------------------------------------------------------------------------------


def Para(p):
	
	global er1, er2, mr1, mr2
	RT = []

	for tet in thetaSplit[p]:
		for i in range(len(lam0)):
		
			r_te, t_te = tmm1d(ER, MR, L, er1[i], mr1[i], er2[i], mr2[i], tet, phi, 1, 0, lam0[i])
			r_tm, t_tm = tmm1d(ER, MR, L, er1[i], mr1[i], er2[i], mr2[i], tet, phi, 0, 1, lam0[i])
			
			RT.append( [r_te[0], t_te[0], r_tm[0], t_tm[0]] )

	return RT


# --------------------------------------------------------------------------------------------------------------------


def tmm2d():

	global R_TE, R_TM, T_TM, T_TE, er1, er2, mr1, mr2
	
	# Create pool for the number of cpu
	pool = Pool(Ncpu)

	# split wavelength for each cpu
	N = int(round(len(theta)/Ncpu))
	Nn = len(theta) % Ncpu
	
	n = 0	
	while n < Ncpu:
		if Nn and n == Ncpu - 1 and theta[n*N : N*(n+1) + Nn].any():
			thetaSplit.append( theta[n*N: N*(n+1) + Nn] )
		elif theta[n*N : N*(n+1)].any():
			thetaSplit.append( theta[n*N: N*(n+1)  ] )
		n += 1

	# Compute R and T in parallel
	RT = pool.map(Para, range(len(thetaSplit)))
		
	pool.close()
	pool.join()
	

	for n in RT:
		for m in n:
			R_TE.append(m[0])
			T_TE.append(m[1])
			R_TM.append(m[2])
			T_TM.append(m[3])
	
	
	R_TE = transpose(reshape(R_TE, (len(theta),len(lam0))))
	T_TE = transpose(reshape(T_TE, (len(theta),len(lam0))))
	R_TM = transpose(reshape(R_TM, (len(theta),len(lam0))))
	T_TM = transpose(reshape(T_TM, (len(theta),len(lam0))))
	
	
	fig, ax = subplots(2,2,figsize=(10,8))
	
	subplot(2,2,1)
	imshow(R_TM, extent = [theta[0], theta[-1], lam0[-1],  lam0[0]], vmin = 0, vmax = 1, cmap=cm.jet, aspect='auto')
	colorbar()
	xlabel("Incidence angle [°]")
	ylabel("Wavelength [m]")
	title("Rp")
	ticklabel_format(style='sci', axis='y', scilimits=(0,3))
	
	subplot(2,2,2)
	imshow(T_TM, extent = [theta[0], theta[-1], lam0[-1],  lam0[0]], vmin = 0, vmax = 1, cmap=cm.jet, aspect='auto')
	colorbar()
	xlabel("Incidence angle [°]")
	ylabel("Wavelength [m]")
	title("Tp")
	ticklabel_format(style='sci', axis='y', scilimits=(0,3))
	
	subplot(2,2,3)
	imshow(R_TE, extent = [theta[0], theta[-1], lam0[-1],  lam0[0]], vmin = 0, vmax = 1, cmap=cm.jet, aspect='auto')
	colorbar()
	xlabel("Incidence angle [°]")
	ylabel("Wavelength [m]")
	title("Rs")
	ticklabel_format(style='sci', axis='y', scilimits=(0,3))
	
	subplot(2,2,4)
	imshow(T_TE, extent = [theta[0], theta[-1], lam0[-1],  lam0[0]], vmin = 0, vmax = 1, cmap=cm.jet, aspect='auto')
	colorbar()
	xlabel("Incidence angle [°]")
	ylabel("Wavelength [m]")
	title("Ts")
	ticklabel_format(style='sci', axis='y', scilimits=(0,3))

	
	fig.subplots_adjust(hspace=.3)
	

	show()


# --------------------------------------------------------------------------------------------------------------------


def main():
	
	global er1, er2, mr1, mr2
	
	er1 	= InterpData(er_1)
	mr1 	= InterpData(mr_1)
	er2 	= InterpData(er_2)
	mr2 	= InterpData(mr_2)
		
	if size(lam0) > 1 and size(theta) == 1:
		tmm1dLam()
	elif size(lam0) == 1 and size(theta) > 1:
		tmm1dAng()
	elif size(lam0) > 1 and size(theta) > 1:
		tmm2d()
		



main()
