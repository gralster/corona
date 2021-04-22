# a program in 1 dimension to simulate the sun's corona
# using numerical equations of hydrodynamics
# Grace Alster

import numpy as np
import math
import matplotlib.pyplot as plt
import time
from options import *

#arraysetup ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃
#Createinitialarraystobefilledthroughouttherun
R=np.zeros((nlevels,maxt))	#position
U=np.zeros((nlevels,maxt))	#velocity
p=np.zeros((nlevels+1,maxt))#pressure
q=np.zeros((nlevels+1,maxt))#pseudoviscouspressure
V=np.zeros((nlevels+1,maxt))#specif icvolume
d=np.zeros((nlevels+1,maxt))#density
E=np.zeros((nlevels+1,maxt))#energy
T=np.zeros((nlevels+1,maxt))#temperature
timesteplog=np.zeros(maxt-1)
innermass=np.zeros(maxt)

# ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃
#Difference equations and  Stepping helper functions
# ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃

#return sdU/dt-optiontoaddaccelerationduetogravityattop
def  stepVelocity(x,t):
	delq=q[x+1][t-1]-q[x][t-1]
	delp=p[x+1][t-1]-p[x][t-1]
	if gravity==True:
		print((G*sm)/(R[x][t-1]**2))
		return (-1*(1/d[x][0])*((delp+delq)/(R[x][0]-R[x-1][0]))*(R[x][t-1]/R[x][0])**2)-(G*sm)/(R[x][t-1]**2)
	else:
		return (-1*(1/d[x][0])*((delp+delq)/(R[x][0]-R[x-1][0]))*(R[x][t-1]/R[x][0])**2)

#SPECif ICVOLUMEdif ferenceequation
def  getSpecVol(x,t):
	return ((1/d[x][0])*((R[x][t]**3-R[x-1][t]**3)/(R[x][0]**3-R[x-1][0]**3)))
#redirectionfunctiontodealwithpressureboundarycondition
def  stepPressure(x,t):
	if  zeroOuterPGrad==True  and  x==nlevels-1:
		return p[nlevels][t]
	else:
		return setPressure(x,t)
#setspressureacc or dingtotheadiabaticequationofstate
#Pressure=const*densityˆ(gamma)
#wheregamma=theadiabaticindex
def  setPressure(x,t):
	#print(x)
	#print(p[x][t-1]) 
	#print(d[x][t])
	#if p[x][t-1] >= maxP:
	#	return maxP
	#else:
	return  k[x]*(d[x][t]**(adiabaticindex))
#Energyequationdif ferenceequation
def stepEnergy(x,t):
	return (E[x][t-1]-(0.5*p[x][t-1]+q[x][t])*(V[x][t]-V[x][t-1]))/(1+(1/3)*(V[x][t]-V[x][t-1]))
#return sE=3PV/2f or settinginitialconditions
def setinitEnergy(x,t):
	return (3/2)*p[x][t]*V[x][t]
#calculatespseudoviscosityif theaccelerationis-ve
def PseudoVisc(x,t):
	delu=U[x][t]-U[x-1][t]
	q=(((2*a**2)/(V[x][t]+V[x][t-1]))*delu**2)
	if delu<0:
		#print("need visc")
		if q<0:
			return q*-1
		else:
			return q
	else:
		return 0
# ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃
# Mainengine - control of timestep and calls difference equations
# ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃

def newTimeStep(t):
	# only occurs if constantTimeStep = False
	#set up quantity gradient arrays 
	deld=[0]*(nlevels-2)
	v=[0]*(nlevels-2)
	delu=[0]*(nlevels-1)
	delp=[0]*(nlevels-1)
	delv=[0]*(nlevels-1)
	#Conform timestep to fastest moving zone
	for i in range(nlevels-2):
		#find speed of each zone
		if adiabaticindex*p[i][t]/d[i][t]>0:
			v[i]=0.5*(U[i][t]+U[i+1][t])+math.sqrt(adiabaticindex*p[i][t]/d[i][t])
			if v[i]<0:
				v[i]=v[i]*-1
		else:
			pass
		deld[i]=(R[i+1][t]-R[i][t]) #find distance between zones
	optionst=np.array(deld)/v		#time = distance/speed
	#print(optionst)
	mins=[min(optionst)]			#choose shortest time
	#Conform timestep to fastest moving rate of quantity input
	if pressurerate>0: # if pressure is being added, add the gradient of pressure change into consideration
		mins.append(p[pressurezone][t]/(pressurerate*p[pressurezone][0]))
	if massrate>0:	# if mass is being added, add the gradient of mass change into consideration
		mins.append(innermass[t]/(massrate*innermass[0]))
	if energyrate>0: #  if energy is being added, add the gradient of energy change into consideration
		mins.append(E[energyzone][t]/(energyrate*E[energyzone][0]))
	minimum=min(mins) #choose smallest time step option
	if minimum<0:#preventscrashbutshouldnotbehappening
		minimum=minimum*-1
		print("MinimumTimestepisnegative!")
	return minimum


def stepForwardCore(x,t,timestep):
	#Updates velocity and position, happens every timestep
	#Core "on" levels
	U[x][t+1]=U[x][t]+stepVelocity(x,t+1)*timestep #velocity
	#findthenextvelocity
	if ShockPiston==True: #Move zone 1 artificially for piston experiment
		U[1][t+1]=U[1][t]
	#wind addition of mass if on
	if massrate>0 and x==masszone:
		delq=q[x+1][t]-q[x][t]
		delp=p[x+1][t]-p[x][t]
		dr=innermass[t+1]
		U[x][t+1]=U[x][t]-((delp+delq)/(dr)*(R[x][t])**2)*timestep

	R[x][t+1]=R[x][t]+U[x][t+1]*timestep #update position
	
def stepForwardPhysics(x,t,timestep):
	#Update specific volume, pseudoviscosity, density, pressure and energy
	#physics "between" levels
	V[x][t+1]=getSpecVol(x,t+1)
	#wind addition of mass if on
	if massrate>0 and (x==masszone):# or x==masszone+1):
		V[x][t+1]=((R[x][t]**3-R[x-1][t]**3)/(innermass[t+1]))
	q[x][t+1]=PseudoVisc(x,t)
	d[x][t+1]=(1/V[x][t+1])
	p[x][t+1]=stepPressure(x,t+1)
	E[x][t+1]=setinitEnergy(x,t+1)
	#wind addition of energy and pressure options
	if useEnergyEq==True and x==energyzone:
		E[x][t+1]=stepEnergy(x,t+1)
		p[x][t+1]=(2/3)*(E[x][t+1]/V[x][t+1])+(energyrate*E[x][0]*timestep)/(R[x][t+1]**3-R[x-1][t+1]**3)
	elif pressurerate>0 and (x==masszone):# or x==masszone+1):
		p[x][t+1]=p[x][t]+p[x][0]*pressurerate*timestep

# ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃
#RUNPREPARATION+CONDITIONSSETUP
# ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃

def fixedBC():
	#Sets boundary conditions from settings in options file
	for t in range(maxt):
		R[nlevels-1][t]=outerposition
		R[0][t]=x0
		d[0][t]=d0
		d[nlevels][t]=outerdensity
		V[0][t]=1/d[0][t]
		V[nlevels][t]=1/d[nlevels][t]
		p[0][t]=pinner
		p[nlevels][t]=setPressure(nlevels,t)
		q[0][t]=0
		T[0][t]=(pinner)/(d0)
		T[nlevels][t]=(outerpressure)/(outerdensity)
	return R,U,V,d,p,q,T

def setUp():
	#sets initial conditions from options file
	#LEVELS
	for x in range(nlevels):
		R[x][0]=x0+x*spacestep#def inefirstposition
		U[x][0]=u0
	#BETWEENLEVELS
	for x in range(nlevels+1):
		d[x][0]=d0
		p[x][0]=setPressure(x,0)
		q[x][0]=0
		V[x][0]=1/d0
		T[x][0]=(p[x][0])/(d[x][0])
		R[nlevels-1][0]=outerposition
		E[energyzone][0]=setinitEnergy(energyzone,0)
		innermass[0]=(4*math.pi*d[masszone][0])*(R[masszone][0]**3-R[
		masszone-1][0]**3)
	return innermass,E,T

# ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃
#Write out all experiment data to file a long with settings used
# ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃
def writeout(outputname,pos,vel,pres,pseu,Vol,den,En,time,ma,Temp):
	f=open(outputname+"settings.txt","w")
	f.write("constantTimestep="+str(constantTimestep)+"\n")
	f.write("constTimestep="+str(constTimestep)+"\n")
	f.write("zeroOuterPGrad="+str(zeroOuterPGrad)+"\n")
	f.write("k="+str(k[0])+"\n")
	f.write("energyrate="+str(energyrate)+"\n")
	f.write("energyzone="+str(energyzone)+"\n")
	f.write("pressurate="+str(pressurerate)+"\n")
	f.write("pressurezone="+str(pressurezone)+"\n")
	f.write("massrate="+str(massrate)+"\n")
	f.write("masszone="+str(masszone)+"\n")
	f.write("\n")
	f.write("adiabaticindex="+str(adiabaticindex)+"\n")
	f.write("kpolytropicconst="+str(k[0])+"\n")
	f.write("initialvelocity="+str(u0)+"\n")
	f.write("initialdensity="+str(d[1][0]*1*10**15*mh)+"\n")
	f.write("innerpressure="+str(p0)+"\n")
	f.write("outerpressure="+str(outerpressure)+"\n")
	f.write("initialvelocity="+str(u0)+"\n")
	f.write("ntimesteps="+str(maxt)+"\n")
	f.write("nlevels="+str(nlevels)+"\n")
	f.write("initiallevelgap="+str(spacestep)+"\n")
	f.write("initialinnerposition="+str(x0)+"\n")
	f.write("initialouterposition="+str(outerposition)+"\n")
	f.close()
	np.savetxt(outputname+"positions.txt",pos,fmt="%.10e")
	np.savetxt(outputname+"velocities.txt",vel,fmt="%.10e")
	np.savetxt(outputname+"pressures.txt",pres,fmt="%.10e")
	np.savetxt(outputname+"viscosities.txt",pseu,fmt="%.10e")
	np.savetxt(outputname+"specificvolumes.txt",Vol,fmt="%.10e")
	np.savetxt(outputname+"densities.txt",den,fmt="%.10e")
	np.savetxt(outputname+"energies.txt",En,fmt="%.10e")
	np.savetxt(outputname+"timesteplog.txt",time,fmt="%.10e")
	np.savetxt(outputname+"masszone.txt",ma,fmt="%.10e")
	np.savetxt(outputname+"temps.txt",Temp,fmt="%.10e")

def main():
	innermass,E,T=setUp() #Setinitialconditions
	R,U,V,d,p,q,T=fixedBC()	#Setboundaryconditions
	#Set up experiments
	if ShockBlast==True:
		p[1][0]=shockscale*p[1][0]
	if ShockPiston==True:
		U[1][0]=pistonspeed
	
	#set timer
	t=0
	timesteplog =[0]*maxt
	# start clock
	while t<(maxt-1):
		print("Step:"+str(t+1))

		#calculate timestep - 0.5 for boundary of error
		if constantTimestep==False:
			timestep=0.5*newTimeStep(t)
		else:
			timestep=0.5*constTimestep
		timesteplog[t]=timestep

		innermass[t+1]=innermass[t]+massrate*innermass[0]*timestep 		#step wind experiment mass
		
		#Calculatenextstep
		for x in range(1,nlevels):
			stepForwardCore(x,t,timestep)
		for x in range(1,nlevels):
			stepForwardPhysics(x,t,timestep)
			T[x][t+1]=(p[x][t+1])/(d[x][t+1])
		
		t+=1
	
	print(" ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃ ̃")
	#Scaleresultstosolarvalues
	mass=innermass
	T=T*(0.5*mh/kb)#*mh**(-2/3)
	vsound0=math.sqrt((adiabaticindex*p[1][0])/d[1][0])
	#Make time record showing evolution of timestep 
	ts=[0]
	for i in range(1,maxt-1):
		ts.append(ts[i-1]+timesteplog[i])
	#Show evolution of grid-gives a good idea of what happened
	for i in range(len(R)):
		plt.plot(ts,R[i][:-1]/sr)
	plt.title("PositionofGridMidpoints")
	plt.xlabel("time(s)")
	plt.ylabel("position(solarradii)")
	plt.show()

	#give option of writing results to file
	towrite=input("Writeoutputfiles?(y/n)")
	if towrite=="y":
		outputname=input("Namethefiles:")
		writeout(outputname,R,U,p,q,V,d,E,ts,mass,T)
		print("Done")

main()
