#options file where all changable model parameters are set

#Physical constants
G =6.67*10**(-11)		#gravity
sr =10					#solar radius
sm=1e13					#solar mass
mh=1.67*10**(-27)		#proton mass
kb=1.38*10**(-23)		#boltzman constant
adiabaticindex=5/3

#model run settings
maxt=5000				#n timesteps
constantTimestep=True	#does the timestep stay the same length or adjust to the circumstance, if F see newTimeStep()
constTimestep=2			#length of timestep if constant
useEnergyEq=False		#use alternative energy equation in the inner zone
zeroOuterPGrad =False	#set 0 pressure gradient between space and outer zone

#Spatial initial conditions
x0=sr*1.1				#position of innermost layer
outerposition=100*sr*5	#position of outermost layer
nlevels=10				#number of layers
k=[1]*(nlevels+1)		#constant at each layer
spacestep=(outerposition-x0)/nlevels 


#boundary conditions
d0=1					#innermost density
outerdensity=1e-10		#outermost density
p0=k[0]*d0**(adiabaticindex) #innermost pressure set by density relation
pinner=p0
outerpressure=k[0]*outerdensity**(adiabaticindex) #outermost pressure same
#other physical options
a=2						#pseudo viscosity (=0 to turn off)
gravity=True			#should gravity exist
u0=0					#do zones start in motion

#Experiment options
ShockBlast=False 		# A one off high pressure perturbation of desired magnitude
shockscale=40
ShockPiston=False       #The 1st zone moves upwards at the desired speed
pistonspeed=1

#Windoptions-specify percentage of additional quantity to be added eachs econd and grid point to add it in.
massrate=0#.001
masszone=1
pressurerate=5e-3
pressurezone=1
energyrate=0
energyzone=1

