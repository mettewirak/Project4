import numpy as np
import matplotlib.pyplot as plt


def read(filename):
	with open(filename) as f:
		lines = f.readlines()
		temperature = [line.split()[0] for line in lines]
		E_ExpectationValues = [line.split()[1] for line in lines]
		Evariance = [line.split()[2] for line in lines]
		M_ExpectationValues = [line.split()[3] for line in lines]
		Mvariance = [line.split()[4] for line in lines]
		Mabs_ExpectationValues = [line.split()[5] for line in lines]
	return temperature,E_ExpectationValues,Evariance,M_ExpectationValues, Mvariance,Mabs_ExpectationValues


fil1="build-arnlaug-Desktop-Debug/Resultater_T_1.txt"
fil2="build-arnlaug-Desktop-Debug/Resultater_T_2_4.txt"
fil3="build-arnlaug-Desktop-Debug/Resultater_T_1_opp.txt"
fil4="build-arnlaug-Desktop-Debug/Resultater_T_2_4_opp.txt"

temperature,E_ExpectationValues,Evariance,M_ExpectationValues, Mvariance,Mabs_ExpectationValues=read(fil1)
temperature_2,E_ExpectationValues_2,Evariance_2,M_ExpectationValues_2, Mvariance_2,Mabs_ExpectationValues_2=read(fil2)
temperature_3,E_ExpectationValues_3,Evariance_3,M_ExpectationValues_3, Mvariance_3,Mabs_ExpectationValues_3=read(fil3)
temperature_4,E_ExpectationValues_4,Evariance_4,M_ExpectationValues_4, Mvariance_4,Mabs_ExpectationValues_4=read(fil4)

n1=np.linspace(0,len(E_ExpectationValues),len(E_ExpectationValues))
n2=np.linspace(0,len(E_ExpectationValues_2),len(E_ExpectationValues_2))
n4=np.linspace(0,len(E_ExpectationValues_4),len(E_ExpectationValues_4))


Energy, [T1,T2] = plt.subplots(2, sharex=True)
Magnetisk, [T1M, T2M] = plt.subplots(2, sharex=True)

T1.plot(n1[10:],E_ExpectationValues[10:],label='random spin orientation ')
T1.plot(n1[10:],E_ExpectationValues_3[10:],label='orderd spin orientation')
T2.plot(n1[10:],E_ExpectationValues_2[10:],label='random spin orientation ')
T2.plot(n4[10:],E_ExpectationValues_4[10:],label='orderd spin orientation')
T2.set_xlabel("Number of MonteCarlo cycles")
T1.legend()
T1.set_title("T=1,Mean Energy")
T2.set_title("T=2.4,Mean Energy")
T1.set_ylim([-1.9955,-1.9981])
T2.set_ylim([-1.2,-1.28])
T2.set_xlabel("Number of MonteCarlo cycles")

Energy.savefig("Energy_MC")


T1M.plot(n1[10:],Mabs_ExpectationValues[10:],label='random spin orientation ')
T1M.plot(n1[10:],Mabs_ExpectationValues_3[10:],label='orderd spin orientation')
T2M.plot(n1[10:],Mabs_ExpectationValues_2[10:],label='random spin orientation ')
T2M.plot(n4[10:],Mabs_ExpectationValues_4[10:],label='orderd spin orientation')
T1M.legend()
T1M.set_title("T=1,Mean absulutvalue of the magnetic moment")
T2M.set_title("T=2.4,Mean absulutvalue of the magnetic moment")
T1M.set_ylim([0.99,1.005])
T2M.set_ylim([0.42,0.5])
T2M.set_xlabel("Number of MonteCarlo cycles")

Magnetisk.savefig("Magnetic_MC")

plt.show()


"""
plt.plot(n1[10:],E_ExpectationValues[10:],label='T=1, random')
plt.plot(n2[10:],E_ExpectationValues_2[10:],label='T=2.4, random')
plt.plot(n1[10:],E_ExpectationValues_3[10:],label='T=1, opp')
plt.plot(n2[10:],E_ExpectationValues_4[10:],label='T=2.4, opp')
"""
"""
plt.plot(n1[10:],Mabs_ExpectationValues[10:],label='T=1, random, magnetisk')
plt.plot(n1[10:],Mabs_ExpectationValues_3[10:],label='T=1, opp, magnetisk')
#plt.plot(n1[10:],Evariance[10:],label='T=1, random, Evariance')
#plt.plot(n1[10:],Evariance_3[10:],label='T=1, opp, Evariance')
plt.legend(loc=5)
plt.xlabel("Number of MC-cylcles")
plt.ylabel("Energy")
#plt.savefig("3d")
plt.show()
"""