
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
path="C:/Users/Mette Wirak/Documents/Faglig/Universitetet i Oslo/Computational Physics/Project4/" #Endre til din

filenames=["4e 40 2.0 2.3 10_7.txt","4e 60 2.0 2.3 10_7.txt", "4e 80 2.0 2.3 10_7.txt", "4e 100 2.0 2.3 10_7.txt" ]

fil_label=["$L=40, n=10^6$", "$L=40, n=10^7$", "$L=60, n=10^7$", "$L=80, n=10^7$"]

i=0



f, [E_exp,E_var] = plt.subplots(2, sharex=True)
f, [M_exp,M_var, M_abs] = plt.subplots(3, sharex=True)

M_abs.set_xlabel('Temperature')
M_exp.set_ylabel('Magnetic Expectation Values')
E_var.set_xlabel('Temperature')
E_exp.set_ylabel('Energy Expectation Values')
M_abs.set_ylabel('Magnetic absolut Values')
#ax.set_ylabel('common ylabel')


for fil in filenames:
	fil=path+fil
	temperature,E_ExpectationValues,Evariance,M_ExpectationValues, Mvariance,Mabs_ExpectationValues=read(fil)
	E_exp.plot(temperature, E_ExpectationValues, label=fil_label[i])
	E_var.plot(temperature, Evariance, label=fil_label[i])	
	M_var.plot(temperature, Mvariance, label=fil_label[i])
	M_exp.plot(temperature, M_ExpectationValues, label=fil_label[i])	
	M_abs.plot(temperature, Mabs_ExpectationValues, label=fil_label[i])

	i=+1
#evar.xlabel('T')
#mabs.xlabel('T')
#mabs.ylabel('absoluttvardien av M')
plt.legend(loc=2)
plt.show()
