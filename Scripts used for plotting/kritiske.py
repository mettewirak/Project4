
import numpy as np
import matplotlib.pyplot as plt


L=np.array([40.0,60.0,80.0,100.0])
T_val=np.array([2.85,2.8,2.75,2.7])

T=T_val

fit= np.polyfit(1/L,T,1)
print fit
L_fit=np.linspace(30,500,1000)
plt.plot(1/L,T, 'o',label='Simulated values')
plt.plot(1/L_fit,fit[0]*1/L_fit+fit[1], label='Best fit')
plt.xlabel("1/N")
plt.ylabel("$T_C$")
plt.legend()
plt.savefig("4f")
plt.show()