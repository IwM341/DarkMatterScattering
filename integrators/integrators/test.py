import numpy as np
from time import process_time 
import scipy.integrate as integrate

mu = 0.01
def f(x):
    return 1.0/np.sqrt(x*x+mu*mu)

realRes = np.log((1+np.sqrt(1+mu*mu))/(2*mu))

t1_start = process_time()
Res = integrate.quad(f, 0, 1)
t1_stop = process_time()

print(Res,realRes,t1_stop-t1_start)