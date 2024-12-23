import numpy as np
import matplotlib.pyplot as plt


## Using the newton-raphson method (can use fsolve as well) to get the cdnc 
## using a fog friendly derivation (LWC != 0 at the surface)

def f_x(n,rho,gamma_eff,tau,r_e,q_0):
    return(((10*rho**(2/3)*gamma_eff*tau)/((972*np.pi*n)**(1/3)))-
            (((4*np.pi*rho*n)**(5/3)*r_e**5)/(3**(5/3)))
            +q_0**(5/3))

def f_x_prime(n,rho,gamma_eff,tau,r_e):
    return(-((10*rho**(2/3)*gamma_eff*tau)/(3*(972*np.pi)**(1/3)*n**(4/3)) + 
              (5*(4*np.pi*rho)**(5/3)*n**(2/3)*r_e**5)/(3**(8/3))))


## using some dummy constants

rho = 997
gamma_eff = .6*2*10**-6
tau = 50
r_e = 10*10**-6
q_0 = 0.05*10**-3 


n = 10**6  ## a guess in #/m^3 to stay in SI units (1 droplet/cc here)
n_list = []
index = []


for i in range(0,1000):
    n_new = n - f_x(n, rho, gamma_eff,tau,r_e,q_0)/f_x_prime(n,rho,gamma_eff,tau,r_e)
    n_list.append(n)
    index.append(i)
    if abs(n-n_new) < 0.00001:  ## setting a tolerance
        break
    n = n_new


print(n/10**6)

steps = index[-1]
plt.figure()
plt.scatter(index,np.array(n_list)*10**-6, s = 20)
plt.ylabel('CDNC (#/cm^3)')
plt.xlabel('Number of Steps')
plt.title(f'Convergance of N_eff (initial guess = {n_list[0]/10**6:.0f} cm-3)')
plt.legend([f'Number of Steps = {steps}'])

plt.savefig('cdnc_convergence')


### now let's show the CDNC for a range of surface LWC values

q_0 = np.array([0,.025,0.05,0.1,0.2,0.3,0.4,0.5])*10**(-3)


n_list = []
index = []


for val in range(len(q_0)):
    n = 10**6  ## a guess in #/m^3 to stay in SI units
    for i in range(0,1000):
        n_new = n - f_x(n, rho, gamma_eff,tau,r_e,q_0[val])/f_x_prime(n,rho,gamma_eff,tau,r_e)
        if abs(n-n_new) < 0.00001:  ## setting a tolerance
            break
        n = n_new
    n_list.append(n*10**-6)
    index.append(i)

print(n/10**6)
plt.figure()
plt.scatter(q_0*10**3,n_list)
plt.ylabel('CDNC (#/cm$^3$)')
plt.xlabel('q$_0$ (g/m^3)')
plt.title('CDNC as a Function of Surface Liquid Water Content')

plt.savefig('cdnc_as_function_of_surface_lwc')




