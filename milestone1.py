from numpy import *
import matplotlib.pyplot as plt


def read_cosmology(filename):
   infile = open(filename, 'r')
   x_list = []
   eta_list = []
   Hp_list = []
   dHp_list = []
   OmegaB_list = []
   OmegaCDM_list = []
   OmegaLambda_list = []
   OmegaR_list = []
   OmegaNu_list =[]
   OmegaK_list = []
   for line in infile:
      values = line.split()
      x = float(values[0])
      eta = float(values[1])
      Hp =  float(values[2])
      dHp =  float(values[3])
      OmegaB = float(values[4])
      OmegaCDM = float(values[5])
      OmegaLambda =  float(values[6])
      OmegaR =  float(values[7])
      OmegaNu=  float(values[8])
      OmegaK =  float(values[9])
      x_list.append(x)
      eta_list.append(eta)
      Hp_list.append(Hp)
      dHp_list.append(dHp)
      OmegaB_list.append(OmegaB)
      OmegaCDM_list.append(OmegaCDM)
      OmegaLambda_list.append(OmegaLambda)
      OmegaR_list.append(OmegaR)
      OmegaNu_list.append(OmegaNu)
      OmegaK_list.append(OmegaK)

   x_array = array(x_list)
   eta_array = array(eta_list)
   Hp_array = array(Hp_list)
   dHp_array = array(dHp_list)
   OmegaB_array = array(OmegaB_list)
   OmegaCDM_array = array(OmegaCDM_list)
   OmegaLambda_array = array(OmegaLambda_list)
   OmegaR_array = array(OmegaR_list)
   OmegaNu_array = array(OmegaNu_list)
   OmegaK_array = array(OmegaK_list)
   infile.close()
   return x_array, eta_array, Hp_array, dHp_array, OmegaB_array, OmegaCDM_array, OmegaLambda_array, OmegaR_array, OmegaNu_array, OmegaK_array

x_array, eta_array, Hp_array, dHp_array, OmegaB_array, OmegaCDM_array, OmegaLambda_array, OmegaR_array, OmegaNu_array, OmegaK_array = read_cosmology("cosmology.txt")
a_array = exp(x_array)
z_array = 1./a_array - 1.0
H_array = Hp_array/a_array # unit 1/s
H_array = H_array*3.08567758e19 # unit km/s/Mpc
eta_array = eta_array/3.08567758e22 #unit Mpc

c = 2.99792458e8
#find a_eq, to mark on the plot in the milestone 4
for i in range(len(a_array)):
   if OmegaB_array[i]+OmegaCDM_array[i] > OmegaR_array[i]-1e-4 and OmegaB_array[i]+OmegaCDM_array[i] < OmegaR_array[i]+1e-4 :
      k_eq = Hp_array[i]/c;
      print i, a_array[i], x_array[i] #4673 0.00018682092878159512 -8.58536
      print k_eq #4.0680809922176227e-25


plt.figure()
plt.title('H(z)', fontsize=15)
plt.plot(z_array, H_array, color='firebrick')
plt.xscale("log")
plt.yscale("log")
plt.xlabel("z", fontsize=14)
plt.ylabel("H [km/s/Mpc]", fontsize=14)
plt.show()

plt.figure()
plt.title('H(x)', fontsize=15)
plt.plot(x_array, H_array, color='darkorchid')
plt.yscale("log")
plt.xlabel("x = ln(a)", fontsize=14)
plt.ylabel("H [km/s/Mpc]", fontsize=14)
plt.show()

plt.figure()
plt.plot(x_array, eta_array, color='teal')
plt.yscale("log")
plt.xlabel("x = ln(a)", fontsize=14)
plt.ylabel("$\eta$ [Mpc]", fontsize=14)
plt.show()

plt.figure()
plt.plot(x_array, OmegaB_array+OmegaCDM_array, label='$\Omega_b + \Omega_{CDM}$', color='dodgerblue')
plt.plot(x_array, OmegaR_array, label='$\Omega_r$', color='salmon')
plt.plot(x_array, OmegaLambda_array, label='$\Omega_{\Lambda}$', color='purple')
plt.legend(fontsize=14)
#plt.yscale("log")
plt.xlabel("x = ln(a)", fontsize=14)
plt.show()



