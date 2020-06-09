from numpy import *
import matplotlib.pyplot as plt


def read_C_ell(filename):
   infile = open(filename, 'r')
   ell_list = []
   C_list = []
   for line in infile:
      values = line.split()
      ell = float(values[0])
      C = float(values[1])
      ell_list.append(ell)
      C_list.append(C)
   ell_array = array(ell_list)
   C_array = array(C_list)
  
   infile.close()
   return ell_array, C_array

ell, C_ell = read_C_ell("cells.txt")


def read_k(filename):
   infile = open(filename, 'r')
   k_list = []
   Pk_list = []
   transfer_f_list_ell1 = []
   transfer_f_list_ell2 = []
   transfer_f_list_ell3 = []
   integrand_list_ell1 = []
   integrand_list_ell2 = []
   integrand_list_ell3 = []
   for line in infile:
      values = line.split()
      k = float(values[0])
      Pk = float(values[1])
      transfer_f_ell1 =  float(values[2])
      transfer_f_ell2 =  float(values[3])
      transfer_f_ell3 =  float(values[4])
      integrand_ell1 =  float(values[5])
      integrand_ell2 =  float(values[6])
      integrand_ell3 =  float(values[7])
      k_list.append(k)
      Pk_list.append(Pk)
      transfer_f_list_ell1.append(transfer_f_ell1)
      transfer_f_list_ell2.append(transfer_f_ell2)
      transfer_f_list_ell3.append(transfer_f_ell3)
      integrand_list_ell1.append(integrand_ell1)
      integrand_list_ell2.append(integrand_ell2)
      integrand_list_ell3.append(integrand_ell3)
   k_array = array(k_list)
   Pk_array = array(Pk_list)
   transfer_f_array_ell1 = array(transfer_f_list_ell1)
   transfer_f_array_ell2 = array(transfer_f_list_ell2)
   transfer_f_array_ell3 = array(transfer_f_list_ell3)
   integrand_array_ell1 = array(integrand_list_ell1)
   integrand_array_ell2 = array(integrand_list_ell2)
   integrand_array_ell3 = array(integrand_list_ell3)

   infile.close()
   return k_array, Pk_array, transfer_f_array_ell1, transfer_f_array_ell2, transfer_f_array_ell3, integrand_array_ell1, integrand_array_ell2, integrand_array_ell3

k, Pk, transfer_f_ell1, transfer_f_ell2, transfer_f_ell3, integrand_ell1, integrand_ell2, integrand_ell3  = read_k("ks.txt")

k_eq_Mpc = 4.0680809922176227e-25*3.08567758e22

plt.figure()
plt.plot(ell, C_ell, color="m")
plt.xscale("log")
plt.yscale("log")
plt.ylabel("$C_{\ell}\ell(\ell+1)/2 \pi$ $[\mu K^2]$", fontsize=14)
plt.xlabel("$\ell$", fontsize=15)
plt.show()

plt.figure()
plt.plot(k, Pk, color="firebrick")
plt.axvline(k_eq_Mpc,0, 1, color="lightcoral")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("k [h/Mpc]", fontsize=14)
plt.ylabel("P(k) [(Mpc/h)$^3$]", fontsize=14)
plt.show()

plt.figure()
plt.plot(k, transfer_f_ell1, label="$\ell = 20$")
plt.plot(k, transfer_f_ell2,label="$\ell = 500$")
plt.plot(k, transfer_f_ell3,label="$\ell = 1200$")
plt.xscale("log")
plt.yscale("log")
plt.ylabel("$\Theta_{\ell}(k)$", fontsize=14)
plt.xlabel("k [h/Mpc]", fontsize=14)
plt.legend(fontsize=14)
plt.show()

plt.figure()
plt.plot(k, integrand_ell1,label="$\ell = 20$")
plt.plot(k, integrand_ell2,label="$\ell = 500$")
plt.plot(k, integrand_ell3,label="$\ell = 1200$")
plt.xscale("log")
plt.yscale("log")
plt.legend(fontsize=14)
plt.ylabel("$\Theta_{\ell}(k)^2/k$", fontsize=14)
plt.xlabel("k [h/Mpc]", fontsize=14)
plt.show()






