from numpy import *
import matplotlib.pyplot as plt


def read_perturbations(filename):
   infile = open(filename, 'r')
   x_list = []
   delta_cdm_list = []
   delta_b_list = []
   v_cdm_list = []
   v_b_list = []
   Phi_list = []
   Theta0_list = []
   Theta1_list = []
   Psi_list = []
   
   
   for line in infile:
      values = line.split()
      x = float(values[0])
      delta_cdm = float(values[1])
      delta_b =  float(values[2])
      v_cdm =  float(values[3])
      v_b = float(values[4])
      Phi = float(values[5])
      Theta0 =  float(values[6])
      Theta1 =  float(values[7])
      Psi =  float(values[8])
      
      x_list.append(x)
      delta_cdm_list.append(delta_cdm)
      delta_b_list.append(delta_b)
      v_cdm_list.append(v_cdm)
      v_b_list.append(v_b)
      Phi_list.append(Phi)
      Theta0_list.append(Theta0)
      Theta1_list.append(Theta1)
      Psi_list.append(Psi)
      

   x_array = array(x_list)
   delta_cdm_array = array(delta_cdm_list)
   delta_b_array = array(delta_b_list)
   v_cdm_array = array(v_cdm_list)
   v_b_array = array(v_b_list)
   Phi_array = array(Phi_list)
   Theta0_array = array(Theta0_list)
   Theta1_array = array(Theta1_list)
   Psi_array = array(Psi_list)
   
   infile.close()
   return x_array, delta_cdm_array, delta_b_array, v_cdm_array, v_b_array, Phi_array, Theta0_array, Theta1_array, Psi_array

#k1 = 0.01/Mpc
x_1, delta_cdm_1, delta_b_1, v_cdm_1, v_b_1, Phi_1, Theta0_1, Theta1_1, Psi_1 = read_perturbations("perturbations_k0.01.txt")

#k2 = 0.001/Mpc
x_2, delta_cdm_2, delta_b_2, v_cdm_2, v_b_2, Phi_2, Theta0_2, Theta1_2, Psi_2 = read_perturbations("perturbations_k0.001.txt")

#k3 = 0.1/Mpc
x_3, delta_cdm_3, delta_b_3, v_cdm_3, v_b_3, Phi_3, Theta0_3, Theta1_3, Psi_3 = read_perturbations("perturbations_k0.1.txt")

#k4 = 0.0005/Mpc
x_4, delta_cdm_4, delta_b_4, v_cdm_4, v_b_4, Phi_4, Theta0_4, Theta1_4, Psi_4 = read_perturbations("perturbations_k0.0005.txt")

#k5 = 0.005/Mpc
x_5, delta_cdm_5, delta_b_5, v_cdm_5, v_b_5, Phi_5, Theta0_5, Theta1_5, Psi_5 = read_perturbations("perturbations_k0.005.txt")

#k6 = 0.15/Mpc
x_6, delta_cdm_6, delta_b_6, v_cdm_6, v_b_6, Phi_6, Theta0_6, Theta1_6, Psi_6 = read_perturbations("perturbations_k0.15.txt")

plt.figure()
plt.title('$\\delta_{cdm}$', fontsize=15)
plt.plot(x_4, delta_cdm_4, label='k = 0.0005/Mpc', color='darksalmon')
plt.plot(x_2, delta_cdm_2, label='k = 0.001/Mpc', color='firebrick')
plt.plot(x_5, delta_cdm_5, label='k = 0.005/Mpc', color='mediumaquamarine')
plt.plot(x_1, delta_cdm_1, label='k = 0.01/Mpc', color='midnightblue')
plt.plot(x_3, delta_cdm_3, label='k = 0.1/Mpc', color='orchid')
plt.plot(x_6, delta_cdm_6, label='k = 0.15/Mpc', color='purple')
plt.yscale("log")
plt.legend(fontsize=14)
plt.xlabel("x = ln(a)", fontsize=14)
plt.show()

plt.figure()
plt.title('$|\\delta_{b}|$', fontsize=15)
plt.plot(x_4, abs(delta_b_4), label='k = 0.0005/Mpc', color='darksalmon')
plt.plot(x_2, abs(delta_b_2), label='k = 0.001/Mpc', color='firebrick')
plt.plot(x_5, abs(delta_b_5), label='k = 0.005/Mpc', color='mediumaquamarine')
plt.plot(x_1, abs(delta_b_1), label='k = 0.01/Mpc', color='midnightblue')
plt.plot(x_3, abs(delta_b_3), label='k = 0.1/Mpc', color='orchid')
plt.plot(x_6, abs(delta_b_6), label='k = 0.15/Mpc', color='purple')
plt.yscale("log")
#plt.legend(fontsize=14)
plt.xlabel("x = ln(a)", fontsize=14)
plt.show()

plt.figure()
plt.title('$v_{cdm}$', fontsize=15)
plt.plot(x_4, v_cdm_4, label='k = 0.0005/Mpc', color='darksalmon')
plt.plot(x_2, v_cdm_2, label='k = 0.001/Mpc', color='firebrick')
plt.plot(x_5, v_cdm_5, label='k = 0.005/Mpc', color='mediumaquamarine')
plt.plot(x_1, v_cdm_1, label='k = 0.01/Mpc', color='midnightblue')
plt.plot(x_3, v_cdm_3, label='k = 0.1/Mpc', color='orchid')
plt.plot(x_6, v_cdm_6, label='k = 0.15/Mpc', color='purple')
plt.yscale("log")
plt.legend(fontsize=14)
plt.xlabel("x = ln(a)", fontsize=14)
plt.show()

plt.figure()
plt.title('$|v_{b}|$', fontsize=15)
plt.plot(x_4, abs(v_b_4), label='k = 0.0005/Mpc', color='darksalmon')
plt.plot(x_2, abs(v_b_2), label='k = 0.001/Mpc', color='firebrick')
plt.plot(x_5, abs(v_b_5), label='k = 0.005/Mpc', color='mediumaquamarine')
plt.plot(x_1, abs(v_b_1), label='k = 0.01/Mpc', color='midnightblue')
plt.plot(x_3, abs(v_b_3), label='k = 0.1/Mpc', color='orchid')
plt.plot(x_6, abs(v_b_6), label='k = 0.15/Mpc', color='purple')
plt.yscale("log")
#plt.legend(fontsize=14)
plt.xlabel("x = ln(a)", fontsize=14)
plt.show()

plt.figure()
plt.title('$\\Phi$', fontsize=15)
plt.plot(x_2, Phi_4, label='k = 0.0005/Mpc', color='darksalmon')
plt.plot(x_2, Phi_2, label='k = 0.001/Mpc', color='firebrick')
plt.plot(x_5, Phi_5, label='k = 0.005/Mpc', color='mediumaquamarine')
plt.plot(x_1, Phi_1, label='k = 0.01/Mpc', color='midnightblue')
plt.plot(x_3, Phi_3, label='k = 0.1/Mpc', color='orchid')
plt.plot(x_6, Phi_6, label='k = 0.15/Mpc', color='purple')
plt.legend(fontsize=13)
plt.xlabel("x = ln(a)", fontsize=14)
plt.show()

plt.figure()
plt.title('$\\delta_{\\gamma} = 4\\Theta_{0}$', fontsize=15)
plt.plot(x_4, 4.*Theta0_4, label='k = 0.0005/Mpc', color='darksalmon')
plt.plot(x_2, 4.*Theta0_2, label='k = 0.001/Mpc', color='firebrick')
plt.plot(x_5, 4.*Theta0_5, label='k = 0.005/Mpc', color='mediumaquamarine')
plt.plot(x_1, 4.*Theta0_1, label='k = 0.01/Mpc', color='midnightblue')
plt.plot(x_3, 4.*Theta0_3, label='k = 0.1/Mpc', color='orchid')
plt.plot(x_6, 4.*Theta0_6, label='k = 0.15/Mpc', color='purple')
plt.legend(fontsize=14)
plt.xlabel("x = ln(a)", fontsize=14)
plt.show()

plt.figure()
plt.title('$v_{\\gamma} = -3\\Theta_{1}$', fontsize=15)
plt.plot(x_4, -3.*Theta1_4, label='k = 0.0005/Mpc', color='darksalmon')
plt.plot(x_2, -3.*Theta1_2, label='k = 0.001/Mpc', color='firebrick')
plt.plot(x_5, -3.*Theta1_5, label='k = 0.005/Mpc', color='mediumaquamarine')
plt.plot(x_1, -3.*Theta1_1, label='k = 0.01/Mpc', color='midnightblue')
plt.plot(x_3, -3.*Theta1_3, label='k = 0.1/Mpc', color='orchid')
plt.plot(x_6, -3.*Theta1_6, label='k = 0.15/Mpc', color='purple')
plt.legend(fontsize=14)
plt.xlabel("x = ln(a)", fontsize=14)
plt.show()

plt.figure()
plt.title('$\\Psi$', fontsize=15)
plt.plot(x_4, Psi_4, label='k = 0.0005/Mpc', color='darksalmon')
plt.plot(x_2, Psi_2, label='k = 0.001/Mpc', color='firebrick')
plt.plot(x_5, Psi_5, label='k = 0.005/Mpc', color='mediumaquamarine')
plt.plot(x_1, Psi_1, label='k = 0.01/Mpc', color='midnightblue')
plt.plot(x_3, Psi_3, label='k = 0.1/Mpc', color='orchid')
plt.plot(x_6, Psi_6, label='k = 0.15/Mpc', color='purple')
#plt.legend(fontsize=14)
plt.xlabel("x = ln(a)", fontsize=14)
plt.show()

"""
plt.figure()
plt.title('$\\Psi + \\Phi$', fontsize=15)
plt.plot(x_2, Psi_2+Phi_2, label='k = 0.001/Mpc')
plt.plot(x_1, Psi_1+Phi_1, label='k = 0.01/Mpc')
plt.plot(x_3, Psi_3+Phi_3, label='k = 0.1/Mpc')
plt.plot(x_3, Phi_3, label='Phi, k = 0.1/Mpc')
plt.legend(fontsize=14)
#plt.axis([-16, 0, 10**(-14), 10**3])
plt.yscale("log")
plt.show()
"""
