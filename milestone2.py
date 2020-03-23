from numpy import *
import matplotlib.pyplot as plt


def read_recombination(filename):
   infile = open(filename, 'r')
   x_list = []
   Xe_list = []
   ne_list = []
   tau_list = []
   dtau_list = []
   ddtau_list = []
   g_list = []
   dg_list = []
   ddg_list =[]
   
   for line in infile:
      values = line.split()
      x = float(values[0])
      Xe = float(values[1])
      ne =  float(values[2])
      tau =  float(values[3])
      dtau = float(values[4])
      ddtau = float(values[5])
      g =  float(values[6])
      dg =  float(values[7])
      ddg =  float(values[8])
      
      x_list.append(x)
      Xe_list.append(Xe)
      ne_list.append(ne)
      tau_list.append(tau)
      dtau_list.append(dtau)
      ddtau_list.append(ddtau)
      g_list.append(g)
      dg_list.append(dg)
      ddg_list.append(ddg)
      

   x_array = array(x_list)
   Xe_array = array(Xe_list)
   ne_array = array(ne_list)
   tau_array = array(tau_list)
   dtau_array = array(dtau_list)
   ddtau_array = array(ddtau_list)
   g_array = array(g_list)
   dg_array = array(dg_list)
   ddg_array = array(ddg_list)
   
   infile.close()
   return x_array, Xe_array, ne_array, tau_array, dtau_array, ddtau_array, g_array, dg_array, ddg_array

x_array, Xe_array, ne_array, tau_array, dtau_array, ddtau_array, g_array, dg_array, ddg_array = read_recombination("recombination.txt")

a_array = exp(x_array)
z_array = 1./a_array - 1.0


for i in range(len(tau_array)):
   if tau_array[i] > 0.9 and tau_array[i] < 1.01: #found two tau's in this range, took the one closer to 1
      last_scattering_x = x_array[i]
      last_scattering_z = z_array[i]
      print "Last scattering:"
      print tau_array[i], last_scattering_x, last_scattering_z #1.00763 -6.98684 1081.2960114059683

for i in range(len(Xe_array)):
   if Xe_array[i] > 0.49 and Xe_array[i] < 0.51: 
      half_rec_x = x_array[i]
      half_rec_z = z_array[i]
      print "Recombination half-way:"
      print Xe_array[i], half_rec_x, half_rec_z #0.496914 -7.1642 1291.327325285989

x_Saha_halfway = -7.23077 #read from terminal after running C++ script
z_Saha_halfway = 1./exp(x_Saha_halfway) - 1.0
print "Saha prediction half-way:"
print x_Saha_halfway, z_Saha_halfway #-7.23077 1380.2856846872003 at Xe = 0.502406

plt.figure()
#plt.title('$X_e (x) $', fontsize=15)
plt.plot(x_array, Xe_array, color='firebrick')
plt.yscale("log")
plt.xlabel("x = ln(a)", fontsize=14)
plt.ylabel('$X_e (x) $', fontsize=14)
plt.axis([-16, 0, 1e-4, 1.5])
plt.show()


plt.figure()
plt.plot(x_array, tau_array, label='$\\tau (x)$', color='dodgerblue')
plt.plot(x_array, -dtau_array, label='$-\\tau \' (x)$', color='salmon')
plt.plot(x_array, ddtau_array,label='$\\tau \'\' (x)$', color='purple')
plt.legend(fontsize=14)
plt.yscale("log")
plt.axis([-16, 0, min(tau_array), max(tau_array)])
plt.xlabel("x = ln(a)", fontsize=14)
plt.show()

plt.figure()
plt.plot(x_array, g_array/max(abs(g_array)), label=r'$\~ g (x)$', color='dodgerblue')
plt.plot(x_array, dg_array/max(abs(dg_array)), label='$\~ g\' (x)$', color='salmon')
plt.plot(x_array, ddg_array/max(abs(ddg_array)),label='$\~ g \'\' (x)$', color='purple',  linestyle='dashed')
plt.legend(fontsize=14)
#plt.yscale("log")
plt.axis([-16, 0, -1.1, 1.1])
plt.xlabel("x = ln(a)", fontsize=14)
plt.show()

print max(abs(g_array)), max(abs(dg_array)), max(abs(ddg_array)) #4.85692 50.83 998.799

