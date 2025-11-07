import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt('T_loop.dat')

# Extract columns of data
Temperature = data[:, 0]
average_magnetisation = data[:, 1]
average_energy = data[:, 2]
cv = data[:, 3]
chi = data[:, 4]
avg_abs_magnetisation = data[:, 5]
binder = data[:,6]

# Plot Monte Carlo steps vs. average energy
plt.figure()
plt.plot(Temperature, average_energy,linestyle='-', color='b')
plt.title('Temperature vs. Average Energy')
plt.xlabel('Temperature')
plt.ylabel('Average Energy')
plt.grid(True)
plt.savefig('Temperature_vs_average_energy1.png')
plt.close()

# Plot Monte Carlo steps vs. average magnetisation
plt.figure()
plt.plot(Temperature, average_magnetisation, linestyle='-', color='g')
plt.title('Temperatures vs. Average Magnetisation')
plt.xlabel('Temperature')
plt.ylabel('Average Magnetisation')
plt.grid(True)
plt.savefig('Temperature_vs_average_magnetisation1.png')
plt.close()

# Plot Monte Carlo steps vs. Cv
plt.figure()
plt.plot(Temperature, cv,  linestyle='-', color='r')
plt.title('Temperature vs. Cv')
plt.xlabel('Temperature')
plt.ylabel('Cv')
plt.grid(True)
plt.savefig('Cv1.png')
plt.close()

# Plot Monte Carlo steps vs. Chi
plt.figure()
plt.plot(Temperature, chi, linestyle='-', color='m')
plt.title('Temperature vs. Chi')
plt.xlabel('Temperature')
plt.ylabel('Chi')
plt.grid(True)
plt.savefig('chi1.png')
plt.close()

# Plot Monte Carlo steps vs. Magnetisation Fluctuations
plt.figure()
plt.plot(Temperature, avg_abs_magnetisation, linestyle='-', color='y')
plt.title('Temperature vs Avgerage Absolute Magnetisation')
plt.xlabel('Temperature')
plt.ylabel('Average Absolute Magnetisation')
plt.grid(True)
plt.savefig('Average_absolute_magnetisation.png')
plt.close()

# Plot Monte Carlo steps vs. Binder Cumulant
plt.figure()
plt.plot(Temperature, binder, linestyle='-', color='c')
plt.title('Temperature vs. Binder Cumulant')
plt.xlabel('Temperature')
plt.ylabel('Binder Cumulant')
plt.grid(True)
plt.savefig('Binder_cumulant.png')
plt.close()

print("Plots saved as PNG files in the current folder.")

