import numpy as np

# Time parameters
t = 0.0
t_max = 100.0

# Rate constant
k_bind = 0.002

# Initial molecular counts
dcas9 = 10
sgrna = 15
complexes = 0
# Data storage for plotting
time_points = [t]
complexes_points = [complexes]
# The Gillespie Engine
while t < t_max:
    # 1. Calculate dynamic propensity for this exact microsecond
    prop_bind = k_bind * dcas9 * sgrna
    a_0 = prop_bind
    
    # If no molecules are left to react, stop the simulation
    if a_0 == 0:
        break
        
    # 2. Generate random floats
    r1, r2 = np.random.random(2)
    
    # 3. Advance the biological clock
    tau = -np.log(r1) / a_0
    t += tau
    # 4. Execute the reaction and update physical counts
    dcas9 -= 1
    sgrna -= 1
    complexes += 1
    
    # 5. Record the new state for plotting
    time_points.append(t)
    complexes_points.append(complexes)

# 6. Plot the stochastic noise (OUTSIDE THE LOOP)
import matplotlib.pyplot as plt

plt.step(time_points, complexes_points, where='post', color='blue', linewidth=2)
plt.xlabel('Time (seconds)')
plt.ylabel('dCas9-sgRNA Complexes')
plt.title('Stochastic Gillespie Simulation of CRISPRi Complex Formation')
plt.grid(True, alpha=0.3)
plt.show()