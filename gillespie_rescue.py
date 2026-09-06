import numpy as np
import matplotlib.pyplot as plt

# 1. Time parameters
t = 0.0
t_max = 300.0

# 2. Mutated Kinetic Rates (The Rescue)
k_form = 0.002
k_on1 = 0.5      # Mutated Target 1: 10x stronger binding affinity
k_off1 = 0.0001  # Mutated Target 1: 10x slower unbinding
k_on2 = 0.05     # Target 2 (Decoy): Original binding
k_off2 = 0.001   # Target 2 (Decoy): Original unbinding

# 3. Initial Molecular Counts
dcas9, sgrna, C = 10, 15, 0
T1, R1 = 2, 0
T2, R2 = 15, 0  # The 15-copy decoy load

# 4. Data Storage
time_points = [t]
T1_points = [R1]
T2_points = [R2]
C_points = [C]

# 5. The Multi-Reaction Gillespie Engine
np.random.seed(42) # Set seed for reproducible stochastic run
while t < t_max:
    # Calculate individual propensities for all 5 reactions
    a1 = k_form * dcas9 * sgrna       # Form complex
    a2 = k_on1 * C * T1               # Bind Target 1 (Strong)
    a3 = k_off1 * R1                  # Unbind Target 1 (Slow)
    a4 = k_on2 * C * T2               # Bind Target 2 (Weak Decoy)
    a5 = k_off2 * R2                  # Unbind Target 2 (Fast)
    
    props = np.array([a1, a2, a3, a4, a5])
    a_0 = np.sum(props)
    
    # If system stalls (no reactions possible), break
    if a_0 <= 0:
        break
        
    # Generate two random floats
    r1, r2 = np.random.random(2)
    
    # Advance the biological clock
    tau = -np.log(r1) / a_0
    t += tau
    
    # Determine which reaction fired using probability partitioning
    reaction_probs = props / a_0
    cumulative_probs = np.cumsum(reaction_probs)
    reaction_index = np.searchsorted(cumulative_probs, r2)
    
    # Execute the winning reaction
    if reaction_index == 0:    # Formation
        dcas9 -= 1
        sgrna -= 1
        C += 1
    elif reaction_index == 1:  # Bind T1
        C -= 1
        T1 -= 1
        R1 += 1
    elif reaction_index == 2:  # Unbind T1
        C += 1
        T1 += 1
        R1 -= 1
    elif reaction_index == 3:  # Bind T2
        C -= 1
        T2 -= 1
        R2 += 1
    elif reaction_index == 4:  # Unbind T2
        C += 1
        T2 += 1
        R2 -= 1
        
    # Record the new discrete state
    time_points.append(t)
    T1_points.append(R1)
    T2_points.append(R2)
    C_points.append(C)

# 6. Plotting the Stochastic Rescue
plt.figure(figsize=(8, 5))
plt.step(time_points, T1_points, where='post', label='Repressed Target 1 (Strong Affinity)', color='#00ff00', linewidth=2.5)
plt.step(time_points, T2_points, where='post', label='Repressed Target 2 (Decoy Load)', color='#ff0000', linewidth=1.5, alpha=0.7)
plt.step(time_points, C_points, where='post', label='Free Complexes', color='blue', linestyle='--', linewidth=1.5)

plt.xlabel('Time (seconds)')
plt.ylabel('Discrete Molecular Counts')
plt.title('Stochastic Rescue: Overcoming Retroactivity via Affinity Tuning')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()