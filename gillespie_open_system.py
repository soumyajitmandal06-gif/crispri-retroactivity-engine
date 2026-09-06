import numpy as np
import matplotlib.pyplot as plt

# 1. Time parameters (extended to watch steady-state formation)
t = 0.0
t_max = 1000.0

# 2. Binding Kinetics (The Rescued Circuit)
k_form = 0.002
k_on1 = 0.5      # Target 1 (Strong)
k_off1 = 0.0001
k_on2 = 0.05     # Target 2 (Weak Decoy)
k_off2 = 0.001

# 3. Dynamic Cell Kinetics (Birth and Death)
k_tx_dcas9 = 0.5  # New dCas9 produced per second
k_tx_sgrna = 0.7  # New sgRNA produced per second
k_deg = 0.01      # Dilution rate due to cell growth

# 4. Initial Molecular Counts (Starting from a completely empty cell)
dcas9, sgrna, C = 0, 0, 0
T1, R1 = 2, 0
T2, R2 = 15, 0

# 5. Data Storage
time_points = [t]
T1_points = [R1]
T2_points = [R2]
C_points = [C]
dcas9_points = [dcas9] # Track free dCas9 production
# 6. The Open-System Gillespie Engine
np.random.seed(42)
while t < t_max:
    # 1-5: Circuit Binding Kinetics
    a1 = k_form * dcas9 * sgrna       # Form complex
    a2 = k_on1 * C * T1               # Bind Target 1
    a3 = k_off1 * R1                  # Unbind Target 1
    a4 = k_on2 * C * T2               # Bind Target 2
    a5 = k_off2 * R2                  # Unbind Target 2
    
    # 6-10: Cell Birth & Death Kinetics
    a6 = k_tx_dcas9                   # Produce dCas9 (Zero-order)
    a7 = k_tx_sgrna                   # Produce sgRNA (Zero-order)
    a8 = k_deg * dcas9                # Degrade free dCas9 (First-order)
    a9 = k_deg * sgrna                # Degrade free sgRNA (First-order)
    a10 = k_deg * C                   # Degrade free Complex (First-order)
    
    props = np.array([a1, a2, a3, a4, a5, a6, a7, a8, a9, a10])
    a_0 = np.sum(props)
    
    if a_0 <= 0:
        break
        
    r1, r2 = np.random.random(2)
    tau = -np.log(r1) / a_0
    t += tau
    
    reaction_probs = props / a_0
    cumulative_probs = np.cumsum(reaction_probs)
    reaction_index = np.searchsorted(cumulative_probs, r2)
    
    # Execute winning reaction
    if reaction_index == 0:    # Form complex
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
    elif reaction_index == 5:  # Produce dCas9
        dcas9 += 1
    elif reaction_index == 6:  # Produce sgRNA
        sgrna += 1
    elif reaction_index == 7:  # Degrade dCas9
        dcas9 -= 1
    elif reaction_index == 8:  # Degrade sgRNA
        sgrna -= 1
    elif reaction_index == 9:  # Degrade Complex
        C -= 1
        
    # Record states
    time_points.append(t)
    T1_points.append(R1)
    T2_points.append(R2)
    C_points.append(C)
    dcas9_points.append(dcas9)

# 7. Plotting the Open System Dynamics (OUTSIDE THE LOOP)
plt.figure(figsize=(9, 5.5))
plt.step(time_points, T1_points, where='post', label='Repressed Target 1', color='#00ff00', linewidth=2.5)
plt.step(time_points, T2_points, where='post', label='Repressed Target 2 (Decoy)', color='#ff0000', linewidth=1.5, alpha=0.7)
plt.step(time_points, C_points, where='post', label='Free Complexes', color='blue', linestyle='--', linewidth=1.5)
plt.step(time_points, dcas9_points, where='post', label='Free dCas9 Pool', color='purple', alpha=0.5)

plt.xlabel('Time (seconds)')
plt.ylabel('Discrete Molecular Counts')
plt.title('Dynamic Open System: CRISPRi Circuit with Cell Growth & Dilution')
plt.legend(loc='upper left')
plt.grid(True, alpha=0.3)
plt.show()