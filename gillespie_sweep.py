import numpy as np
import matplotlib.pyplot as plt

# 1. Parameter Sweep Configuration
decoy_loads = [5, 15, 30]  # Testing mild, moderate, and severe retroactivity load
t_max = 600.0

plt.figure(figsize=(9, 5.5))

for T2_initial in decoy_loads:
    # Time & State Initialization
    t = 0.0
    dcas9, sgrna, C = 0, 0, 0
    T1, R1 = 2, 0
    T2, R2 = T2_initial, 0
    
    time_points = [t]
    T1_points = [R1]
    
    # Kinetic Rates
    k_form = 0.002
    k_on1, k_off1 = 0.5, 0.0001
    k_on2, k_off2 = 0.05, 0.001
    k_tx_dcas9, k_tx_sgrna, k_deg = 0.5, 0.7, 0.01
    
    # Execution Loop
    while t < t_max:
        a1 = k_form * dcas9 * sgrna
        a2 = k_on1 * C * T1
        a3 = k_off1 * R1
        a4 = k_on2 * C * T2
        a5 = k_off2 * R2
        a6 = k_tx_dcas9
        a7 = k_tx_sgrna
        a8 = k_deg * dcas9
        a9 = k_deg * sgrna
        a10 = k_deg * C
        
        props = np.array([a1, a2, a3, a4, a5, a6, a7, a8, a9, a10])
        a_0 = np.sum(props)
        if a_0 <= 0:
            break
            
        r1, r2 = np.random.random(2)
        tau = -np.log(r1) / a_0
        t += tau
        
        reaction_index = np.searchsorted(np.cumsum(props / a_0), r2)
        
        if reaction_index == 0:
            dcas9 -= 1; sgrna -= 1; C += 1
        elif reaction_index == 1:
            C -= 1; T1 -= 1; R1 += 1
        elif reaction_index == 2:
            C += 1; T1 += 1; R1 -= 1
        elif reaction_index == 3:
            C -= 1; T2 -= 1; R2 += 1
        elif reaction_index == 4:
            C += 1; T2 += 1; R2 -= 1
        elif reaction_index == 5:
            dcas9 += 1
        elif reaction_index == 6:
            sgrna += 1
        elif reaction_index == 7:
            dcas9 -= 1
        elif reaction_index == 8:
            sgrna -= 1
        elif reaction_index == 9:
            C -= 1
            
        time_points.append(t)
        T1_points.append(R1)

    # Plot each load profile
    plt.step(time_points, T1_points, where='post', label=f'Decoy Load = {T2_initial} copies', linewidth=2)

# 2. Formatting the Sweep Output
plt.xlabel('Time (seconds)')
plt.ylabel('Target 1 Repression Count')
plt.title('Robustness Analysis: Circuit Performance Across Variable Decoy Loads')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()