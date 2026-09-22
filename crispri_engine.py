import numpy as np
import matplotlib.pyplot as plt

class CRISPRiEngine:
    def __init__(self, decoy_load=15, t_max=3600.0):
        """Initializes the engine with biologically accurate decoupled degradation."""
        self.t_max = t_max
        
        # Kinetic Rates (The Rescued Open System)
        self.k_form = 0.01
        self.k_on1, self.k_off1 = 0.05, 0.0001   # Stochastic k_on, SpCas9 baseline k_off
        self.k_on2, self.k_off2 = 0.005, 0.001   # Weak decoy affinity
        self.k_tx_dcas9, self.k_tx_sgrna = 0.05, 0.1 
        
        # --- BIOLOGICALLY CORRECT DEGRADATION SPLIT ---
        self.k_dilution = 0.000385      # 30-minute cell doubling time (dCas9 & Complex)
        self.grna_deg_rate = 0.005      # Matches app.py slider! Fast RNase degradation.
        
        # Target constraints
        self.plasmid_copy_number = 2    # Matches app.py slider! Primary target baseline.
        self.T2_init = decoy_load
        
        self.reset_system()

    def reset_system(self):
        """Resets the biological environment to absolute zero."""
        self.t = 0.0
        self.dcas9, self.sgrna, self.C = 0, 0, 0
        
        # Dynamically pull from the UI sliders
        self.T1, self.R1 = self.plasmid_copy_number, 0
        self.T2, self.R2 = self.T2_init, 0
        
        self.time_points = [self.t]
        self.T1_points = [self.R1]
        self.T2_points = [self.R2]
        self.C_points = [self.C]

    def run(self, shock_time=None, decoy_spike=0):
        shock_applied = False
        
        while self.t < self.t_max:
            # Inject the competitive retroactivity stress if the timer hits the threshold
            if shock_time is not None and self.t >= shock_time and not shock_applied:
                self.T2 += decoy_spike
                shock_applied = True
                
            # Pure Python scalar propensities
            a0_prop = self.k_form * self.dcas9 * self.sgrna
            a1_prop = self.k_on1 * self.C * self.T1
            a2_prop = self.k_off1 * self.R1
            a3_prop = self.k_on2 * self.C * self.T2
            a4_prop = self.k_off2 * self.R2
            a5_prop = self.k_tx_dcas9
            a6_prop = self.k_tx_sgrna
            
            # --- APPLYING DECOUPLED DEGRADATION ---
            a7_prop = self.k_dilution * self.dcas9       # dCas9 dilutes via cell division
            a8_prop = self.grna_deg_rate * self.sgrna    # sgRNA actively degrades (UI Controlled)
            a9_prop = self.k_dilution * self.C           # Complex dilutes via cell division
            
            props = (a0_prop, a1_prop, a2_prop, a3_prop, a4_prop, a5_prop, a6_prop, a7_prop, a8_prop, a9_prop)
            a_sum = sum(props)
            
            if a_sum <= 0: break
            
            # Time step calculation
            r1, r2 = np.random.random(2)
            self.t += -np.log(r1) / a_sum
            
            # Fast pure Python selection to avoid np.cumsum overhead
            target = r2 * a_sum
            current = 0.0
            reaction_index = 0
            for i, p in enumerate(props):
                current += p
                if current > target:
                    reaction_index = i
                    break
            
            # Execute Winning Reaction
            if reaction_index == 0: self.dcas9 -= 1; self.sgrna -= 1; self.C += 1
            elif reaction_index == 1: self.C -= 1; self.T1 -= 1; self.R1 += 1
            elif reaction_index == 2: self.C += 1; self.T1 += 1; self.R1 -= 1
            elif reaction_index == 3: self.C -= 1; self.T2 -= 1; self.R2 += 1
            elif reaction_index == 4: self.C += 1; self.T2 += 1; self.R2 -= 1
            elif reaction_index == 5: self.dcas9 += 1
            elif reaction_index == 6: self.sgrna += 1
            elif reaction_index == 7: self.dcas9 -= 1
            elif reaction_index == 8: self.sgrna -= 1
            elif reaction_index == 9: self.C -= 1
            
            self.time_points.append(self.t)
            self.T1_points.append(self.R1)
            self.T2_points.append(self.R2)
            self.C_points.append(self.C)

# --- LOCAL TESTING BLOCK ---
if __name__ == "__main__":
    sim = CRISPRiEngine(decoy_load=25, t_max=3600)
    plt.figure(figsize=(9, 5.5))
    
    print("Executing Single Cell Trajectory for local verification...")
    sim.run(shock_time=1800, decoy_spike=50) 
    
    plt.step(sim.time_points, sim.T1_points, where='post', color='#00ff00', label='Primary Target')
    plt.step(sim.time_points, sim.T2_points, where='post', color='#ff0000', label='Decoy Load')
    plt.xlabel("Time (seconds)")
    plt.ylabel("Discrete Molecular Counts")
    plt.title("CRISPRi Circuit Dynamics (Decoupled RNA/Protein Degradation)")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()