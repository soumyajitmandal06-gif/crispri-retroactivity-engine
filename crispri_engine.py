import numpy as np
import matplotlib.pyplot as plt

class CRISPRiEngine:
    def __init__(self, decoy_load=15, t_max=600.0):
        """Initializes the engine with user-defined cellular stress constraints."""
        self.t_max = t_max
        
        # Kinetic Rates (The Rescued Open System)
        self.k_form = 0.002
        self.k_on1, self.k_off1 = 0.5, 0.0001   # Strong primary affinity
        self.k_on2, self.k_off2 = 0.05, 0.001   # Weak decoy affinity
        self.k_tx_dcas9, self.k_tx_sgrna = 0.5, 0.7
        self.k_deg = 0.01
        
        # Target constraints
        self.T1_init = 2
        self.T2_init = decoy_load
        
        self.reset_system()

    def reset_system(self):
        """Resets the biological environment to absolute zero."""
        self.t = 0.0
        self.dcas9, self.sgrna, self.C = 0, 0, 0
        self.T1, self.R1 = self.T1_init, 0
        self.T2, self.R2 = self.T2_init, 0
        
        self.time_points = [self.t]
        self.T1_points = [self.R1]
        self.T2_points = [self.R2]
        self.C_points = [self.C]

    def run(self):
        """Executes the 10-reaction stochastic Gillespie matrix."""
        while self.t < self.t_max:
            # Propensity Calculations
            props = np.array([
                self.k_form * self.dcas9 * self.sgrna,  # 0: Form complex
                self.k_on1 * self.C * self.T1,          # 1: Bind T1
                self.k_off1 * self.R1,                  # 2: Unbind T1
                self.k_on2 * self.C * self.T2,          # 3: Bind T2
                self.k_off2 * self.R2,                  # 4: Unbind T2
                self.k_tx_dcas9,                        # 5: Produce dCas9
                self.k_tx_sgrna,                        # 6: Produce sgRNA
                self.k_deg * self.dcas9,                # 7: Degrade dCas9
                self.k_deg * self.sgrna,                # 8: Degrade sgRNA
                self.k_deg * self.C                     # 9: Degrade Complex
            ])
            
            a_0 = np.sum(props)
            if a_0 <= 0: break
            
            r1, r2 = np.random.random(2)
            self.t += -np.log(r1) / a_0
            reaction_index = np.searchsorted(np.cumsum(props / a_0), r2)
            
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

# --- USER INTERACTION BLOCK ---
# --- USER INTERACTION BLOCK ---
if __name__ == "__main__":
    # Initialize the engine
    sim = CRISPRiEngine(decoy_load=25, t_max=400)
    
    num_cells = 50
    plt.figure(figsize=(9, 5.5))
    
    print(f"Executing Monte Carlo Ensemble: {num_cells} cells...")
    
    # Run 50 independent stochastic simulations (The Monte Carlo Ensemble)
    for i in range(num_cells):
        sim.reset_system()
        sim.run()
        
        # Plot each cell's trajectory with low opacity to build the noise cloud
        if i == 0:
            # Add labels only for the first loop so the legend stays clean
            plt.step(sim.time_points, sim.T1_points, where='post', color='#00ff00', alpha=0.15, label='Primary Target (Rescued)')
            plt.step(sim.time_points, sim.T2_points, where='post', color='#ff0000', alpha=0.15, label='Decoy Load (25 Copies)')
        else:
            plt.step(sim.time_points, sim.T1_points, where='post', color='#00ff00', alpha=0.15)
            plt.step(sim.time_points, sim.T2_points, where='post', color='#ff0000', alpha=0.15)

    plt.xlabel('Time (seconds)')
    plt.ylabel('Discrete Molecular Counts')
    plt.title(f'Monte Carlo Stochastic Ensemble ({num_cells} Simulated E. coli Cells)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.show()
    print("Ensemble complete.")
