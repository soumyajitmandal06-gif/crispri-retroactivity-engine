import psycopg2
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import sys

# ---> PASTE YOUR EXACT NEON URL HERE <---
NEON_URL = "postgresql://neondb_owner:npg_gDUX5atMH2bz@ep-cool-butterfly-adfbjv0g-pooler.c-2.us-east-1.aws.neon.tech/neondb?sslmode=require&channel_binding=require"

def fetch_parameters():
    try:
        conn = psycopg2.connect(NEON_URL)
        cur = conn.cursor()
        cur.execute("SELECT parameter_type, min_value, max_value FROM kinetic_constants WHERE part_id = 'genome_ecoli_mg1655'")
        rows = cur.fetchall()
        params = {row[0]: (float(row[1]), float(row[2])) for row in rows}
        cur.close()
        conn.close()
        return params
    except Exception as e:
        print(f"--> [DATABASE ERROR] {e}")
        sys.exit(1)

def retroactivity_dynamics(t, y, k_tx1, k_tx2, k_tl_per_ribo, d_m, d_p, R_total, Kd):
    mRNA1, Protein1, mRNA2, Protein2 = y
    
    dmRNA1_dt = k_tx1 - d_m * mRNA1
    dmRNA2_dt = k_tx2 - d_m * mRNA2
    
    allocation_1 = (mRNA1 / Kd) / (1.0 + (mRNA1 / Kd) + (mRNA2 / Kd))
    allocation_2 = (mRNA2 / Kd) / (1.0 + (mRNA1 / Kd) + (mRNA2 / Kd))
    
    dProtein1_dt = (k_tl_per_ribo * R_total * allocation_1) - (d_p * Protein1)
    dProtein2_dt = (k_tl_per_ribo * R_total * allocation_2) - (d_p * Protein2)
    
    return [dmRNA1_dt, dProtein1_dt, dmRNA2_dt, dProtein2_dt]

if __name__ == "__main__":
    print("--> Fetching biological parameters from Neon Cloud...")
    bio_params = fetch_parameters()

    half_life_min = np.mean(bio_params['mRNA Half-Life'])
    ribosome_rate = np.mean(bio_params['Ribosome Translation Rate'])
    R_total = np.mean(bio_params['Total Free Ribosomes'])

    protein_length = 300 
    k_tl_per_ribo = ribosome_rate / protein_length 
    d_m = np.log(2) / (half_life_min * 60)
    d_p = 0.0003 
    Kd = 50.0 

    k_tx1 = 0.05  
    k_tx2 = 0.50  

    print("--> Initializing SciPy Retroactivity Solver (Radau Method)...")
    t_span = (0, 7200) 
    initial_conditions = [0, 0, 0, 0] 

    solution = solve_ivp(
        retroactivity_dynamics, 
        t_span, 
        initial_conditions, 
        args=(k_tx1, k_tx2, k_tl_per_ribo, d_m, d_p, R_total, Kd), 
        dense_output=True,
        method='Radau' 
    )

    print("-" * 50)
    print(f"RETROACTIVITY METRICS (At T=2 Hours):")
    print(f"  Primary Circuit Protein:   {solution.y[1][-1]:,.2f} molecules")
    print(f"  Parasitic Circuit Protein: {solution.y[3][-1]:,.2f} molecules")
    print("-" * 50)

    # --- NEW MATPLOTLIB VISUALIZATION ---
    print("--> Generating Visualization...")
    t_minutes = solution.t / 60 

    plt.figure(figsize=(10, 6))
    plt.plot(t_minutes, solution.y[1], label="Primary Circuit Protein", color="blue", linewidth=2.5)
    plt.plot(t_minutes, solution.y[3], label="Parasitic Circuit Protein", color="red", linestyle="--", linewidth=2.5)
    
    plt.title("Resource Competition & Retroactivity in E. coli", fontsize=14, fontweight="bold")
    plt.xlabel("Time (Minutes)", fontsize=12)
    plt.ylabel("Protein Molecule Count", fontsize=12)
    plt.legend(loc="upper left")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    plt.savefig("retroactivity_graph.png", dpi=300)
    print("--> [SUCCESS] Graph saved as 'retroactivity_graph.png'")
    plt.show()
    import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# 1. Kinetic Parameters
k_form = 0.002  # Complex formation rate
k_on1, k_off1 = 0.05, 0.001  # Target 1 binding/unbinding
k_on2, k_off2 = 0.05, 0.001  # Target 2 binding/unbinding (identical affinity)

# 2. The Differential Equation System
def retroactivity_system(t, y):
    dcas9, sgrna, C, T1, R1, T2, R2 = y
    
    # Reaction velocities
    formation = k_form * dcas9 * sgrna
    bind1 = k_on1 * C * T1 - k_off1 * R1
    bind2 = k_on2 * C * T2 - k_off2 * R2
    
    # ODEs
    d_dcas9 = -formation
    d_sgrna = -formation
    d_C = formation - bind1 - bind2
    d_T1 = -bind1
    d_R1 = bind1
    d_T2 = -bind2
    d_R2 = bind2
    
    return [d_dcas9, d_sgrna, d_C, d_T1, d_R1, d_T2, d_R2]

# 3. Initial States: 10 dCas9, 15 sgRNA, 2 Target 1, 15 Target 2 (The Load)
y0 = [10, 15, 0, 2, 0, 15, 0] 

# 4. Execute the Radau Solver
sol = solve_ivp(retroactivity_system, [0, 300], y0, method='Radau', dense_output=True)
t = np.linspace(0, 300, 500)
y = sol.sol(t)

# 5. Visualize the Circuit Failure
plt.figure(figsize=(8, 5))
plt.plot(t, y[4], label='Repressed Target 1 (Intended)', color='#00ff00', linewidth=2.5)
plt.plot(t, y[6], label='Repressed Target 2 (Decoy Load)', color='#ff0000', linewidth=2.5)
plt.plot(t, y[2], label='Free Complexes', color='blue', linestyle='--', linewidth=1.5)

plt.xlabel('Time (seconds)')
plt.ylabel('Molecular Counts')
plt.title('CRISPRi Retroactivity: Circuit Failure Under Decoy Load')
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()