import psycopg2
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import sys

# ---> PASTE YOUR EXACT NEON URL HERE <---
NEON_URL = "postgresql://neondb_owner:npg_gDUX5atMH2bz@ep-cool-butterfly-adfbjv0g-pooler.c-2.us-east-1.aws.neon.tech/neondb?sslmode=require&channel_binding=require"

def fetch_all_parameters():
    try:
        conn = psycopg2.connect(NEON_URL)
        cur = conn.cursor()
        
        # Fetch E. coli baseline
        cur.execute("SELECT parameter_type, min_value, max_value FROM kinetic_constants WHERE part_id = 'genome_ecoli_mg1655'")
        ecoli_rows = cur.fetchall()
        params = {row[0]: np.mean([float(row[1]), float(row[2])]) for row in ecoli_rows}
        
        # Fetch CRISPRi parameters
        cur.execute("SELECT parameter_type, min_value, max_value FROM kinetic_constants WHERE part_id IN ('dcas9_protein', 'sgrna_transcript')")
        crispri_rows = cur.fetchall()
        params.update({row[0]: np.mean([float(row[1]), float(row[2])]) for row in crispri_rows})
        
        cur.close()
        conn.close()
        return params
    except Exception as e:
        print(f"--> [DATABASE ERROR] {e}")
        sys.exit(1)

def crispri_dynamics(t, y, k_tx_primary, d_m, d_p, k_on, k_off, k_tx_sgrna, k_tl_dcas9, Kd_repress):
    m_p, P_p, dCas9, sgRNA, Complex = y
    
    formation_rate = (k_on * dCas9 * sgRNA) - (k_off * Complex)
    
    d_dCas9_dt = k_tl_dcas9 - (d_p * dCas9) - formation_rate
    d_sgRNA_dt = k_tx_sgrna - (d_m * sgRNA) - formation_rate
    d_Complex_dt = formation_rate - (d_p * Complex)
    
    repression_factor = 1.0 / (1.0 + (Complex / Kd_repress)**2)
    
    d_m_p_dt = (k_tx_primary * repression_factor) - (d_m * m_p)
    d_P_p_dt = (10.0 * m_p) - (d_p * P_p) 
    
    return [d_m_p_dt, d_P_p_dt, d_dCas9_dt, d_sgRNA_dt, d_Complex_dt]

if __name__ == "__main__":
    print("--> Fetching all biological parameters from Neon Cloud...")
    p = fetch_all_parameters()

    d_m = np.log(2) / (p['mRNA Half-Life'] * 60)
    d_p = 0.0003 
    
    k_on = p['Complex Association (kon)']
    k_off = p['Complex Dissociation (koff)']
    k_tx_sgrna = p['Transcription Rate']
    k_tl_dcas9 = 0.5 
    Kd_repress = 100.0 
    k_tx_primary = 0.5 

    print("--> Initializing SciPy CRISPRi Solver (Radau Method)...")
    t_span = (0, 7200) 
    initial_conditions = [0, 0, 0, 0, 0] 

    solution = solve_ivp(
        crispri_dynamics, t_span, initial_conditions, 
        args=(k_tx_primary, d_m, d_p, k_on, k_off, k_tx_sgrna, k_tl_dcas9, Kd_repress), 
        dense_output=True, method='Radau'
    )

    print("-" * 50)
    print("CRISPRi STEADY-STATE METRICS (At T=2 Hours):")
    print(f"  Active Repressor Complex: {solution.y[4][-1]:,.2f} molecules")
    print(f"  Primary Target Protein:   {solution.y[1][-1]:,.2f} molecules")
    print("-" * 50)

    # --- MATPLOTLIB VISUALIZATION ---
    print("--> Generating CRISPRi Repression Visualization...")
    t_minutes = solution.t / 60 

    plt.figure(figsize=(10, 6))
    plt.plot(t_minutes, solution.y[1], label="Primary Target Protein", color="purple", linewidth=2.5)
    plt.plot(t_minutes, solution.y[4], label="Active dCas9-sgRNA Complex", color="teal", linestyle="--", linewidth=2.5)
    
    plt.title("CRISPR Interference & Target Repression Dynamics", fontsize=14, fontweight="bold")
    plt.xlabel("Time (Minutes)", fontsize=12)
    plt.ylabel("Molecule Count", fontsize=12)
    plt.legend(loc="upper right")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    plt.savefig("crispri_graph.png", dpi=300)
    print("--> [SUCCESS] Graph saved as 'crispri_graph.png'")
    plt.show()