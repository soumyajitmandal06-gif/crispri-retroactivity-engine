import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from crispri_engine import CRISPRiEngine

# 1. UI Configuration & Specifications
st.set_page_config(page_title="RetroSim", layout="wide")
st.title("RetroSim")

st.markdown("""
**Engine Specifications:**
* **Mathematical Framework:** 10-reaction Gillespie Stochastic Simulation Algorithm (SSA).
* **Biological Scope:** Quantifies resource competition, kinetic rescue, and target flatlining under **competitive retroactivity (sequestration)**.
* **Architecture:** Compares standard nucleases against orthogonal and AI-designed variants.
---
""")

# 1.5 System Architecture (Sidebar)
st.sidebar.header("0. System Architecture")
cas_variant = st.sidebar.selectbox(
    "Select Cas Enzyme Profile",
    ("SpCas9 (Baseline)", "Nme1Cas9 (High-Affinity/AI-Designed)", "Miniature Cas12f")
)
plasmid_copy_number = st.sidebar.slider("Plasmid Copy Number", min_value=1, max_value=100, value=25)
grna_deg_rate = st.sidebar.slider("gRNA Degradation Rate (1/s)", min_value=0.001, max_value=0.100, value=0.010, format="%.3f")

# 2. User Input Panel (Sidebar)
st.sidebar.header("1. Sequestration Load")
user_decoy_load = st.sidebar.slider("Decoy Target Load (Copies)", min_value=0, max_value=100, value=25)
user_t_max = st.sidebar.number_input("Simulation Time (seconds)", min_value=600, max_value=7200, value=3600, step=600)

st.sidebar.header("2. Engineering Strategy")
strategy = st.sidebar.radio(
    "Circuit Kinetics", 
    ["Wild-Type (Baseline Affinity)", "Kinetic Rescue (10x Affinity)"]
)
st.sidebar.markdown("---")
st.sidebar.header("2.5 Dynamic Perturbation")
shock_time = st.sidebar.slider("Time to Inject Shock (s)", 0, int(user_t_max), int(user_t_max/2))
decoy_spike = st.sidebar.slider("Decoy Target Spike (Copies)", 0, 100, 50)
st.sidebar.header("3. Execution")
sim_mode = st.sidebar.radio("Simulation Mode", ["Single Cell Trajectory", "Monte Carlo Ensemble (50 cells)"])

# 3. Execution Engine
if st.sidebar.button("Run Simulation", type="primary"):
    # Initialize the core engine with the base parameters
    sim = CRISPRiEngine(decoy_load=user_decoy_load, t_max=user_t_max)
    
    # Wire the Cas variant dropdown to the kinetic parameters
    if cas_variant == "SpCas9 (Baseline)":
        sim.k_on1 = 0.05
        sim.k_off1 = 0.0001      # Matches the methodology PDF
    elif cas_variant == "Nme1Cas9 (High-Affinity/AI-Designed)":
        sim.k_on1 = 0.1          
        sim.k_off1 = 0.00001     # Pushed below k_deg boundary
    elif cas_variant == "Miniature Cas12f":
        sim.k_on1 = 0.01 
        sim.k_off1 = 0.005
        
    # Wire the UI sliders directly into the engine
    sim.plasmid_copy_number = plasmid_copy_number
    sim.grna_deg_rate = grna_deg_rate
    
    # CRITICAL: Force a reset so the engine absorbs the new UI values before running
    sim.reset_system()
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    if sim_mode == "Single Cell Trajectory":
        with st.spinner("Calculating stochastic matrix..."):
            sim.run(shock_time=shock_time, decoy_spike=decoy_spike)
            ax.step(sim.time_points, sim.T1_points, where='post', color='#00ff00', linewidth=2.5, label='Primary Target')
            ax.step(sim.time_points, sim.T2_points, where='post', color='#ff0000', linewidth=2, label=f'Decoy Load ({user_decoy_load})')
    
    elif sim_mode == "Monte Carlo Ensemble (50 cells)":
        with st.spinner(f"Executing 50-cell Monte Carlo ensemble under {strategy}..."):
            for i in range(50):
                sim.reset_system()
                sim.run(shock_time=shock_time, decoy_spike=decoy_spike)
                if i == 0:
                    ax.step(sim.time_points, sim.T1_points, where='post', color='#00ff00', alpha=0.15, label='Primary Target')
                    ax.step(sim.time_points, sim.T2_points, where='post', color='#ff0000', alpha=0.15, label=f'Decoy Load ({user_decoy_load})')
                else:
                    ax.step(sim.time_points, sim.T1_points, where='post', color='#00ff00', alpha=0.15)
                    ax.step(sim.time_points, sim.T2_points, where='post', color='#ff0000', alpha=0.15)

    # Format Output Graph
    ax.set_xlabel("Time (seconds)")
    ax.set_ylabel("Discrete Molecular Counts")
    ax.set_title(f"CRISPRi Circuit Dynamics: {strategy}")
    ax.legend(loc="upper left")
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-1, max(user_decoy_load + decoy_spike + 5, 10))
    
    # 4. Render Output
    st.pyplot(fig)
    
    # --- Quantitative Metric Callouts ---
    st.divider()
    st.subheader("Simulation Metrics")
    
    col1, col2, col3 = st.columns(3)
    
    final_target_level = sim.T1_points[-1] 
    
    with col1:
        st.metric(label="Final Target Repression", value=f"{final_target_level:.1f} copies")
    with col2:
        st.metric(label="Sequestration Load", value=f"{user_decoy_load} decoys")
    with col3:
        st.metric(label="Cas Variant Profile", value=cas_variant.split()[0])

    # --- Data Export (.csv) ---
    st.write("")
    
    df_results = pd.DataFrame({
        "Time (s)": sim.time_points,
        "Primary Target Count": sim.T1_points,
        "Decoy Count": sim.T2_points
    })
    
    st.download_button(
        label="Download Trajectory Data (.csv)",
        data=df_results.to_csv(index=False).encode('utf-8'),
        file_name="retrosim_data.csv",
        mime="text/csv",
        type="primary"
    )