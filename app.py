
import streamlit as st
import matplotlib.pyplot as plt
from crispri_engine import CRISPRiEngine

# 1. UI Configuration
st.set_page_config(page_title="CRISPRi Retroactivity Engine", layout="wide")
st.title("In Silico CRISPRi Retroactivity Simulator")
st.markdown("Test the robustness of a CRISPRi genetic circuit under metabolic stress, and apply thermodynamic affinity tuning to rescue the system.")

# 2. User Input Panel (Sidebar)
st.sidebar.header("1. Environmental Stress")
user_decoy_load = st.sidebar.slider("Decoy Target Load (Copies)", min_value=0, max_value=100, value=25)
user_t_max = st.sidebar.number_input("Simulation Time (seconds)", min_value=100, max_value=2000, value=400, step=100)

st.sidebar.header("2. Engineering Strategy")
strategy = st.sidebar.radio(
    "Circuit Kinetics", 
    ["Wild-Type (Baseline Affinity)", "Kinetic Rescue (10x Affinity)"]
)

st.sidebar.header("3. Execution")
sim_mode = st.sidebar.radio("Simulation Mode", ["Single Cell Trajectory", "Monte Carlo Ensemble (50 cells)"])

# 3. Execution Engine
if st.sidebar.button("Run Simulation", type="primary"):
    # Initialize the core engine
    sim = CRISPRiEngine(decoy_load=user_decoy_load, t_max=user_t_max)
    
    # Inject the user's chosen kinetic strategy into the engine
    if strategy == "Wild-Type (Baseline Affinity)":
        sim.k_on1 = 0.05
        sim.k_off1 = 0.001
    else:
        sim.k_on1 = 0.5      # 10x stronger binding
        sim.k_off1 = 0.0001  # 10x slower unbinding
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    if sim_mode == "Single Cell Trajectory":
        with st.spinner("Calculating stochastic matrix..."):
            sim.run()
            ax.step(sim.time_points, sim.T1_points, where='post', color='#00ff00', linewidth=2.5, label='Primary Target')
            ax.step(sim.time_points, sim.T2_points, where='post', color='#ff0000', linewidth=2, label=f'Decoy Load ({user_decoy_load})')
    
    elif sim_mode == "Monte Carlo Ensemble (50 cells)":
        with st.spinner(f"Executing 50-cell Monte Carlo ensemble under {strategy}..."):
            for i in range(50):
                sim.reset_system()
                sim.run()
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
    ax.set_ylim(-1, max(user_decoy_load + 5, 10))
    
    # 4. Render Output
    st.pyplot(fig)