# In Silico Analysis of CRISPRi Circuit Robustness 

**Author:** Soumyajit Mandal
**Domain:** Computational Synthetic Biology / Systems Biology

## Overview
This repository contains a fully stochastic biophysical engine designed to model resource competition (retroactivity) in CRISPR interference (CRISPRi) genetic circuits. Built entirely in Python, the engine transitions from deterministic mass-action kinetics to a discrete multi-reaction Gillespie algorithm to accurately capture intrinsic molecular noise at low intracellular copy numbers.

The primary objective of this computational pipeline is to map the exact failure point of a primary genetic circuit under a heavy decoy load, and mathematically prove that the circuit can be rescued via kinetic affinity tuning without increasing the metabolic burden on the host cell.

## Biophysical Architecture

The engine simulates a dynamic, open-system *E. coli* environment computing 10 competing reactions simultaneously:
1. **Transcription & Translation:** Zero-order continuous production of dCas9 and sgRNA.
2. **Cellular Dilution:** First-order degradation and dilution representing cell growth.
3. **Complex Formation:** Stochastic binding of dCas9 and sgRNA.
4. **Primary Target Binding:** Reversible binding kinetics ($k_{on1}$, $k_{off1}$) to the intended circuit.
5. **Decoy Target Binding:** Reversible binding kinetics ($k_{on2}$, $k_{off2}$) to secondary retroactive loads.

## The Retroactivity Problem & Kinetic Rescue
In synthetic biology, retroactivity acts as a massive metabolic sink. When a dCas9-sgRNA pool is strictly constrained (e.g., 10 molecules), introducing high-copy decoy binding sites mathematically guarantees primary circuit failure. The decoys act as a molecular sponge.

This engine demonstrates **Kinetic Rescue**: by computationally mutating the primary target to possess a 10x stronger thermodynamic binding affinity ($k_{on1} = 0.5$, $k_{off1} = 0.0001$), the primary circuit successfully hijacks the limited repressor pool back from the decoys, restoring full repression without requiring the synthesis of additional dCas9 proteins.

## Object-Oriented Implementation

The core biophysics are encapsulated in a robust Python class, allowing for dynamic parameter sweeps and stress-testing across variable cellular environments.

### Usage
```python
from crispri_engine import CRISPRiEngine
import matplotlib.pyplot as plt

# Initialize the cell environment with a severe 25-copy decoy load
sim = CRISPRiEngine(decoy_load=25, t_max=400)

# Execute the 10-reaction Gillespie matrix
sim.run()

# The system automatically tracks time_points, T1_points, and T2_points for plotting