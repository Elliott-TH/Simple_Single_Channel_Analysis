# Single Channel Analysis (SCA)

**Author:** Elliott Hirko  
**Date:** February 2026

---

## Introduction

This code performs a simple steady-state single channel thermal-hydraulic analysis of a fuel pin in a Light Water Reactor (LWR).  
It includes:

- A water property table  
- An input Excel file for user-specified operating and geometric parameters  

---

## Simple Use

A user can specify the following inputs in the Excel file:

### Input Parameters

| Symbol | Description | Units |
|--------|-------------|-------|
| Type | Light Water Reactor type (PWR or BWR) | -- |
| $L$ | Channel length | m |
| $L_e$ | Extrapolated channel length | m |
| $D$ | Fuel pin outer diameter | m |
| Pitch | Lattice pitch (center-to-center pin spacing) | m |
| $T_{m,\text{in}}$ | Coolant inlet temperature | K |
| $\dot{m}$ | Coolant mass flow rate | kg/s |
| $P_{\text{nom}}$ | Nominal system pressure | Pa |
| $q_0$ | Peak linear heat generation rate | W/m |
| $D_{ci}$ | Cladding inner diameter | m |
| $D_{co}$ | Cladding outer diameter | m |
| SCB | Subchannel boiling model toggle (on/off) | -- |

---

## Output

Upon running the code, a CSV file will be generated containing the following outputs:

### Output Parameters

| Symbol | Description | Units |
|--------|-------------|-------|
| Cell | Axial computational cell index | -- |
| $z$ | Axial position along channel | m |
| $T_m$ | Coolant bulk (mixture) temperature | K |
| $T_{co}$ | Cladding outer surface temperature | K |
| $T_{ci}$ | Cladding inner surface temperature | K |
| $T_{fo}$ | Fuel outer surface temperature | K |
| $T_{\text{max}}$ | Maximum fuel centerline temperature | K |
| $x$ | Thermodynamic vapor quality | -- |
| $x_e$ | Equilibrium vapor quality | -- |
| CHFR | Critical Heat Flux Ratio | -- |
| $\Delta P$ | Pressure drop along channel | Pa |
