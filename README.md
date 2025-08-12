# Rectangular Microstrip Patch Antenna Calculator

## Author
**Sanju S Pillai** (`Sherlock-191b`)

---

## 📌 Purpose
This project calculates the **physical** and **electrical** parameters of a rectangular microstrip patch antenna based on user-provided design specifications.  
It is useful for **RF/microwave engineers** and **researchers** working in antenna design.

---

## 📥 Inputs
- Operating frequency (in MHz or GHz)
- Dielectric constant (εr) of the substrate
- Substrate thickness (in mm or cm)
- Desired characteristic impedance (R₀)
- Selection of desired unit for final output (m, cm, mm)

---

## ⚙️ Features
- GUI interface for intuitive user interaction
- Calculates:
  - Patch width (W)  
  - Effective dielectric constant (ε_eff)  
  - Effective and actual patch length (L_eff, L)  
  - Fringing length extension (ΔL)  
  - Ground plane dimensions (Lg, Wg)  
  - Radiation conductances (G₁, G₁₂)  
  - Input resistance (Rin) and optimal feed location (y₀)
- Unit customization for final results output (m/cm/mm)
- Clean, text-based results display
- Educational introduction screen with design equations

---

## 🖥️ Requirements
- Python 3.x
- `tkinter`
- `numpy`
- `scipy`

Install dependencies:
```bash
pip install numpy scipy
