"""
====================================================================================
        RECTANGULAR MICROSTRIP PATCH ANTENNA DESIGN CALCULATOR
====================================================================================

    Author  : SANJU S PILLAI
    Purpose : To calculate the physical and electrical parameters of a rectangular 
              microstrip patch antenna based on input design specifications.
              Useful for RF/microwave engineers and researchers working in antenna design.

    Inputs  :
        - Operating frequency (in MHz or GHz)
        - Dielectric constant (εr) of the substrate
        - Substrate thickness (in mm or cm)
        - Desired characteristic impedance (R₀)
        - Selection of desired unit for final output (m, cm, mm)

    Features:
        - GUI interface for intuitive user interaction
        - Calculates patch width (W), effective dielectric constant (ε_eff)
        - Calculates effective and actual patch length (L_eff, L)
        - Computes fringing length extension (ΔL)
        - Calculates ground plane dimensions (Lg, Wg)
        - Computes radiation conductance terms (G₁, G₁₂)
        - Computes input resistance (Rin) and optimal feed location (y₀)
        - Unit customization for final results output (m/cm/mm)
        - Option to review all results in a clean text-based display
        - Educational introduction screen with design equations

    Version : 1.0
    Date    : June 2025
    License : Open and free to use, modify, and improve — for learning and research.

====================================================================================
"""


import tkinter as tk
from tkinter import ttk
import numpy as np
from scipy.constants import c, mu_0
from scipy.integrate import quad

# Helper Functions
def calc_patch_antenna_params(frequency, freq_unit, epsilon_r, h, h_unit, R_0):
    freq = frequency * 1e9 if freq_unit == "GHz" else frequency * 1e6
    h = h / 100.0 if h_unit == "cm" else h / 1000.0

    W = (c / (2 * freq)) * np.sqrt(2 / (epsilon_r + 1))
    epsilon_eff = (epsilon_r + 1)/2 + (epsilon_r - 1)/2 * (1 / np.sqrt(1 + 12 * h / W))
    delta_L = h * 0.412 * ((epsilon_eff + 0.3)*(W/h + 0.264)) / ((epsilon_eff - 0.258)*(W/h + 0.8))
    L_eff = c / (2 * freq * np.sqrt(epsilon_eff))
    L = L_eff - 2 * delta_L
    Lg = L + 6 * h
    Wg = W + 6 * h

    k0 = 2 * np.pi * freq / c

    def integrand_G1(theta):
        numerator = np.sin(k0 * W * np.cos(theta) / 2)**2
        denominator = np.cos(theta)**2
        return (numerator / denominator) * np.sin(theta)

    def integrand_G12(theta):
        term1 = np.sin(k0 * W * np.cos(theta) / 2)**2 / (np.cos(theta)**2)
        term2 = np.cos(k0 * L * np.sin(theta))
        return term1 * term2 * np.sin(theta)

    epsilon = 1e-4
    G1_1, _ = quad(integrand_G1, 0, np.pi/2 - epsilon)
    G1_2, _ = quad(integrand_G1, np.pi/2 + epsilon, np.pi)
    G1 = (1 / (120 * np.pi**2)) * (G1_1 + G1_2)

    G12_1, _ = quad(integrand_G12, 0, np.pi/2 - epsilon)
    G12_2, _ = quad(integrand_G12, np.pi/2 + epsilon, np.pi)
    G12 = (1 / (120 * np.pi**2)) * (G12_1 + G12_2)

    Rin = 1 / (2 * G1 + 2 * G12)
    ratio = R_0 / Rin
    y0 = (L / np.pi) * np.arccos(np.sqrt(ratio)) if 0 <= ratio <= 1 else float('nan')

    return {
        "Patch Width (W)": W, "Effective Dielectric Constant (ε_eff)": epsilon_eff, "Effective Length (L_eff)": L_eff, 
        "Length Extension (ΔL)": delta_L, "Patch Length (L)": L, "Ground Length (Lg)": Lg, "Ground Width (Wg)": Wg,
        "G1": G1, "G12": G12, "Input Resistance (Rin)": Rin, "Feed Location (y₀)": y0
    }

# GUI Code
class PatchAntennaGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Microstrip Patch Antenna Calculator")
        self.show_intro()

    def show_intro(self):
        self.clear_frame()
        tk.Label(self.root, text="Microstrip Patch Antenna Designer", font=("Arial", 16, "bold")).pack(pady=10)

        description = tk.Text(self.root, wrap=tk.WORD, height=25, font=("Arial", 10))
        description.insert(tk.END, "This tool calculates parameters of a rectangular microstrip patch antenna based on input frequency, dielectric constant, substrate height, and impedance.\n\n")
        description.insert(tk.END, "Inputs:\n- Frequency (GHz or MHz)\n- Dielectric Constant (εr)\n- Substrate Height (mm or cm)\n- Characteristic Impedance (R₀)\n\n")
        description.insert(tk.END, "Outputs:\n- Patch Width and Length\n- Effective Dielectric Constant\n- Ground Plane Dimensions\n- Integral values (G1 and G12)\n- Input Resistance and Feed Location\n\n")
        description.insert(tk.END, "Equations:\n")
        description.insert(tk.END, "W = (c / 2f) * sqrt(2 / (εr + 1))\n")
        description.insert(tk.END, "ε_eff = (εr + 1)/2 + (εr - 1)/2 * 1 / sqrt(1 + 12h/W)\n")
        description.insert(tk.END, "ΔL = h * 0.412 * [(ε_eff + 0.3)(W/h + 0.264)] / [(ε_eff - 0.258)(W/h + 0.8)]\n")
        description.insert(tk.END, "L_eff = c / (2f * sqrt(ε_eff)), L = L_eff - 2ΔL\n")
        description.insert(tk.END, "y₀ = (L / π) * arccos(sqrt(R₀ / Rin))\n")
        description.config(state=tk.DISABLED)
        description.pack(padx=20)

        tk.Button(self.root, text="Continue", command=self.show_input_form).pack(pady=10)
        self.footer()

    def show_input_form(self):
        self.clear_frame()
        tk.Label(self.root, text="Enter Design Parameters:", font=("Arial", 12, "bold")).pack(pady=10)

        self.inputs = {
            'frequency': tk.DoubleVar(), 'freq_unit': tk.StringVar(value="GHz"),
            'epsilon_r': tk.DoubleVar(), 'h': tk.DoubleVar(), 'h_unit': tk.StringVar(value="mm"),
            'R_0': tk.DoubleVar(value=50.0)
        }
        self.output_unit = tk.StringVar(value="m")

        form = tk.Frame(self.root)
        form.pack()
        self.create_input(form, "Frequency:", self.inputs['frequency'], self.inputs['freq_unit'], ["GHz", "MHz"])
        self.create_input(form, "Dielectric Constant (εr):", self.inputs['epsilon_r'])
        self.create_input(form, "Substrate Height:", self.inputs['h'], self.inputs['h_unit'], ["mm", "cm"])
        self.create_input(form, "Characteristic Impedance (R₀):", self.inputs['R_0'])

        tk.Button(self.root, text="Calculate", command=self.calculate).pack(pady=10)
        self.footer()

    def create_input(self, frame, label, variable, unit_var=None, unit_options=None):
        row = tk.Frame(frame)
        row.pack(pady=5)
        tk.Label(row, text=label, width=30, anchor='w').pack(side=tk.LEFT)
        tk.Entry(row, textvariable=variable, width=15).pack(side=tk.LEFT)
        if unit_options:
            ttk.Combobox(row, textvariable=unit_var, values=unit_options, width=6).pack(side=tk.LEFT)

    def calculate(self):
        self.results_raw = calc_patch_antenna_params(
            self.inputs['frequency'].get(), self.inputs['freq_unit'].get(),
            self.inputs['epsilon_r'].get(), self.inputs['h'].get(),
            self.inputs['h_unit'].get(), self.inputs['R_0'].get()
        )
        self.show_results()

    def show_results(self):
        self.clear_frame()
        tk.Label(self.root, text="Calculated Parameters:", font=("Arial", 12, "bold")).pack(pady=10)

        result_frame = tk.Frame(self.root)
        result_frame.pack()

        tk.Label(result_frame, text="Select Output Unit:").pack()
        unit_selector = ttk.Combobox(result_frame, values=["m", "cm", "mm"], textvariable=self.output_unit)
        unit_selector.pack()
        unit_selector.bind("<<ComboboxSelected>>", lambda e: self.update_results())

        self.result_box = tk.Text(result_frame, height=25, width=70, font=("Consolas", 10))
        self.result_box.pack()
        self.update_results()
        self.footer()

    def update_results(self):
        unit = self.output_unit.get()
        factor = 1 if unit == "m" else (100 if unit == "cm" else 1000)

        self.result_box.config(state=tk.NORMAL)
        self.result_box.delete("1.0", tk.END)

        def fmt(val): return f"{val:.6f}"

        for key, val in self.results_raw.items():
            if any(term in key for term in ["Width", "Length", "y₀", "ΔL"]):
                self.result_box.insert(tk.END, f"{key} [{unit}]: {fmt(val * factor)}\n")
            else:
                self.result_box.insert(tk.END, f"{key}: {fmt(val)}\n")

        self.result_box.config(state=tk.DISABLED)

    def clear_frame(self):
        for widget in self.root.winfo_children():
            widget.destroy()

    def footer(self):
        tk.Label(self.root, text="Created by sherlock191b", font=("Arial", 8), anchor="se").pack(side=tk.BOTTOM, anchor="se")

if __name__ == '__main__':
    root = tk.Tk()
    gui = PatchAntennaGUI(root)
    root.mainloop()