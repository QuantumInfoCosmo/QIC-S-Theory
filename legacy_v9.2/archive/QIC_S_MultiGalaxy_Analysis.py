import numpy as np
import matplotlib.pyplot as plt
import os

"""
QIC-S Theory: Zero-Parameter Analysis (Definitive Edition)
Target: 8 Galaxies (Standard, LSB, Dwarf)
Author: Sasada (Independent Researcher)
Date: 2025-12-29 (Ver 4.1)
"""

# ==========================================
# 1. PHYSICS CONSTANTS & QIC-S PARAMETERS
# ==========================================
# Zero-Parameter Constraint (FIXED)
A0_THEORETICAL = 1.23e-10  # m/s^2 (Derived from cH0/2pi)
ML_DISK_FIXED = 0.5        # Standard stellar population
ML_GAS_FIXED = 1.0         # Gas mass factor

# Conversions
KPC_TO_M = 3.086e19
KMS_TO_MS = 1000.0

def qics_acceleration(g_bar_si):
    """
    QIC-S Formula: g_obs = g_bar / (1 - exp(-sqrt(g_bar/a0)))
    Returns g_tot in m/s^2
    """
    # Avoid zero division
    g_bar_safe = np.where(g_bar_si < 1e-15, 1e-15, g_bar_si)
    x = g_bar_safe / A0_THEORETICAL
    
    # Core QIC-S Logic (Interpolation function)
    g_tot_si = g_bar_safe / (1 - np.exp(-np.sqrt(x)))
    return g_tot_si

def load_galaxy_data(filename):
    """
    Loads .dat file (SPARC format assumed).
    Checks both current directory and 'data/' subdirectory.
    """
    # Check possible paths
    paths_to_check = [filename, os.path.join("data", filename)]
    target_path = None
    
    for p in paths_to_check:
        if os.path.exists(p):
            target_path = p
            break
            
    if target_path is None:
        print(f"Warning: {filename} not found in current folder or 'data/' folder.")
        return None
    
    try:
        data = np.loadtxt(target_path)
        # Assuming format: Rad(kpc), Vobs, Verr, Vgas, Vdisk, Vbul
        # Using np.abs to handle potential negative flags in raw data
        return {
            "r": data[:, 0],
            "v_obs": data[:, 1],
            "v_err": data[:, 2],
            "v_gas": np.abs(data[:, 3]), 
            "v_disk": np.abs(data[:, 4]),
            "v_bul": np.abs(data[:, 5]) if data.shape[1] > 5 else np.zeros_like(data[:, 0])
        }
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return None

def analyze_single_galaxy(galaxy_name, filename):
    data = load_galaxy_data(filename)
    if data is None: return None
    
    r = data["r"]
    
    # 1. Baryonic Acceleration (Zero Param: Fixed M/L)
    # V_bar^2 = V_disk^2 * ML + V_gas^2 * ML_gas + V_bul^2 * ML
    v_bar_sq = (data["v_disk"]**2 * ML_DISK_FIXED) + \
               (data["v_gas"]**2 * ML_GAS_FIXED) + \
               (data["v_bul"]**2 * ML_DISK_FIXED)
               
    g_bar = v_bar_sq / r # (km/s)^2 / kpc
    g_bar_si = g_bar * (KMS_TO_MS**2) / KPC_TO_M
    
    # 2. QIC-S Prediction
    g_tot_si = qics_acceleration(g_bar_si)
    g_tot = g_tot_si * KPC_TO_M / (KMS_TO_MS**2)
    v_qics = np.sqrt(g_tot * r)
    
    # 3. Statistics
    # Calculate deviation in the flat part (outer half of radius)
    mid_idx = len(r) // 2
    if mid_idx == 0: mid_idx = 0
    
    # Percentage difference
    # We use sum of differences to get a weighted feel or mean of percentages
    percentage_diff = (v_qics[mid_idx:] - data["v_obs"][mid_idx:]) / data["v_obs"][mid_idx:] * 100
    avg_deviation = np.mean(percentage_diff)
    
    rms = np.sqrt(np.mean((data["v_obs"] - v_qics)**2))
    
    return {
        "name": galaxy_name,
        "deviation": avg_deviation,
        "rms": rms,
        "type": "" # Placeholder
    }

# ==========================================
# MAIN EXECUTION
# ==========================================
def main():
    # Target List with Galaxy Types
    targets = [
        ("NGC 3198", "NGC3198_rotmod.dat", "Spiral"),
        ("NGC 7331", "NGC7331_rotmod.dat", "Massive Spiral"),
        ("NGC 2841", "NGC2841_rotmod.dat", "Massive Spiral"),
        ("NGC 6946", "NGC6946_rotmod.dat", "Spiral"),
        ("UGC 128",  "UGC00128_rotmod.dat", "LSB"),
        ("NGC 2903", "NGC2903_rotmod.dat", "Massive Spiral"),
        ("DDO 154",  "DDO154_rotmod.dat",   "Dwarf (gas-rich)"),
        ("IC 2574",  "IC2574_rotmod.dat",   "Dwarf (irregular)")
    ]
    
    results = []
    
    print("--- Starting QIC-S Zero-Parameter Analysis (Ver 4.1) ---")
    print(f"Constraints: a0={A0_THEORETICAL:.2e}, M/L_disk={ML_DISK_FIXED}")
    
    for name, fname, gtype in targets:
        print(f"Analyzing {name}...")
        res = analyze_single_galaxy(name, fname)
        if res:
            res["type"] = gtype
            results.append(res)
            print(f"  -> Deviation: {res['deviation']:+.1f}% | RMS: {res['rms']:.1f} km/s")

    # ==========================================
    # PLOTTING: SPECTRUM CHART
    # ==========================================
    if not results:
        print("No results to plot. Please check data files.")
        return

    # Sort results by deviation
    results.sort(key=lambda x: x["deviation"])
    
    names = [r["name"] + f"\n({r['type']})" for r in results]
    devs = [r["deviation"] for r in results]
    
    # Assign colors based on zones
    colors = []
    for d in devs:
        if d < 7.0: colors.append('#2ecc71')       # Green (Standard)
        elif d < 25.0: colors.append('#f39c12')    # Orange (Intermediate)
        else: colors.append('#e74c3c')             # Red (Signal)

    plt.figure(figsize=(12, 8))
    bars = plt.barh(names, devs, color=colors, edgecolor='black', alpha=0.8)
    
    # Add text labels
    for bar, d in zip(bars, devs):
        width = bar.get_width()
        x_pos = width + 1 if width >= 0 else width - 5
        plt.text(x_pos, bar.get_y() + bar.get_height()/2, 
                 f"{d:+.1f}%", va='center', fontweight='bold', fontsize=12)

    # Zones
    plt.axvline(x=7.0, color='green', linestyle='--', alpha=0.5)
    plt.text(3.5, -0.8, "STANDARD\nZONE", color='green', fontweight='bold', ha='center', fontsize=10)
    
    plt.axvline(x=30.0, color='red', linestyle='--', alpha=0.5)
    plt.text(35, -0.8, "30% SIGNAL\nZONE", color='darkred', fontweight='bold', ha='center', fontsize=10)
    
    plt.title("QIC-S Deff Emergence Spectrum\nZero-Parameter: a0=1.23e-10, M/L=0.5 (ALL FIXED)", fontsize=16, fontweight='bold')
    plt.xlabel("Model Deviation from Observation [%] (Positive = Model is faster)", fontsize=12)
    plt.grid(axis='x', alpha=0.3)
    plt.xlim(0, 50)
    
    plt.tight_layout()
    plt.savefig("QICS_Spectrum_AutoGenerated.png")
    print("\nSpectrum chart saved as 'QICS_Spectrum_AutoGenerated.png'")
    # plt.show() # Uncomment if running in a windowed environment

if __name__ == "__main__":
    main()
