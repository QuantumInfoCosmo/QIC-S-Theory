"""
QIC-S Universal Scaling Law Plotter (Ver 8.3)
Author: Yoshiaki Sasada
Description: Visualizes the renormalization group flow (scaling law) of the effective 
             transport coefficient D_eff across galactic and filament scales.
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob
from scipy.stats import linregress

# ==========================================
# Configuration
# ==========================================
DATA_DIR = 'data'
CSV_FILE = 'results/QIC_S_Result_N170.csv'
OUTPUT_IMAGE = 'Fig4_Universal_Scaling_Law.png'
MIN_POINTS = 5

# Filament Data (From Tudorache et al., 2025)
# Format: (Label, Radius_kpc, Velocity_km_s)
FILAMENT_DATA = [
    ("Filament Core", 50, 110),       # 0.05 Mpc
    ("HI Structure", 1700, 110),      # 1.7 Mpc
    ("Full Filament", 15000, 110)     # 15 Mpc
]

def calculate_galaxy_properties(filepath):
    """Calculates R_max, V_last, and D_eff for a galaxy."""
    try:
        df = pd.read_csv(filepath, sep=r'\s+', engine='python', header=None, comment='#')
        data = df.apply(pd.to_numeric, errors='coerce').dropna()
        r = data.iloc[:, 0].values
        v = data.iloc[:, 1].values
        
        # Physical filter
        mask = (r > 0) & (v > 0)
        r, v = r[mask], v[mask]
        
        if len(r) < MIN_POINTS: return None
        
        # Representative Scale and Dynamics
        # Using the last observed point (largest scale)
        r_max = np.max(r)
        v_at_rmax = v[np.argmax(r)] # Velocity at max radius
        
        # D_eff ≈ R * v
        d_eff = r_max * v_at_rmax
        
        return r_max, d_eff
    except:
        return None

def main():
    # 1. Load Galaxy Phase Info
    if not os.path.exists(CSV_FILE):
        print(f"Error: {CSV_FILE} not found. Please run phase_analysis.py first.")
        return
        
    m_df = pd.read_csv(CSV_FILE)
    galaxy_map = dict(zip(m_df['Galaxy'], m_df['M']))
    
    # 2. Collect Galaxy Scaling Data
    scaling_data = []
    files = glob.glob(os.path.join(DATA_DIR, "*_rotmod.dat"))
    
    print("Processing galaxy data...")
    for f in files:
        fname = os.path.basename(f)
        if fname in galaxy_map:
            props = calculate_galaxy_properties(f)
            if props:
                r, d = props
                m = galaxy_map[fname]
                scaling_data.append({'R': r, 'D_eff': d, 'M': m, 'Type': 'Galaxy'})

    df = pd.DataFrame(scaling_data)
    
    # 3. Add Filament Data
    filament_points = []
    for label, r, v in FILAMENT_DATA:
        d_eff = r * v
        filament_points.append({'Label': label, 'R': r, 'D_eff': d_eff})
        
    # 4. Visualization (Log-Log Plot)
    plt.figure(figsize=(12, 8))
    sns.set_style("whitegrid")
    
    # Plot Galaxies (Colored by Phase M)
    sc = plt.scatter(df['R'], df['D_eff'], c=df['M'], cmap='coolwarm', 
                     s=50, alpha=0.7, edgecolors='k', label='SPARC Galaxies (N=170)')
    plt.colorbar(sc, label='Phase Metric M (Order < 0.5 < Chaos)')
    
    # Plot Filament Points
    for pt in filament_points:
        plt.scatter(pt['R'], pt['D_eff'], color='gold', marker='*', s=300, edgecolors='black', zorder=10)
        plt.text(pt['R']*1.1, pt['D_eff'], f"  {pt['Label']}\n  ({pt['R']/1000} Mpc)", 
                 fontsize=11, fontweight='bold', color='darkorange', va='center')

    # 5. Fit Scaling Law (Power Law)
    # Combine all data for fitting
    all_r = np.concatenate([df['R'].values, [p['R'] for p in filament_points]])
    all_d = np.concatenate([df['D_eff'].values, [p['D_eff'] for p in filament_points]])
    
    slope, intercept, r_value, p_value, std_err = linregress(np.log10(all_r), np.log10(all_d))
    
    # Plot Trend Line
    x_range = np.linspace(min(all_r)*0.8, max(all_r)*1.5, 100)
    y_fit = 10**intercept * x_range**slope
    plt.plot(x_range, y_fit, 'k--', linewidth=1.5, alpha=0.8, 
             label=f'Universal Scaling Law\n(Slope = {slope:.2f}, R²={r_value**2:.3f})')

    # Formatting
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Characteristic Scale R [kpc]', fontsize=14)
    plt.ylabel('Effective Transport Coefficient D_eff [kpc km/s]', fontsize=14)
    plt.title('Universal Scaling of Hamiltonian Dynamics: From Galaxies to Filaments', fontsize=18)
    plt.legend(fontsize=12, loc='upper left')
    plt.grid(True, which="both", ls="--", alpha=0.4)
    
    plt.savefig(OUTPUT_IMAGE, dpi=300, bbox_inches='tight')
    print(f"\n[Success] Scaling Law Plot generated: {OUTPUT_IMAGE}")
    print(f"Scaling Slope: {slope:.3f} (Theoretical expectation ≈ 1.0 for flat rotation)")

if __name__ == "__main__":
    main()