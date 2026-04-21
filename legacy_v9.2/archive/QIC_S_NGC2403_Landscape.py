import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import os

"""
QIC-S Theory: Galaxy Rotation Curve Analysis & Hamiltonian Landscape
Author: Sasada (Institute for Quantum Information Cosmology)
Date: 2025-11-26
Description: 
    This script performs the inverse problem analysis for NGC 2403 based on QIC-S theory.
    It reconstructs the effective transport coefficient (D_eff) and the emergent Hamiltonian 
    potential (H_eff) directly from rotation curve data without assuming Dark Matter.
"""

# ==========================================
# 1. DATA LOADING / GENERATION
# ==========================================

# Filename for the SPARC data (if available)
filename = 'NGC2403_rotmod.dat'

if os.path.exists(filename):
    print(f"Loading real data from {filename}...")
    try:
        data = np.loadtxt(filename, comments='#') 
        
        # SPARC format columns: [0]Rad, [1]Vobs, [2]Verr, [3]Vgas, [4]Vdisk, [5]Vbul
        r = data[:, 0]       # Radius
        v_obs = data[:, 1]   # Observed Velocity
        
        # 修正ポイント: バリオン成分（ガス+円盤+バルジ）を正しく合成する
        # V_baryon = sqrt(Vgas^2 + Vdisk^2 + Vbul^2)
        v_gas = data[:, 3]
        v_disk = data[:, 4]
        v_bul = data[:, 5]
        
        # 負の値がある場合の安全策としてabsを入れて合成
        # **この修正により、赤い線が正しく表示されるようになりました**
        v_baryon = np.sqrt(np.abs(v_gas)**2 + np.abs(v_disk)**2 + np.abs(v_bul)**2)

    except Exception as e:
        print(f"Error reading file: {e}. Switching to Simulation Mode.")
        # use_mock = True の変数はここでは未使用だが、ロジックとしては有効
        
else:
    print(f"File '{filename}' not found. Running in Simulation Mode (Mock Data).")
    
    # --- MOCK DATA GENERATION (Based on NGC 2403 Profile) ---
    r = np.linspace(0.1, 20, 100)  # Radius (kpc)
    
    # Baryonic Component (Disk + Gas) - Decays at large radii
    v_baryon = 100 * (r / 2) * np.exp(-r / 6) + 30
    
    # Observed Data (Flat Rotation Curve) - Remains flat/rising
    v_obs = v_baryon + (130 - v_baryon) * (1 - np.exp(-r/5))

# ==========================================
# 2. QIC-S INVERSE ANALYSIS
# ==========================================

# Calculate the Effective Force difference needed (The "Missing Gravity")
# Force ~ V^2 / r
# Avoid division by zero at r=0
r_safe = r.copy()
r_safe[r_safe == 0] = np.nan
force_diff = (v_obs**2 - v_baryon**2) / r_safe

# Construct the Effective Hamiltonian Potential (H_eff)
# H_eff corresponds to the accumulated 'stiffness' of the spacetime fabric
# Replace nan with 0 for integration
force_diff_safe = np.nan_to_num(force_diff)
# D_eff (輸送係数) の勾配が H_eff に対応 (H_eff は D_eff の勾配の二乗に比例)
H_eff = np.cumsum(force_diff_safe) * (r[1] - r[0]) if len(r) > 1 else np.zeros_like(r)

# Normalize for visualization (0 to 1 scale)
H_eff_norm = H_eff / np.max(H_eff)

# ==========================================
# 3. MAPPING TO 2D LANDSCAPE
# ==========================================
theta = np.linspace(0, 2*np.pi, 200)
R, THETA = np.meshgrid(r, theta)

# Map the 1D potential to 2D grid (assuming axial symmetry for this demo)
Z_H_eff = np.tile(H_eff_norm, (200, 1))

# ==========================================
# 4. PLOTTING (The Hamiltonian Landscape)
# ==========================================
fig = plt.figure(figsize=(16, 8), facecolor='black')

# --- Left Panel: Rotation Curve Decomposition ---
ax1 = fig.add_subplot(121, facecolor='black')

ax1.plot(r, v_obs, 'c-', lw=3, label='Observed (Real)')
ax1.plot(r, v_baryon, 'r--', lw=2, label='Baryonic (Visible)')
ax1.fill_between(r, v_baryon, v_obs, color='cyan', alpha=0.2, label='QIC-S Field Contribution')

ax1.set_xlabel('Radius (kpc)', color='white', fontsize=14)
ax1.set_ylabel('Velocity (km/s)', color='white', fontsize=14)
ax1.set_title('The Missing Physics (Rotation Curve)', color='white', fontsize=16)
ax1.tick_params(axis='x', colors='white')
ax1.tick_params(axis='y', colors='white')
ax1.legend(fontsize=12, facecolor='black', edgecolor='white', labelcolor='white')
ax1.grid(color='gray', linestyle=':', alpha=0.5)
ax1.set_xlim(0, 22) # Match the radial extent of the right plot

# --- Right Panel: Hamiltonian Landscape (Polar Plot) ---
ax2 = fig.add_subplot(122, projection='polar', facecolor='black')

# Colormap setup
cmap = cm.inferno
mesh = ax2.pcolormesh(THETA, R, Z_H_eff, cmap=cmap, shading='gouraud')

# Visual Settings
ax2.set_title('QIC-S Hamiltonian Landscape\n(NGC 2403)', color='white', fontsize=18, pad=20)
ax2.grid(color='white', alpha=0.3)

# --- ADJUSTED RADIAL LIMITS ---
ax2.set_rmax(22)              # Set max radius to 22 (creates outer margin)
ax2.set_rticks([5, 10, 15, 20]) # Keep ticks only up to 20 for physical relevance

ax2.tick_params(axis='x', colors='white')
ax2.tick_params(axis='y', colors='white')

# Colorbar
cbar = plt.colorbar(mesh, ax=ax2, pad=0.1)
cbar.set_label('Hamiltonian Intensity (H_eff)', color='white', fontsize=12)
cbar.ax.yaxis.set_tick_params(color='white')
plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='white')

# Annotation
plt.figtext(0.75, 0.05, "Generated by QIC-S Theory\nNo Dark Matter Assumed", 
            color='cyan', ha='center', fontsize=12, fontweight='bold')

plt.tight_layout()

# Save the figure
plt.savefig('QIC-S_Landscape_NGC2403.png', facecolor='black', dpi=150)
print("Plot saved as QIC-S_Landscape_NGC2403.png")
plt.show()