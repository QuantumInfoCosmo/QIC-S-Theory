# -*- coding: utf-8 -*-
"""
================================================================================
QIC-S Theory: Hamiltonian Landscape Analyzer
Version: 2.1 (Scientific Update - Fixed Metric)
Author: Yoshiaki Sasada
License: MIT

Description:
    This script analyzes galaxy rotation curves using the QIC-S theoretical framework.
    
    [UPDATE Ver 2.1 Fixed]: 
    - Aligned Phase Metric calculation with phase_analysis.py (v^2/r).
    - M = Var(log(v^2/r)).
    - Landscape visualization now represents the Log-Gradient field directly.

Usage:
    python qics_analyzer.py --file1 data/NGC0100_rotmod.dat --output Fig1.png
================================================================================
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import argparse
import os
import sys

# ==========================================
# 1. PHYSICS CONSTANTS
# ==========================================
A0 = 1.23e-10
ML_DISK = 0.5
ML_GAS = 1.0
ML_BULGE = 0.7
KPC_TO_M = 3.086e19
KMS_TO_MS = 1000.0

# ==========================================
# 2. COLOR SCHEME
# ==========================================
COLORS = {
    'observed_mature': 'cyan',      # Phase 5: Mature galaxy observation
    'observed_germinating': 'red',  # Phase 4: Germinating galaxy observation
    'prediction': 'white',          # QIC-S theoretical prediction
    'entropic_release': 'red',      # Energy excess region
    'background': 'black',
    'text': 'white',
    'title_mature': 'cyan',
    'title_germinating': 'orange'
}

CMAP_ORDER = 'viridis'
CMAP_CHAOS = 'inferno'

# ==========================================
# 3. QICS ANALYZER CLASS
# ==========================================
class QICSAnalyzer:
    def __init__(self, ml_disk=ML_DISK, ml_gas=ML_GAS, ml_bulge=ML_BULGE):
        self.ml_disk = ml_disk
        self.ml_gas = ml_gas
        self.ml_bulge = ml_bulge
        
    def load_sparc_file(self, filename):
        if not os.path.exists(filename):
            print(f"[ERROR] File not found: {filename}")
            return None
            
        try:
            data = np.loadtxt(filename, comments='#')
            r = data[:, 0]
            v_obs = data[:, 1]
            
            # Simple baryon estimation (if columns exist)
            v_gas = np.abs(data[:, 3]) if data.shape[1] > 3 else np.zeros_like(r)
            v_disk = np.abs(data[:, 4]) if data.shape[1] > 4 else np.zeros_like(r)
            v_bul = np.abs(data[:, 5]) if data.shape[1] > 5 else np.zeros_like(r)
            
            v_baryon = np.sqrt(
                (v_gas**2 * self.ml_gas) + 
                (v_disk**2 * self.ml_disk) + 
                (v_bul**2 * self.ml_bulge)
            )
            
            galaxy_name = os.path.basename(filename).replace('_rotmod.dat', '').replace('.dat', '')
            
            return {
                'name': galaxy_name,
                'radius': r,
                'v_obs': v_obs,
                'v_baryon': v_baryon,
                'filename': filename
            }
        except Exception as e:
            print(f"[ERROR] Failed to load {filename}: {e}")
            return None
    
    def compute_qics_prediction(self, r, v_baryon):
        g_bar = (v_baryon**2) / r * (KMS_TO_MS**2) / KPC_TO_M
        g_bar = np.where(g_bar < 1e-15, 1e-15, g_bar)
        x = g_bar / A0
        g_tot = g_bar / (1 - np.exp(-np.sqrt(x)))
        v_pred = np.sqrt(g_tot * r * KPC_TO_M) / KMS_TO_MS
        return v_pred
    
    def compute_hamiltonian_landscape(self, r, v_obs):
        """
        Compute Phase Metric M using strict definition from phase_analysis.py.
        Metric: Log-Variance of the Gradient (v^2/r).
        """
        # --- CORRECT METRIC CALCULATION ---
        with np.errstate(divide='ignore', invalid='ignore'):
            # Total field gradient (v_obs^2 / r)
            grad_H = (v_obs ** 2) / r
            grad_H = np.nan_to_num(grad_H, nan=0.0, posinf=0.0, neginf=0.0)
            
        epsilon = 1e-10
        # Log-Gradient (The quantity whose variance is M)
        log_grad = np.log(np.abs(grad_H) + epsilon)
        
        # Calculate Metric
        phase_metric = np.var(log_grad)
        
        # --- VISUALIZATION DATA ---
        # Visualize the Log-Gradient itself so the plot matches the metric visually.
        # Normalize for plotting contrast
        if np.max(np.abs(log_grad)) > 0:
            # Shift to make it positive for visualization or just normalize range
             landscape_data = (log_grad - np.min(log_grad)) / (np.max(log_grad) - np.min(log_grad) + epsilon)
        else:
             landscape_data = np.zeros_like(log_grad)

        return landscape_data, phase_metric
    
    def classify_phase(self, phase_metric):
        PHYSICAL_THRESHOLD = 0.5 
        if phase_metric >= PHYSICAL_THRESHOLD:
            return 4, "Phase 4: Germinating (Chaos)"
        else:
            return 5, "Phase 5: Mature (Order)"
    
    def analyze(self, data):
        if data is None: return None
        r = data['radius']
        v_obs = data['v_obs']
        v_baryon = data['v_baryon']
        
        v_pred = self.compute_qics_prediction(r, v_baryon)
        landscape_data, phase_metric = self.compute_hamiltonian_landscape(r, v_obs)
        phase, phase_label = self.classify_phase(phase_metric)
        energy_excess = np.maximum(0, v_obs - v_pred)
        
        return {
            **data,
            'v_pred': v_pred,
            'landscape_data': landscape_data, 
            'phase_metric': phase_metric,
            'phase': phase,
            'phase_label': phase_label,
            'energy_excess': energy_excess
        }

# ==========================================
# 4. VISUALIZATION FUNCTIONS
# ==========================================
def plot_single_galaxy(results, ax_curve, ax_landscape):
    if results is None: return
    
    r = results['radius']
    v_obs = results['v_obs']
    v_pred = results['v_pred']
    landscape_data = results['landscape_data']
    
    phase = results['phase']
    
    if phase == 5:
        obs_color = COLORS['observed_mature']
        title_color = COLORS['title_mature']
        cmap = CMAP_ORDER
        landscape_title = "Hamiltonian Log-Gradient (Order)"
        show_excess = False
    else:
        obs_color = COLORS['observed_germinating']
        title_color = COLORS['title_germinating']
        cmap = CMAP_CHAOS
        landscape_title = "Hamiltonian Log-Gradient (Chaos)"
        show_excess = True
    
    # --- Curve ---
    ax_curve.set_facecolor(COLORS['background'])
    ax_curve.plot(r, v_obs, color=obs_color, lw=2, label='Observed')
    ax_curve.plot(r, v_pred, 'w--', lw=1.5, alpha=0.9, label='QIC-S Pred')
    
    if show_excess:
        ax_curve.fill_between(r, v_pred, v_obs, where=(v_obs > v_pred),
                              color=COLORS['entropic_release'], alpha=0.3)
    
    title = f"{results['name']} (Phase {phase})"
    ax_curve.set_title(title, color=title_color, fontsize=14, fontweight='bold')
    ax_curve.set_xlabel("Radius (kpc)", color=COLORS['text'])
    ax_curve.set_ylabel("Velocity (km/s)", color=COLORS['text'])
    ax_curve.legend(facecolor='black', labelcolor='white', loc='lower right')
    ax_curve.tick_params(colors=COLORS['text'])
    ax_curve.grid(alpha=0.3, color='gray')
    
    # --- Landscape ---
    ax_landscape.set_facecolor(COLORS['background'])
    theta = np.linspace(0, 2*np.pi, 100)
    R, THETA = np.meshgrid(r, theta)
    Z = np.tile(landscape_data, (100, 1))
    
    ax_landscape.pcolormesh(THETA, R, Z, cmap=cmap, shading='auto')
    ax_landscape.set_title(landscape_title, color=COLORS['text'], fontsize=12, pad=10)
    ax_landscape.set_xticklabels([])
    ax_landscape.tick_params(colors=COLORS['text'])
    
    ax_landscape.text(0, 0, f"M = {results['phase_metric']:.3f}", 
                     color='white', ha='center', va='center', 
                     fontsize=12, fontweight='bold',
                     bbox=dict(facecolor='black', alpha=0.6, edgecolor='white'))

def create_comparison_figure(results_list, output_filename):
    n_galaxies = len(results_list)
    fig = plt.figure(figsize=(10*n_galaxies, 12), facecolor=COLORS['background'])
    gs = GridSpec(2, n_galaxies, figure=fig, hspace=0.25, wspace=0.15)
    
    for i, results in enumerate(results_list):
        if results is None: continue
        ax_curve = fig.add_subplot(gs[0, i], facecolor=COLORS['background'])
        ax_landscape = fig.add_subplot(gs[1, i], projection='polar', facecolor=COLORS['background'])
        plot_single_galaxy(results, ax_curve, ax_landscape)
    
    plt.suptitle("QIC-S Theory: Hamiltonian Landscape", color=COLORS['text'], fontsize=18, fontweight='bold', y=0.98)
    plt.savefig(output_filename, dpi=150, facecolor=COLORS['background'], bbox_inches='tight', pad_inches=0.2)
    print(f"[SUCCESS] Figure saved: {output_filename}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--file1', required=True)
    parser.add_argument('--output', default='Fig1.png')
    args = parser.parse_args()
    
    analyzer = QICSAnalyzer()
    results = analyzer.analyze(analyzer.load_sparc_file(args.file1))
    
    if results:
        create_comparison_figure([results], args.output)

if __name__ == "__main__":
    main()