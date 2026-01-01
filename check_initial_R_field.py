#!/usr/bin/env python3
"""Check initial R field structure across grids"""
import pandas as pd
import numpy as np

grids = {
    '64×64': 'output/20251219_031821_defect_localization_64x64/N_100/observables.csv',
    '128×128': 'output/20251219_031829_defect_localization/N_100/observables.csv',
    '256×256': 'output/20251219_031834_defect_localization_256x256/N_100/observables.csv',
}

print("Initial R field statistics:")
print()

for label, path in grids.items():
    df = pd.read_csv(path)
    print(f"{label}:")
    print(f"  R_avg(0) = {df['R_avg'].iloc[0]:.6f}")
    print(f"  R_min(0) = {df['R_min'].iloc[0]:.6f}")
    print(f"  R_max(0) = {df['R_max'].iloc[0]:.6f}")
    print(f"  R_var(0) = {df['R_var'].iloc[0]:.6f}")
    print()

