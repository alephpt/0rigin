#!/usr/bin/env python
"""
Run all field theory demonstrations and generate validation report.
"""

import sys
import os
import subprocess
import time

def run_demo(script_name, description):
    """Run a demo script and capture results."""
    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"Script: {script_name}")
    print('='*60)

    try:
        result = subprocess.run(
            [sys.executable, script_name],
            capture_output=True,
            text=True,
            timeout=30
        )

        # Print first 100 lines of output
        lines = result.stdout.split('\n')
        for i, line in enumerate(lines[:100]):
            print(line)
        if len(lines) > 100:
            print(f"\n... ({len(lines)-100} more lines) ...")

        if result.returncode != 0:
            print(f"\nError output:\n{result.stderr[:500]}")
            return False
        return True
    except subprocess.TimeoutExpired:
        print(f"Timeout: Script took more than 30 seconds")
        return False
    except Exception as e:
        print(f"Failed to run: {e}")
        return False

def main():
    """Run all demonstrations."""
    print("="*60)
    print("FIELD THEORY PROTOTYPES - FULL DEMONSTRATION")
    print("="*60)

    # Change to examples directory
    os.chdir('examples/field_theory')

    demos = [
        ('validate_prototypes.py', 'Prototype Validation Tests'),
        ('hamiltonian_demo.py', 'Hamiltonian Kuramoto Demonstration'),
        ('grid_test.py', 'Spatial Grid Infrastructure Test'),
        ('pde_integration_test.py', 'PDE Solver Integration Test'),
    ]

    results = {}

    for script, description in demos:
        if os.path.exists(script):
            results[script] = run_demo(script, description)
        else:
            print(f"Script not found: {script}")
            results[script] = False

    # Summary
    print("\n" + "="*60)
    print("DEMONSTRATION SUMMARY")
    print("="*60)

    for script, success in results.items():
        status = "✓ PASS" if success else "✗ FAIL"
        print(f"  {status}: {script}")

    passed = sum(1 for v in results.values() if v)
    total = len(results)

    print(f"\nTotal: {passed}/{total} demonstrations completed successfully")

    if passed == total:
        print("\n✓ All field theory prototypes validated successfully!")
    else:
        print(f"\n✗ {total-passed} demonstration(s) failed")

    return passed == total

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)