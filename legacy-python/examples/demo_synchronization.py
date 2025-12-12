"""
Kuramoto model synchronization demonstration.

This script demonstrates the complete functionality of the Kuramoto
model implementation, including:
1. Multiple frequency distributions
2. Order parameter analysis
3. Synchronization transitions
4. Visualization tools
"""

import sys
sys.path.insert(0, '/home/persist/neotec/0rigin/src')

import numpy as np
import matplotlib.pyplot as plt

from kuramoto import KuramotoModel
from kuramoto.distributions import (
    LorentzianDistribution,
    GaussianDistribution,
    UniformDistribution
)
from kuramoto.analysis import OrderParameter
from kuramoto.visualization import (
    plot_phases,
    plot_order_parameter,
    plot_bifurcation_diagram
)


def demo_lorentzian_synchronization():
    """
    Demonstrate synchronization with Lorentzian distribution.

    Shows the three regimes: subcritical, critical, and supercritical.
    """
    print("=" * 60)
    print("DEMO 1: Lorentzian Distribution Synchronization")
    print("=" * 60)

    N = 100
    gamma = 1.0
    Kc = 2 * gamma

    print(f"N = {N} oscillators")
    print(f"Lorentzian: γ = {gamma}, Kc = {Kc:.2f}\n")

    # Test three coupling values
    K_values = [0.5 * Kc, Kc, 2 * Kc]
    labels = ['Subcritical', 'Critical', 'Supercritical']

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    for idx, (K, label) in enumerate(zip(K_values, labels)):
        print(f"Simulating {label}: K = {K:.2f}")

        # Create model
        dist = LorentzianDistribution(center=0, width=gamma)
        model = KuramotoModel(N=N, coupling=K, frequencies=dist)

        # Simulate
        solution = model.evolve((0, 50), solver='rk45')

        # Compute theoretical steady state
        R_theory = dist.steady_state_order_parameter(K)

        # Plot
        ax = axes[idx]
        ax.plot(solution['t'], solution['R'], 'b-', linewidth=2)
        ax.axhline(R_theory, color='r', linestyle='--',
                   label=f'Theory: {R_theory:.3f}')
        ax.set_xlabel('Time (t)')
        ax.set_ylabel('Order Parameter R')
        ax.set_title(f'{label}\nK = {K:.2f}')
        ax.set_ylim([0, 1])
        ax.legend()
        ax.grid(True, alpha=0.3)

        final_R = np.mean(solution['R'][-100:])
        print(f"  Final R = {final_R:.3f} (Theory: {R_theory:.3f})")

    plt.tight_layout()
    plt.savefig('/home/persist/neotec/0rigin/lorentzian_regimes.png', dpi=150)
    print(f"  Saved: lorentzian_regimes.png\n")
    plt.close()


def demo_distribution_comparison():
    """
    Compare synchronization for different frequency distributions.
    """
    print("=" * 60)
    print("DEMO 2: Distribution Comparison")
    print("=" * 60)

    N = 150
    K = 3.0

    distributions = {
        'Lorentzian': LorentzianDistribution(center=0, width=1.0),
        'Gaussian': GaussianDistribution(mean=0, std=1.0),
        'Uniform': UniformDistribution(low=-1.5, high=1.5),
    }

    print(f"N = {N} oscillators, K = {K:.2f}\n")

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    for idx, (name, dist) in enumerate(distributions.items()):
        print(f"Simulating {name} distribution...")

        # Create model
        model = KuramotoModel(N=N, coupling=K, frequencies=dist)

        # Simulate
        solution = model.evolve((0, 50), solver='rk45')

        # Plot
        ax = axes[idx]
        ax.plot(solution['t'], solution['R'], linewidth=2)
        ax.set_xlabel('Time (t)')
        ax.set_ylabel('Order Parameter R')
        ax.set_title(name)
        ax.set_ylim([0, 1])
        ax.grid(True, alpha=0.3)

        final_R = np.mean(solution['R'][-100:])
        print(f"  Final R = {final_R:.3f}")

        # Get critical coupling if available
        Kc = dist.critical_coupling()
        if Kc is not None:
            print(f"  Critical coupling Kc = {Kc:.3f}")

    plt.tight_layout()
    plt.savefig('/home/persist/neotec/0rigin/distribution_comparison.png', dpi=150)
    print(f"  Saved: distribution_comparison.png\n")
    plt.close()


def demo_bifurcation_diagram():
    """
    Create bifurcation diagram showing R vs K transition.
    """
    print("=" * 60)
    print("DEMO 3: Bifurcation Diagram")
    print("=" * 60)

    N = 200
    gamma = 1.0
    Kc = 2 * gamma

    # Scan coupling values
    K_values = np.linspace(0, 4 * Kc, 30)
    R_steady = []
    R_theory = []

    print(f"N = {N} oscillators")
    print(f"Scanning K from 0 to {4*Kc:.1f}...\n")

    dist = LorentzianDistribution(center=0, width=gamma)

    for K in K_values:
        # Create and simulate
        model = KuramotoModel(N=N, coupling=K, frequencies=dist)
        solution = model.evolve((0, 100), solver='rk45')

        # Extract steady-state R
        R_final = np.mean(solution['R'][-200:])
        R_steady.append(R_final)

        # Theoretical value
        R_th = dist.steady_state_order_parameter(K)
        R_theory.append(R_th)

    # Plot
    fig, ax = plt.subplots(figsize=(10, 6))
    plot_bifurcation_diagram(
        K_values, np.array(R_steady),
        R_theory=np.array(R_theory),
        Kc=Kc,
        ax=ax,
        title='Kuramoto Bifurcation Diagram (Lorentzian)'
    )

    plt.tight_layout()
    plt.savefig('/home/persist/neotec/0rigin/bifurcation_diagram.png', dpi=150)
    print(f"Saved: bifurcation_diagram.png\n")
    plt.close()


def demo_phase_visualization():
    """
    Visualize phase distributions on unit circle.
    """
    print("=" * 60)
    print("DEMO 4: Phase Visualization")
    print("=" * 60)

    N = 100
    gamma = 1.0

    K_values = [1.0, 2.0, 4.0]
    labels = ['Incoherent', 'Critical', 'Synchronized']

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    for idx, (K, label) in enumerate(zip(K_values, labels)):
        print(f"Simulating {label} state: K = {K:.2f}")

        dist = LorentzianDistribution(center=0, width=gamma)
        model = KuramotoModel(N=N, coupling=K, frequencies=dist)

        # Evolve to steady state
        solution = model.evolve((0, 50), solver='rk45')
        final_phases = solution['phases'][-1, :]

        # Plot phases
        ax = axes[idx]
        plot_phases(final_phases, model.frequencies, ax=ax, title=label)

    plt.tight_layout()
    plt.savefig('/home/persist/neotec/0rigin/phase_distributions.png', dpi=150)
    print(f"Saved: phase_distributions.png\n")
    plt.close()


def demo_order_parameter_analysis():
    """
    Demonstrate order parameter analysis tools.
    """
    print("=" * 60)
    print("DEMO 5: Order Parameter Analysis")
    print("=" * 60)

    N = 150
    K = 3.0
    gamma = 1.0

    print(f"N = {N} oscillators, K = {K:.2f}\n")

    # Create and simulate
    dist = LorentzianDistribution(center=0, width=gamma)
    model = KuramotoModel(N=N, coupling=K, frequencies=dist)
    solution = model.evolve((0, 50), solver='rk45')

    # Analyze with OrderParameter class
    op = OrderParameter(solution['phases'])

    print("Analysis Results:")
    print(f"  Mean amplitude: <R> = {op.mean_amplitude():.3f}")
    print(f"  Steady-state amplitude: R_∞ = {op.steady_state_amplitude():.3f}")
    print(f"  Synchronized: {op.is_synchronized(threshold=0.5)}")

    conv_time = op.convergence_time(t=solution['t'], threshold=0.01)
    if conv_time is not None:
        print(f"  Convergence time: t_conv ≈ {conv_time:.2f}")

    locked_frac = op.phase_locked_fraction(tolerance=0.1)
    print(f"  Phase-locked fraction: {locked_frac:.2f}")

    print()


def main():
    """Run all demonstrations."""
    print("\n")
    print("╔" + "═" * 58 + "╗")
    print("║" + " " * 58 + "║")
    print("║" + "  Kuramoto Model: Complete Implementation Demo".center(58) + "║")
    print("║" + " " * 58 + "║")
    print("╚" + "═" * 58 + "╝")
    print("\n")

    # Run demos
    demo_lorentzian_synchronization()
    demo_distribution_comparison()
    demo_bifurcation_diagram()
    demo_phase_visualization()
    demo_order_parameter_analysis()

    print("=" * 60)
    print("All demonstrations complete!")
    print("=" * 60)
    print("\nGenerated files:")
    print("  - lorentzian_regimes.png")
    print("  - distribution_comparison.png")
    print("  - bifurcation_diagram.png")
    print("  - phase_distributions.png")
    print("\nImplementation features:")
    print("  ✓ Multiple frequency distributions (Lorentzian, Gaussian, Uniform)")
    print("  ✓ Order parameter computation and analysis")
    print("  ✓ Multiple solvers (RK4, RK45, Euler)")
    print("  ✓ Comprehensive visualization tools")
    print("  ✓ Clean module organization")
    print("  ✓ Type hints and docstrings throughout")


if __name__ == "__main__":
    main()
