"""
Basic Kuramoto model synchronization example.

This script demonstrates:
1. Creating a Kuramoto model
2. Running simulations
3. Computing order parameters
4. Visualizing results
5. Finding critical coupling
"""

import numpy as np
import matplotlib.pyplot as plt

# These imports will work once the package is implemented
# from kuramoto import KuramotoModel
# from kuramoto.distributions import LorentzianDistribution
# from kuramoto.analysis import BifurcationAnalysis
# from kuramoto.visualization import plot_order_parameter, plot_phases


def example_basic_synchronization():
    """
    Basic example of Kuramoto synchronization.

    Shows transition from incoherent to partially synchronized state
    as coupling strength increases past critical value.
    """
    # Set parameters
    N = 100  # Number of oscillators
    gamma = 1.0  # Lorentzian width
    Kc = 2 * gamma  # Theoretical critical coupling

    print(f"Kuramoto Model Synchronization Example")
    print(f"N = {N} oscillators")
    print(f"Lorentzian distribution: γ = {gamma}")
    print(f"Theoretical Kc = {Kc}\n")

    # Test three regimes: subcritical, critical, supercritical
    K_values = [0.5 * Kc, Kc, 2 * Kc]
    labels = ['Subcritical (K < Kc)', 'Critical (K = Kc)', 'Supercritical (K > Kc)']

    # Would be replaced with actual implementation:
    """
    fig, axes = plt.subplots(3, 2, figsize=(12, 10))

    for idx, (K, label) in enumerate(zip(K_values, labels)):
        print(f"Simulating {label}: K = {K:.2f}")

        # Create model
        dist = LorentzianDistribution(center=0, width=gamma)
        model = KuramotoModel(
            N=N,
            coupling=K,
            frequencies=dist,
            initial_phases=None  # Random initial conditions
        )

        # Simulate
        t_span = (0, 50)
        solution = model.evolve(t_span, solver='rk45')

        # Extract results
        t = solution['t']
        R = solution['R']
        Psi = solution['Psi']
        phases = solution['phases']

        # Theoretical steady state (from Ott-Antonsen)
        R_theory = dist.steady_state_order_parameter(K)

        # Plot time series
        ax1 = axes[idx, 0]
        ax1.plot(t, R, 'b-', label='R(t)')
        ax1.axhline(R_theory, color='r', linestyle='--', label=f'Theory: R={R_theory:.3f}')
        ax1.set_xlabel('Time (t)')
        ax1.set_ylabel('Order Parameter R')
        ax1.set_title(f'{label}')
        ax1.set_ylim([0, 1])
        ax1.legend()
        ax1.grid(True, alpha=0.3)

        # Plot final phase distribution on unit circle
        ax2 = axes[idx, 1]
        ax2 = plt.subplot(3, 2, 2*idx + 2, projection='polar')
        final_phases = phases[-1, :]

        # Plot oscillators
        ax2.scatter(final_phases, np.ones_like(final_phases),
                   c=model.frequencies, cmap='coolwarm', s=20, alpha=0.6)

        # Plot order parameter vector
        final_R, final_Psi = R[-1], Psi[-1]
        ax2.arrow(0, 0, final_Psi, final_R, color='red',
                 width=0.05, head_width=0.1, head_length=0.05,
                 label=f'R={final_R:.3f}')

        ax2.set_ylim([0, 1])
        ax2.set_title(f'Phase Distribution (t={t_span[1]})')

        print(f"  Final R = {final_R:.3f} (Theory: {R_theory:.3f})")
        print(f"  Steady-state reached: {np.std(R[-100:]) < 0.01}\n")

    plt.tight_layout()
    plt.savefig('synchronization_regimes.png', dpi=150)
    plt.show()
    """

    # DEMONSTRATION OUTPUT (actual visualization requires matplotlib)
    print("Simulation would show:")
    print("- Subcritical: R → 0 (incoherent)")
    print("- Critical: R fluctuates near 0")
    print("- Supercritical: R → finite value (partial synchronization)")


def example_bifurcation_diagram():
    """
    Create bifurcation diagram showing R vs K.

    Demonstrates the continuous phase transition at Kc.
    """
    print("\nBifurcation Diagram")
    print("-" * 40)

    # Would be replaced with actual implementation:
    """
    # Parameters
    N = 500
    gamma = 1.0
    Kc = 2 * gamma
    K_range = np.linspace(0, 4 * Kc, 50)

    # Storage for results
    R_steady = []
    R_theory = []

    print(f"Scanning K from 0 to {4*Kc:.1f}")

    for K in K_range:
        # Create and simulate model
        dist = LorentzianDistribution(center=0, width=gamma)
        model = KuramotoModel(N=N, coupling=K, frequencies=dist)

        # Simulate to steady state
        solution = model.evolve((0, 100), solver='rk45')

        # Extract steady-state R (average last 20% of simulation)
        R_final = np.mean(solution['R'][-int(len(solution['R'])*0.2):])
        R_steady.append(R_final)

        # Theoretical value
        R_th = dist.steady_state_order_parameter(K)
        R_theory.append(R_th)

    # Create bifurcation plot
    plt.figure(figsize=(10, 6))

    plt.plot(K_range, R_steady, 'bo', markersize=4, label='Simulation')
    plt.plot(K_range, R_theory, 'r-', linewidth=2, label='Ott-Antonsen theory')
    plt.axvline(Kc, color='gray', linestyle='--', alpha=0.5, label=f'Kc = {Kc}')

    plt.xlabel('Coupling Strength K')
    plt.ylabel('Order Parameter R')
    plt.title('Bifurcation Diagram: Kuramoto Synchronization Transition')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.xlim([0, max(K_range)])
    plt.ylim([0, 1])

    # Add text annotation
    plt.text(Kc * 0.5, 0.5, 'Incoherent\nPhase', ha='center', fontsize=12)
    plt.text(Kc * 2.5, 0.5, 'Partially\nSynchronized', ha='center', fontsize=12)

    plt.savefig('bifurcation_diagram.png', dpi=150)
    plt.show()

    # Find numerical critical coupling
    analyzer = BifurcationAnalysis(model)
    Kc_numerical = analyzer.find_critical_coupling()
    print(f"\nCritical coupling:")
    print(f"  Theoretical: Kc = {Kc:.3f}")
    print(f"  Numerical:   Kc = {Kc_numerical:.3f}")
    """

    # DEMONSTRATION calculation
    gamma = 1.0
    Kc = 2 * gamma
    K_test = 3.0
    R_theory = np.sqrt(1 - (Kc/K_test)**2) if K_test > Kc else 0

    print(f"γ = {gamma}, Kc = {Kc}")
    print(f"For K = {K_test}: R_theory = {R_theory:.3f}")
    print("\nBifurcation shows continuous transition at Kc")
    print("R = 0 for K < Kc (second-order phase transition)")
    print("R ∝ √(K - Kc) near critical point")


def example_network_coupling():
    """
    Example with network-coupled oscillators.

    Shows how coupling topology affects synchronization.
    """
    print("\nNetwork-Coupled Oscillators")
    print("-" * 40)

    # Would be replaced with actual implementation:
    """
    import networkx as nx
    from kuramoto import NetworkCoupling

    # Create different network topologies
    N = 50
    networks = {
        'All-to-all': nx.complete_graph(N),
        'Ring': nx.cycle_graph(N),
        'Small-world': nx.watts_strogatz_graph(N, 4, 0.3),
        'Scale-free': nx.barabasi_albert_graph(N, 2)
    }

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()

    for idx, (name, G) in enumerate(networks.items()):
        print(f"Testing {name} network...")

        # Get adjacency matrix
        A = nx.adjacency_matrix(G).todense()

        # Create model with network coupling
        coupling = NetworkCoupling(A, strength=5.0)
        model = KuramotoModel(
            N=N,
            coupling=coupling,
            frequencies='gaussian'  # Gaussian distributed frequencies
        )

        # Simulate
        solution = model.evolve((0, 50))

        # Plot results
        ax = axes[idx]
        ax.plot(solution['t'], solution['R'])
        ax.set_xlabel('Time (t)')
        ax.set_ylabel('Order Parameter R')
        ax.set_title(f'{name} Network (N={N})')
        ax.set_ylim([0, 1])
        ax.grid(True, alpha=0.3)

        final_R = solution['R'][-1]
        print(f"  Final R = {final_R:.3f}")

    plt.tight_layout()
    plt.savefig('network_synchronization.png', dpi=150)
    plt.show()
    """

    print("Network topology strongly affects synchronization:")
    print("- All-to-all: Easiest to synchronize (mean-field)")
    print("- Ring: Hardest to synchronize (local coupling only)")
    print("- Small-world: Intermediate (shortcuts help)")
    print("- Scale-free: Hubs facilitate synchronization")


def example_custom_analysis():
    """
    Example of custom analysis and metrics.

    Shows how to extract additional information from simulations.
    """
    print("\nCustom Analysis Example")
    print("-" * 40)

    # DEMONSTRATION for actual implementation
    print("Analysis capabilities:")
    print("1. Phase coherence: ρ = |⟨e^(iθ)⟩|")
    print("2. Phase diffusion: D = ⟨(θ - ⟨θ⟩)²⟩")
    print("3. Cluster detection: Find phase-locked groups")
    print("4. Stability analysis: Lyapunov exponents")
    print("5. Spectral analysis: Frequency entrainment")

    # Example calculation
    N = 100
    phases = np.random.uniform(0, 2*np.pi, N)

    # Phase coherence
    z = np.mean(np.exp(1j * phases))
    R = np.abs(z)

    # Phase variance
    mean_phase = np.angle(z)
    phase_variance = np.var(np.angle(np.exp(1j * (phases - mean_phase))))

    print(f"\nExample metrics for N={N} oscillators:")
    print(f"  Phase coherence R = {R:.3f}")
    print(f"  Phase variance = {phase_variance:.3f}")
    print(f"  Effective coupling = K_eff = K·R")


if __name__ == "__main__":
    # Run examples
    example_basic_synchronization()
    example_bifurcation_diagram()
    example_network_coupling()
    example_custom_analysis()

    print("\n" + "="*50)
    print("Examples complete!")
    print("This demonstrates the architecture's capabilities:")
    print("- Flexible model configuration")
    print("- Various frequency distributions")
    print("- Network coupling options")
    print("- Comprehensive analysis tools")
    print("- Extensible design for future enhancements")