"""
Demonstration of Self-consistent Mean Field Theory (SMFT) system.

Shows the complete integration of discrete oscillators with continuous
mediator fields, illustrating the bridge between particle and field descriptions.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Add parent to path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from src.kuramoto.field_theory import SMFTSystem


def demo_basic_evolution():
    """Basic SMFT system evolution."""
    print("=" * 60)
    print("SMFT System: Basic Evolution")
    print("=" * 60)

    # Create system
    system = SMFTSystem(
        grid_shape=(50, 50),
        N_oscillators=100,
        coupling='local',
        mediator_mass=10.0
    )

    print(f"\nSystem: {system}")
    print(f"Grid: {system.grid}")
    print(f"Oscillators: N={system.N}")

    # Evolve
    print("\nEvolving system...")
    solution = system.evolve(
        t_span=(0, 30),
        dt=0.01,
        diffusion_coeff=0.02,
        kernel_width=0.1,
        store_interval=100
    )

    # Analyze
    print(f"\nFinal time: t={solution['t'][-1]:.2f}")
    print(f"Initial R: {solution['R'][0]:.3f}")
    print(f"Final R: {solution['R'][-1]:.3f}")
    print(f"Energy conservation: {solution['energy'][-1]/solution['energy'][0]:.3f}")

    # Plot results
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Order parameter evolution
    axes[0, 0].plot(solution['t'], solution['R'], 'b-', linewidth=2)
    axes[0, 0].set_xlabel('Time (t)')
    axes[0, 0].set_ylabel('Global R')
    axes[0, 0].set_title('Synchronization Dynamics')
    axes[0, 0].grid(True, alpha=0.3)

    # Energy evolution
    axes[0, 1].plot(solution['t'], solution['energy'], 'r-', linewidth=2)
    axes[0, 1].set_xlabel('Time (t)')
    axes[0, 1].set_ylabel('Energy')
    axes[0, 1].set_title('System Energy')
    axes[0, 1].grid(True, alpha=0.3)

    # Final sync field
    im1 = axes[1, 0].imshow(
        solution['sync_field'][-1].T,
        origin='lower',
        cmap='viridis',
        aspect='auto'
    )
    axes[1, 0].set_title('Final Sync Field R(x,y)')
    axes[1, 0].set_xlabel('x')
    axes[1, 0].set_ylabel('y')
    plt.colorbar(im1, ax=axes[1, 0])

    # Final mediator field
    im2 = axes[1, 1].imshow(
        solution['mediator_field'][-1].T,
        origin='lower',
        cmap='plasma',
        aspect='auto'
    )
    axes[1, 1].set_title('Final Mediator Field σ(x,y)')
    axes[1, 1].set_xlabel('x')
    axes[1, 1].set_ylabel('y')
    plt.colorbar(im2, ax=axes[1, 1])

    plt.tight_layout()
    plt.savefig('smft_basic_evolution.png', dpi=150, bbox_inches='tight')
    print("\nPlot saved: smft_basic_evolution.png")

    return solution


def demo_mass_scaling():
    """Demonstrate M→∞ limit recovery of Kuramoto model."""
    print("\n" + "=" * 60)
    print("SMFT System: Mass Scaling Test")
    print("=" * 60)

    # Test different masses
    M_values = [1.0, 5.0, 10.0, 50.0, 100.0]
    R_final_values = []

    print("\nTesting mediator mass scaling...")
    for M in M_values:
        system = SMFTSystem(
            grid_shape=(30, 30),
            N_oscillators=80,
            coupling='local',
            mediator_mass=M,
            oscillator_frequencies=np.random.normal(0, 1, 80)
        )

        solution = system.evolve((0, 20), dt=0.01, store_interval=100)
        R_final = solution['R'][-1]
        R_final_values.append(R_final)

        print(f"M={M:6.1f}: Final R = {R_final:.4f}")

    # Plot convergence
    plt.figure(figsize=(8, 6))
    plt.semilogx(M_values, R_final_values, 'o-', linewidth=2, markersize=8)
    plt.xlabel('Mediator Mass M')
    plt.ylabel('Steady-state R')
    plt.title('Heavy Mass Limit: M→∞ Recovery')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('smft_mass_scaling.png', dpi=150, bbox_inches='tight')
    print("\nPlot saved: smft_mass_scaling.png")

    return M_values, R_final_values


def demo_local_vs_global():
    """Compare local vs global coupling."""
    print("\n" + "=" * 60)
    print("SMFT System: Local vs Global Coupling")
    print("=" * 60)

    N = 100
    frequencies = np.random.normal(0, 1, N)

    # Local coupling
    print("\nLocal coupling system...")
    system_local = SMFTSystem(
        grid_shape=(40, 40),
        N_oscillators=N,
        coupling='local',
        mediator_mass=10.0,
        oscillator_frequencies=frequencies
    )

    sol_local = system_local.evolve((0, 25), dt=0.01, store_interval=100)

    # Global coupling
    print("Global coupling system...")
    system_global = SMFTSystem(
        grid_shape=(40, 40),
        N_oscillators=N,
        coupling='global',
        mediator_mass=10.0,
        oscillator_frequencies=frequencies
    )

    sol_global = system_global.evolve((0, 25), dt=0.01, store_interval=100)

    # Compare
    print(f"\nLocal final R: {sol_local['R'][-1]:.4f}")
    print(f"Global final R: {sol_global['R'][-1]:.4f}")

    # Spatial structure
    local_variance = np.var(sol_local['sync_field'][-1])
    global_variance = np.var(sol_global['sync_field'][-1])

    print(f"\nLocal field variance: {local_variance:.4f}")
    print(f"Global field variance: {global_variance:.4f}")

    # Plot comparison
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Local R(t)
    axes[0, 0].plot(sol_local['t'], sol_local['R'], 'b-', linewidth=2)
    axes[0, 0].set_title('Local Coupling: R(t)')
    axes[0, 0].set_xlabel('Time (t)')
    axes[0, 0].set_ylabel('R')
    axes[0, 0].grid(True, alpha=0.3)

    # Global R(t)
    axes[0, 1].plot(sol_global['t'], sol_global['R'], 'r-', linewidth=2)
    axes[0, 1].set_title('Global Coupling: R(t)')
    axes[0, 1].set_xlabel('Time (t)')
    axes[0, 1].set_ylabel('R')
    axes[0, 1].grid(True, alpha=0.3)

    # Local field
    im1 = axes[1, 0].imshow(
        sol_local['sync_field'][-1].T,
        origin='lower',
        cmap='viridis',
        aspect='auto'
    )
    axes[1, 0].set_title('Local: R(x,y)')
    plt.colorbar(im1, ax=axes[1, 0])

    # Global field
    im2 = axes[1, 1].imshow(
        sol_global['sync_field'][-1].T,
        origin='lower',
        cmap='viridis',
        aspect='auto'
    )
    axes[1, 1].set_title('Global: R(x,y)')
    plt.colorbar(im2, ax=axes[1, 1])

    plt.tight_layout()
    plt.savefig('smft_local_vs_global.png', dpi=150, bbox_inches='tight')
    print("\nPlot saved: smft_local_vs_global.png")


def demo_effective_mass():
    """Demonstrate effective mass field computation."""
    print("\n" + "=" * 60)
    print("SMFT System: Effective Mass Field")
    print("=" * 60)

    system = SMFTSystem(
        grid_shape=(50, 50),
        N_oscillators=150,
        coupling='local',
        mediator_mass=10.0
    )

    print("\nEvolving to create spatial structure...")
    solution = system.evolve((0, 20), dt=0.01, store_interval=100)

    # Compute effective mass
    print("Computing effective mass field...")
    m_eff = system.compute_effective_mass()

    print(f"Mean effective mass: {np.mean(m_eff):.2f}")
    print(f"Mass range: [{np.min(m_eff):.2f}, {np.max(m_eff):.2f}]")

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    # Sync field
    im1 = axes[0].imshow(
        system.sync_field.values.T,
        origin='lower',
        cmap='viridis',
        aspect='auto'
    )
    axes[0].set_title('Sync Field R(x,y)')
    plt.colorbar(im1, ax=axes[0])

    # Mediator field
    im2 = axes[1].imshow(
        system.mediator_field.values.T,
        origin='lower',
        cmap='plasma',
        aspect='auto'
    )
    axes[1].set_title('Mediator σ(x,y)')
    plt.colorbar(im2, ax=axes[1])

    # Effective mass
    im3 = axes[2].imshow(
        m_eff.T,
        origin='lower',
        cmap='coolwarm',
        aspect='auto'
    )
    axes[2].set_title('Effective Mass m_eff(x,y)')
    plt.colorbar(im3, ax=axes[2])

    plt.tight_layout()
    plt.savefig('smft_effective_mass.png', dpi=150, bbox_inches='tight')
    print("\nPlot saved: smft_effective_mass.png")


def main():
    """Run all demonstrations."""
    print("\n" + "=" * 60)
    print("SMFT SYSTEM DEMONSTRATIONS")
    print("=" * 60)

    # Run demos
    demo_basic_evolution()
    demo_mass_scaling()
    demo_local_vs_global()
    demo_effective_mass()

    print("\n" + "=" * 60)
    print("All demonstrations complete!")
    print("=" * 60)
    print("\nGenerated files:")
    print("  - smft_basic_evolution.png")
    print("  - smft_mass_scaling.png")
    print("  - smft_local_vs_global.png")
    print("  - smft_effective_mass.png")


if __name__ == "__main__":
    main()
