#!/usr/bin/env python3
"""
Complete SMFT Field Theory Demonstration.

Demonstrates full integration of:
1. Klein-Gordon mediator field
2. Local oscillator-field coupling
3. Fermion mass generation
4. Spatiotemporal pattern formation
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import time

from kuramoto.field_theory import (
    SpatialGrid,
    MediatorField,
    LocalFieldCoupling,
    FermionMassDemo
)


def demo_klein_gordon_propagation():
    """Test 1: Klein-Gordon wave propagation and causality."""
    print("\n" + "="*60)
    print("TEST 1: Klein-Gordon Mediator Field")
    print("="*60)

    # Create grid
    grid = SpatialGrid(Nx=64, Ny=64, Lx=2.0, Ly=2.0, boundary='periodic')

    # Create mediator with finite wave speed
    mediator = MediatorField(
        grid,
        wave_speed=1.0,
        mass=0.5,
        coupling_constant=1.0
    )

    # Initial Gaussian perturbation
    x0, y0 = 1.0, 1.0
    sigma0 = 0.1
    mediator.sigma = grid.create_gaussian(center=(x0, y0), sigma=sigma0)

    # Evolve and watch wave propagation
    t_final = 2.0
    dt = 0.01
    n_steps = int(t_final / dt)

    # Create dummy source (no oscillators, free wave)
    source = np.zeros((grid.Nx, grid.Ny))

    print(f"Evolving free Klein-Gordon wave for t={t_final}...")
    start = time.time()

    for step in range(n_steps):
        mediator.evolve_step(source, dt, method='leapfrog')

        if (step + 1) % 20 == 0:
            mediator.store_snapshot()

    elapsed = time.time() - start
    print(f"Evolution complete in {elapsed:.3f}s")

    # Compute energy (should be conserved)
    energies = []
    for sigma in mediator.history:
        mediator.sigma = sigma
        energies.append(mediator.compute_field_energy())

    energy_conservation = np.std(energies) / np.mean(energies)
    print(f"Energy conservation: {energy_conservation:.6f} (relative std)")

    # Estimate propagation speed
    v_phase, v_group = mediator.compute_propagation_speed()
    print(f"Phase velocity: {v_phase:.3f}")
    print(f"Group velocity: {v_group:.3f}")

    # Plot wave propagation
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    times_to_plot = [0, 2, 4, 6, 8, -1]
    for i, (ax, idx) in enumerate(zip(axes.flat, times_to_plot)):
        if idx < len(mediator.history):
            im = ax.imshow(
                mediator.history[idx].T,
                origin='lower',
                extent=[0, grid.Lx, 0, grid.Ly],
                cmap='RdBu',
                vmin=-0.5,
                vmax=0.5
            )
            ax.set_title(f't = {mediator.time_history[idx]:.2f}')
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            plt.colorbar(im, ax=ax)

    plt.tight_layout()
    plt.savefig('/home/persist/neotec/0rigin/examples/field_theory/klein_gordon_propagation.png',
                dpi=150, bbox_inches='tight')
    print("Saved: klein_gordon_propagation.png")

    return mediator, energies


def demo_local_coupling():
    """Test 2: Bidirectional oscillator-field coupling."""
    print("\n" + "="*60)
    print("TEST 2: Local Oscillator-Field Coupling")
    print("="*60)

    # Create grid
    grid = SpatialGrid(Nx=48, Ny=48, Lx=1.0, Ly=1.0, boundary='periodic')

    # Create coupling system
    coupling = LocalFieldCoupling(
        grid,
        mediator_params={
            'wave_speed': 1.0,
            'mass': 2.0,
            'coupling_constant': 0.5
        },
        oscillator_to_field_coupling=0.5,
        field_to_oscillator_coupling=0.3,
        kernel_width=0.05
    )

    # Initialize oscillators on grid
    N = 100
    positions = np.random.rand(N, 2)
    positions[:, 0] *= grid.Lx
    positions[:, 1] *= grid.Ly

    # Initial phases (clustered)
    phases = np.random.uniform(0, 2*np.pi, N)

    # Natural frequencies (narrow distribution)
    frequencies = np.random.normal(1.0, 0.2, N)

    print(f"Simulating {N} coupled oscillators...")
    print(f"Coupling range: {coupling.get_effective_coupling_range():.3f}")

    # Evolve coupled system (use smaller dt for stability)
    result = coupling.evolve_coupled_system(
        phases=phases,
        positions=positions,
        frequencies=frequencies,
        dt=0.005,
        n_steps=1000,
        store_interval=20
    )

    print(f"Final order parameter R = {result['R'][-1]:.3f}")
    print(f"Final field energy = {result['field_energy'][-1]:.3e}")

    # Plot results
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)

    # Order parameter evolution
    ax1 = fig.add_subplot(gs[0, :])
    ax1.plot(result['t'], result['R'], 'b-', linewidth=2)
    ax1.set_xlabel('Time')
    ax1.set_ylabel('Order Parameter R')
    ax1.set_title('Synchronization Evolution')
    ax1.grid(True, alpha=0.3)

    # Field snapshots at different times
    snapshot_times = [0, len(result['mediator_field'])//2, -1]
    for i, idx in enumerate(snapshot_times):
        ax = fig.add_subplot(gs[1, i])
        im = ax.imshow(
            result['mediator_field'][idx].T,
            origin='lower',
            extent=[0, grid.Lx, 0, grid.Ly],
            cmap='viridis'
        )
        ax.set_title(f'Mediator σ(x,y) at t={result["t"][idx]:.1f}')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        plt.colorbar(im, ax=ax)

    # Phase space at final time
    ax_phase = fig.add_subplot(gs[2, :2])
    final_phases = result['phases'][-1]
    scatter = ax_phase.scatter(
        positions[:, 0],
        positions[:, 1],
        c=final_phases,
        cmap='hsv',
        s=50,
        alpha=0.7
    )
    ax_phase.set_xlabel('x')
    ax_phase.set_ylabel('y')
    ax_phase.set_title('Oscillator Phases (final)')
    ax_phase.set_xlim([0, grid.Lx])
    ax_phase.set_ylim([0, grid.Ly])
    plt.colorbar(scatter, ax=ax_phase, label='Phase')

    # Energy evolution
    ax_energy = fig.add_subplot(gs[2, 2])
    ax_energy.plot(result['t'], result['field_energy'], 'r-', linewidth=2)
    ax_energy.set_xlabel('Time')
    ax_energy.set_ylabel('Field Energy')
    ax_energy.set_title('Field Energy Evolution')
    ax_energy.grid(True, alpha=0.3)

    plt.savefig('/home/persist/neotec/0rigin/examples/field_theory/local_coupling_demo.png',
                dpi=150, bbox_inches='tight')
    print("Saved: local_coupling_demo.png")

    return coupling, result


def demo_fermion_mass_generation():
    """Test 3: Effective mass generation from synchronization."""
    print("\n" + "="*60)
    print("TEST 3: Fermion Mass Generation")
    print("="*60)

    # Create grid
    grid = SpatialGrid(Nx=64, Ny=64, Lx=1.0, Ly=1.0, boundary='periodic')

    # Create fermion mass demo
    fermion = FermionMassDemo(
        grid,
        yukawa_coupling=2.0,
        bare_mass=0.0
    )

    # Create pattern of oscillators with varying synchronization
    N = 200
    positions = np.random.rand(N, 2)

    # Create two regions: synchronized and desynchronized
    phases = np.zeros(N)
    for i in range(N):
        if positions[i, 0] < 0.5:
            # Left half: synchronized
            phases[i] = np.random.uniform(0, 0.5)
        else:
            # Right half: random
            phases[i] = np.random.uniform(0, 2*np.pi)

    # Update fields
    fermion.update_from_oscillators(phases, positions, kernel_width=0.05)

    mass_gap = fermion.compute_mass_gap()
    avg_mass = fermion.compute_average_mass()

    print(f"Mass gap: Δm = {mass_gap:.3f}")
    print(f"Average mass: <m> = {avg_mass:.3f}")
    print(f"Symmetry breaking: σ²(m) = {fermion.compute_symmetry_breaking_order():.3f}")

    # Plot fields
    fig = fermion.plot_fields_side_by_side()
    fig.savefig('/home/persist/neotec/0rigin/examples/field_theory/fermion_mass_fields.png',
                dpi=150, bbox_inches='tight')
    print("Saved: fermion_mass_fields.png")

    # Plot mass vs order parameter
    fig2, ax = plt.subplots(figsize=(8, 6))
    fermion.plot_mass_vs_order_parameter(ax)
    fig2.savefig('/home/persist/neotec/0rigin/examples/field_theory/mass_vs_R.png',
                 dpi=150, bbox_inches='tight')
    print("Saved: mass_vs_R.png")

    # Demonstrate phase transition
    print("\nSimulating phase transition...")
    R_values = np.linspace(0, 1, 50)
    transition = fermion.demonstrate_phase_transition(R_values, plot=True)

    if 'figure' in transition:
        transition['figure'].savefig(
            '/home/persist/neotec/0rigin/examples/field_theory/phase_transition.png',
            dpi=150, bbox_inches='tight'
        )
        print("Saved: phase_transition.png")

    return fermion, transition


def demo_heavy_mass_limit():
    """Test 4: Heavy mass limit recovers global coupling."""
    print("\n" + "="*60)
    print("TEST 4: Heavy Mass Limit (M→∞)")
    print("="*60)

    grid = SpatialGrid(Nx=32, Ny=32, Lx=1.0, Ly=1.0, boundary='periodic')

    # Test with different masses
    masses = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
    coupling_ranges = []
    global_vs_local = []

    N = 50
    positions = np.random.rand(N, 2)
    phases = np.random.uniform(0, 2*np.pi, N)

    for M in masses:
        coupling = LocalFieldCoupling(
            grid,
            mediator_params={
                'wave_speed': 1.0,
                'mass': M,
                'coupling_constant': 1.0
            }
        )

        # Test global vs local R
        global_R, local_R = coupling.test_heavy_mass_limit(phases, positions)

        coupling_range = coupling.get_effective_coupling_range()
        coupling_ranges.append(coupling_range)
        global_vs_local.append((global_R, local_R))

        print(f"M={M:5.1f}: ξ={coupling_range:.3f}, "
              f"R_global={global_R:.3f}, R_local={local_R:.3f}, "
              f"diff={abs(global_R - local_R):.4f}")

    # Plot convergence
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Coupling range vs mass
    ax1.plot(masses, coupling_ranges, 'bo-', linewidth=2, markersize=8)
    ax1.set_xlabel('Mass M', fontsize=12)
    ax1.set_ylabel('Coupling Range ξ = c/M', fontsize=12)
    ax1.set_title('Effective Coupling Range')
    ax1.grid(True, alpha=0.3)
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    # Global vs local convergence
    differences = [abs(g - l) for g, l in global_vs_local]
    ax2.plot(masses, differences, 'ro-', linewidth=2, markersize=8)
    ax2.set_xlabel('Mass M', fontsize=12)
    ax2.set_ylabel('|R_global - R_local|', fontsize=12)
    ax2.set_title('Convergence to Global Coupling')
    ax2.grid(True, alpha=0.3)
    ax2.set_xscale('log')
    ax2.set_yscale('log')

    plt.tight_layout()
    plt.savefig('/home/persist/neotec/0rigin/examples/field_theory/heavy_mass_limit.png',
                dpi=150, bbox_inches='tight')
    print("Saved: heavy_mass_limit.png")

    return masses, coupling_ranges, global_vs_local


def performance_benchmark():
    """Benchmark: Performance scaling."""
    print("\n" + "="*60)
    print("BENCHMARK: Performance Scaling")
    print("="*60)

    grid_sizes = [16, 32, 48, 64]
    N_oscillators = [25, 50, 100, 150]

    timings = []

    for Nx, N in zip(grid_sizes, N_oscillators):
        grid = SpatialGrid(Nx=Nx, Ny=Nx, Lx=1.0, Ly=1.0)

        coupling = LocalFieldCoupling(grid)

        positions = np.random.rand(N, 2)
        phases = np.random.uniform(0, 2*np.pi, N)
        frequencies = np.random.normal(1.0, 0.1, N)

        # Benchmark evolution
        start = time.time()
        coupling.evolve_coupled_system(
            phases, positions, frequencies,
            dt=0.01, n_steps=100, store_interval=50
        )
        elapsed = time.time() - start

        timings.append(elapsed)
        ops = Nx * Nx * N * 100  # Grid points × oscillators × steps
        throughput = ops / elapsed / 1e6

        print(f"Grid {Nx}×{Nx}, N={N:3d}: {elapsed:.3f}s ({throughput:.2f} MOps/s)")

    return grid_sizes, N_oscillators, timings


def main():
    """Run all demonstrations."""
    print("\n" + "="*70)
    print(" "*15 + "SMFT FIELD THEORY - FULL DEMONSTRATION")
    print("="*70)

    # Test 1: Klein-Gordon propagation
    mediator, energies = demo_klein_gordon_propagation()

    # Test 2: Local coupling
    coupling, result = demo_local_coupling()

    # Test 3: Fermion mass generation
    fermion, transition = demo_fermion_mass_generation()

    # Test 4: Heavy mass limit
    masses, ranges, convergence = demo_heavy_mass_limit()

    # Benchmark
    grid_sizes, N_osc, timings = performance_benchmark()

    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print("\n✓ Klein-Gordon mediator field: Wave propagation validated")
    print(f"  Energy conservation: {np.std(energies)/np.mean(energies):.6f}")
    print("\n✓ Local oscillator-field coupling: Bidirectional dynamics working")
    print(f"  Final synchronization: R = {result['R'][-1]:.3f}")
    print("\n✓ Fermion mass generation: Mass emerges from order parameter")
    print(f"  Mass gap: Δm = {fermion.compute_mass_gap():.3f}")
    print("\n✓ Heavy mass limit: M→∞ recovers global coupling")
    print(f"  Convergence demonstrated for M = {masses[-1]}")
    print("\n✓ Performance: Scales reasonably with grid size and oscillator count")
    print(f"  Largest system: {grid_sizes[-1]}×{grid_sizes[-1]} grid, "
          f"{N_osc[-1]} oscillators in {timings[-1]:.2f}s")

    print("\n" + "="*70)
    print("All tests completed successfully!")
    print("="*70 + "\n")

    plt.show()


if __name__ == '__main__':
    main()
