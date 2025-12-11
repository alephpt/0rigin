"""
Simple test script to verify implementation without pytest.
"""

import sys
sys.path.insert(0, '/home/persist/neotec/0rigin/src')

import numpy as np

def test_distributions():
    """Test all frequency distributions."""
    print("Testing Distributions...")

    from kuramoto.distributions import (
        LorentzianDistribution,
        GaussianDistribution,
        UniformDistribution
    )

    # Test Lorentzian
    dist = LorentzianDistribution(center=0, width=1.0)
    freqs = dist.sample(100, seed=42)
    assert len(freqs) == 100
    assert dist.critical_coupling() == 2.0
    pdf_vals = dist.pdf(np.array([0.0, 1.0]))
    assert pdf_vals.shape == (2,)
    print("  ✓ LorentzianDistribution")

    # Test Gaussian
    dist = GaussianDistribution(mean=0, std=1.0)
    freqs = dist.sample(100, seed=42)
    assert len(freqs) == 100
    Kc = dist.critical_coupling()
    assert abs(Kc - 1.596) < 0.01
    pdf_vals = dist.pdf(np.array([0.0, 1.0]))
    assert pdf_vals.shape == (2,)
    print("  ✓ GaussianDistribution")

    # Test Uniform
    dist = UniformDistribution(low=-1, high=1)
    freqs = dist.sample(100, seed=42)
    assert len(freqs) == 100
    assert np.all((freqs >= -1) & (freqs <= 1))
    assert dist.critical_coupling() is None
    pdf_vals = dist.pdf(np.array([-0.5, 0.5, 2.0]))
    assert pdf_vals[0] > 0 and pdf_vals[1] > 0 and pdf_vals[2] == 0
    print("  ✓ UniformDistribution")


def test_kuramoto_model():
    """Test KuramotoModel with all distributions."""
    print("\nTesting KuramotoModel...")

    from kuramoto import KuramotoModel
    from kuramoto.distributions import LorentzianDistribution

    # Test with string distribution
    model = KuramotoModel(N=50, coupling=2.0, frequencies='lorentzian')
    assert model.N == 50
    print("  ✓ Model initialization (string)")

    # Test with distribution object
    dist = LorentzianDistribution(center=0, width=1.0)
    model = KuramotoModel(N=50, coupling=2.0, frequencies=dist)
    assert model.N == 50
    print("  ✓ Model initialization (distribution)")

    # Test with array
    freqs = np.random.normal(0, 1, 50)
    model = KuramotoModel(N=50, coupling=2.0, frequencies=freqs)
    assert model.N == 50
    print("  ✓ Model initialization (array)")

    # Test all distribution types
    for dist_name in ['lorentzian', 'gaussian', 'uniform']:
        model = KuramotoModel(N=30, coupling=2.0, frequencies=dist_name)
        solution = model.evolve((0, 10), solver='rk45')
        assert 't' in solution
        assert 'phases' in solution
        assert 'R' in solution
        assert 'Psi' in solution
        assert len(solution['t']) > 0
        print(f"  ✓ Simulation with {dist_name}")


def test_order_parameter():
    """Test order parameter calculation."""
    print("\nTesting Order Parameter...")

    from kuramoto import KuramotoModel
    from kuramoto.analysis import OrderParameter

    # Run simulation
    model = KuramotoModel(N=50, coupling=3.0, frequencies='lorentzian')
    solution = model.evolve((0, 20), solver='rk45')

    # Test single snapshot
    phases = solution['phases'][0, :]
    R, Psi = model.compute_order_parameter(phases)
    assert 0 <= R <= 1
    assert -np.pi <= Psi <= np.pi
    print("  ✓ Single snapshot order parameter")

    # Test OrderParameter class
    op = OrderParameter(solution['phases'])
    R_series, Psi_series = op.time_series()
    assert len(R_series) == len(solution['t'])
    assert np.all((R_series >= 0) & (R_series <= 1))
    print("  ✓ OrderParameter time series")

    # Test analysis methods
    mean_R = op.mean_amplitude()
    assert 0 <= mean_R <= 1
    print("  ✓ Mean amplitude")

    steady_R = op.steady_state_amplitude(fraction=0.2)
    assert 0 <= steady_R <= 1
    print("  ✓ Steady-state amplitude")

    is_sync = op.is_synchronized(threshold=0.3)
    assert isinstance(is_sync, (bool, np.bool_))
    print("  ✓ Synchronization check")


def test_solvers():
    """Test different solvers."""
    print("\nTesting Solvers...")

    from kuramoto import KuramotoModel

    model = KuramotoModel(N=30, coupling=2.0, frequencies='lorentzian')

    for solver_name in ['rk45', 'rk4', 'euler']:
        solution = model.evolve((0, 10), solver=solver_name, dt=0.1)
        assert len(solution['t']) > 0
        assert solution['phases'].shape[0] == len(solution['t'])
        assert solution['phases'].shape[1] == 30
        print(f"  ✓ {solver_name.upper()} solver")


def test_visualization():
    """Test visualization imports."""
    print("\nTesting Visualization...")

    try:
        from kuramoto.visualization import (
            plot_phases,
            plot_order_parameter,
            plot_bifurcation_diagram
        )
        print("  ✓ Visualization imports")
    except ImportError as e:
        print(f"  ⚠ Visualization requires matplotlib: {e}")


def test_metrics():
    """Test analysis metrics."""
    print("\nTesting Analysis Metrics...")

    from kuramoto.analysis import (
        phase_coherence,
        phase_variance,
        metastability,
        frequency_entrainment
    )

    # Test phase coherence
    phases = np.random.uniform(0, 2*np.pi, 100)
    rho = phase_coherence(phases)
    assert 0 <= rho <= 1
    print("  ✓ Phase coherence")

    # Test phase variance
    var = phase_variance(phases, circular=True)
    assert 0 <= var <= 1
    print("  ✓ Phase variance")

    # Test metastability
    R_series = np.random.uniform(0.3, 0.7, 100)
    M = metastability(R_series)
    assert M >= 0
    print("  ✓ Metastability")

    # Test frequency entrainment
    freqs = np.random.normal(0, 0.5, 100)
    ent = frequency_entrainment(freqs, threshold=0.1)
    assert 0 <= ent <= 1
    print("  ✓ Frequency entrainment")


def test_synchronization_regimes():
    """Test synchronization in different regimes."""
    print("\nTesting Synchronization Regimes...")

    from kuramoto import KuramotoModel
    from kuramoto.distributions import LorentzianDistribution

    N = 100
    gamma = 1.0
    Kc = 2 * gamma

    dist = LorentzianDistribution(center=0, width=gamma)

    # Subcritical: should have low R
    model = KuramotoModel(N=N, coupling=0.5*Kc, frequencies=dist)
    solution = model.evolve((0, 50), solver='rk45')
    R_steady = np.mean(solution['R'][-100:])
    assert R_steady < 0.3, f"Subcritical R too high: {R_steady}"
    print(f"  ✓ Subcritical regime (R = {R_steady:.3f})")

    # Supercritical: should have high R
    model = KuramotoModel(N=N, coupling=2*Kc, frequencies=dist)
    solution = model.evolve((0, 50), solver='rk45')
    R_steady = np.mean(solution['R'][-100:])
    R_theory = dist.steady_state_order_parameter(2*Kc)
    assert R_steady > 0.5, f"Supercritical R too low: {R_steady}"
    # Allow larger tolerance due to finite-size effects
    assert abs(R_steady - R_theory) < 0.2, f"R mismatch: {R_steady} vs {R_theory}"
    print(f"  ✓ Supercritical regime (R = {R_steady:.3f}, theory = {R_theory:.3f})")


def main():
    """Run all tests."""
    print("=" * 60)
    print("Kuramoto Model Implementation Tests")
    print("=" * 60)

    try:
        test_distributions()
        test_kuramoto_model()
        test_order_parameter()
        test_solvers()
        test_metrics()
        test_synchronization_regimes()
        test_visualization()

        print("\n" + "=" * 60)
        print("All tests passed! ✓")
        print("=" * 60)
        return 0

    except Exception as e:
        print(f"\n✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    exit(main())
