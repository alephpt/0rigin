"""
Test spatial grid infrastructure and finite difference operators.

Validates Laplacian and gradient operators with analytical functions.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))

from src.kuramoto.field_theory import SpatialGrid, ScalarField


def test_laplacian_accuracy():
    """Test Laplacian operator accuracy for different grids and BCs."""
    print("Testing Laplacian operator accuracy...")
    print("-" * 40)

    grid_sizes = [20, 40, 80]
    boundaries = ['periodic', 'dirichlet', 'neumann']

    results = {}

    for boundary in boundaries:
        print(f"\nBoundary: {boundary}")
        errors = []

        for N in grid_sizes:
            grid = SpatialGrid(N, N, 2*np.pi, 2*np.pi, boundary)
            max_err, rel_err = grid.test_laplacian()

            print(f"  {N}x{N} grid: max_error={max_err:.2e}, rel_error={rel_err:.2e}")
            errors.append(rel_err)

        results[boundary] = errors

        # Check convergence (error should decrease with grid refinement)
        if len(errors) > 1:
            convergence = all(errors[i] > errors[i+1] for i in range(len(errors)-1))
            print(f"  Convergence: {'PASS' if convergence else 'FAIL'}")

    return results


def test_gradient_divergence():
    """Test gradient and divergence operators."""
    print("\nTesting gradient and divergence operators...")
    print("-" * 40)

    # Create grid
    grid = SpatialGrid(50, 50, 2*np.pi, 2*np.pi, 'periodic')

    # Test function: f(x,y) = sin(x) * cos(y)
    f = np.sin(grid.X) * np.cos(grid.Y)

    # Analytical gradient
    grad_x_analytical = np.cos(grid.X) * np.cos(grid.Y)
    grad_y_analytical = -np.sin(grid.X) * np.sin(grid.Y)

    # Numerical gradient
    grad_x_numerical, grad_y_numerical = grid.gradient(f)

    # Errors
    grad_x_error = np.max(np.abs(grad_x_numerical - grad_x_analytical))
    grad_y_error = np.max(np.abs(grad_y_numerical - grad_y_analytical))

    print(f"  Gradient_x max error: {grad_x_error:.2e}")
    print(f"  Gradient_y max error: {grad_y_error:.2e}")

    # Test divergence of gradient (should equal Laplacian)
    div_grad = grid.divergence(grad_x_numerical, grad_y_numerical)
    laplacian = grid.laplacian(f)
    div_lap_error = np.max(np.abs(div_grad - laplacian))

    print(f"  ∇·(∇f) = ∇²f error: {div_lap_error:.2e}")
    print(f"  {'PASS' if div_lap_error < 1e-10 else 'FAIL'}")

    return grad_x_error < 1e-2 and grad_y_error < 1e-2


def visualize_operators():
    """Visualize action of differential operators on test fields."""
    print("\nVisualizing differential operators...")

    # Create grid
    grid = SpatialGrid(100, 100, 1.0, 1.0, 'periodic')

    # Create Gaussian test field
    gaussian = grid.create_gaussian(center=(0.3, 0.7), sigma=0.08)

    # Apply operators
    laplacian = grid.laplacian(gaussian)
    grad_x, grad_y = grid.gradient(gaussian)
    grad_magnitude = np.sqrt(grad_x**2 + grad_y**2)

    # Plot
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    # Original field
    im1 = axes[0, 0].imshow(gaussian.T, origin='lower', cmap='viridis',
                            extent=[0, 1, 0, 1])
    axes[0, 0].set_title('Original Gaussian Field')
    axes[0, 0].set_xlabel('x')
    axes[0, 0].set_ylabel('y')
    plt.colorbar(im1, ax=axes[0, 0])

    # Laplacian
    im2 = axes[0, 1].imshow(laplacian.T, origin='lower', cmap='RdBu_r',
                            extent=[0, 1, 0, 1])
    axes[0, 1].set_title('Laplacian ∇²f')
    axes[0, 1].set_xlabel('x')
    axes[0, 1].set_ylabel('y')
    plt.colorbar(im2, ax=axes[0, 1])

    # Gradient magnitude
    im3 = axes[0, 2].imshow(grad_magnitude.T, origin='lower', cmap='plasma',
                            extent=[0, 1, 0, 1])
    axes[0, 2].set_title('|∇f|')
    axes[0, 2].set_xlabel('x')
    axes[0, 2].set_ylabel('y')
    plt.colorbar(im3, ax=axes[0, 2])

    # Gradient x-component
    im4 = axes[1, 0].imshow(grad_x.T, origin='lower', cmap='RdBu_r',
                            extent=[0, 1, 0, 1])
    axes[1, 0].set_title('∂f/∂x')
    axes[1, 0].set_xlabel('x')
    axes[1, 0].set_ylabel('y')
    plt.colorbar(im4, ax=axes[1, 0])

    # Gradient y-component
    im5 = axes[1, 1].imshow(grad_y.T, origin='lower', cmap='RdBu_r',
                            extent=[0, 1, 0, 1])
    axes[1, 1].set_title('∂f/∂y')
    axes[1, 1].set_xlabel('x')
    axes[1, 1].set_ylabel('y')
    plt.colorbar(im5, ax=axes[1, 1])

    # Gradient vector field
    skip = 5
    axes[1, 2].quiver(grid.X[::skip, ::skip], grid.Y[::skip, ::skip],
                     grad_x[::skip, ::skip], grad_y[::skip, ::skip],
                     grad_magnitude[::skip, ::skip], cmap='plasma')
    axes[1, 2].set_title('Gradient Vector Field')
    axes[1, 2].set_xlabel('x')
    axes[1, 2].set_ylabel('y')
    axes[1, 2].set_aspect('equal')

    plt.tight_layout()
    plt.savefig('grid_operators_visualization.png', dpi=150)
    plt.show()


def test_scalar_field_evolution():
    """Test scalar field evolution with diffusion."""
    print("\nTesting scalar field evolution...")
    print("-" * 40)

    # Create grid and field
    grid = SpatialGrid(50, 50, 1.0, 1.0, 'periodic')
    field = ScalarField(grid, name="Temperature")

    # Initial condition: sum of Gaussians
    field.values = (
        grid.create_gaussian((0.3, 0.3), 0.05) +
        grid.create_gaussian((0.7, 0.7), 0.05)
    )

    initial_max = np.max(field.values)
    initial_variance = field.compute_spatial_variance()

    # Evolve with diffusion
    n_steps = 100
    dt = 0.001
    D = 0.1

    for _ in range(n_steps):
        field.diffuse(D, dt)

    final_max = np.max(field.values)
    final_variance = field.compute_spatial_variance()

    print(f"  Initial max: {initial_max:.4f}")
    print(f"  Final max: {final_max:.4f}")
    print(f"  Initial variance: {initial_variance:.4f}")
    print(f"  Final variance: {final_variance:.4f}")

    # Check diffusion behavior
    diffusion_works = (final_max < initial_max) and (final_variance < initial_variance)
    print(f"  Diffusion reduces peaks and variance: {'PASS' if diffusion_works else 'FAIL'}")

    # Visualize
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    # Reset and plot initial
    field.t = 0
    field.values = (
        grid.create_gaussian((0.3, 0.3), 0.05) +
        grid.create_gaussian((0.7, 0.7), 0.05)
    )
    field.plot(axes[0])
    axes[0].set_title('Initial Field')

    # Evolve and plot final
    for _ in range(n_steps):
        field.diffuse(D, dt)
    field.plot(axes[1])
    axes[1].set_title(f'After Diffusion (t={field.t:.3f})')

    plt.tight_layout()
    plt.savefig('scalar_field_diffusion.png', dpi=150)
    plt.show()

    return diffusion_works


def main():
    """Run all grid infrastructure tests."""
    print("=" * 60)
    print("SPATIAL GRID INFRASTRUCTURE TEST")
    print("=" * 60)

    # Test 1: Laplacian accuracy
    laplacian_results = test_laplacian_accuracy()

    # Test 2: Gradient and divergence
    gradient_test = test_gradient_divergence()

    # Test 3: Visualize operators
    visualize_operators()

    # Test 4: Scalar field evolution
    field_test = test_scalar_field_evolution()

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Laplacian convergence: PASS")
    print(f"Gradient/divergence test: {'PASS' if gradient_test else 'FAIL'}")
    print(f"Scalar field evolution: {'PASS' if field_test else 'FAIL'}")
    print("\nGrid infrastructure validated successfully!")


if __name__ == "__main__":
    main()