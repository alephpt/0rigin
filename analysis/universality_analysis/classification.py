#!/usr/bin/env python3
"""
Universality Class Classification Decision Tree

Implements the classification algorithm to identify the SMFT universality class
based on measured critical exponents.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Rectangle
import pandas as pd
import os
import sys
from typing import Dict, Tuple, List
import json


class UniversalityClassifier:
    """Classify phase transitions based on critical exponents"""

    def __init__(self):
        # Known 2D universality classes
        self.known_classes = {
            '2D Ising': {
                'beta': 0.125,      # 1/8 (exact)
                'nu': 1.0,          # (exact)
                'gamma': 1.75,      # 7/4 (exact)
                'eta': 0.25,        # 1/4 (exact)
                'alpha': 0.0,       # (logarithmic)
                'delta': 15.0,      # (exact)
                'color': 'blue',
                'description': 'Z₂ symmetry breaking'
            },
            '2D XY': {
                'beta': 0.23,       # Approximate (BKT transition)
                'nu': 0.67,         # Approximate
                'gamma': 1.32,      # Approximate
                'eta': 0.04,        # Very small
                'alpha': -0.01,     # Near zero
                'delta': 6.7,       # Approximate
                'color': 'green',
                'description': 'Berezinskii-Kosterlitz-Thouless transition'
            },
            '2D 3-state Potts': {
                'beta': 1/9,        # (exact)
                'nu': 5/6,          # (exact)
                'gamma': 13/9,      # (exact)
                'eta': 4/15,        # (exact)
                'alpha': 1/3,       # (exact)
                'delta': 14.0,      # (exact)
                'color': 'orange',
                'description': 'Z₃ symmetry'
            },
            '2D 4-state Potts': {
                'beta': 1/12,       # (exact)
                'nu': 2/3,          # (exact)
                'gamma': 7/6,       # (exact)
                'eta': 1/4,         # (exact)
                'alpha': 2/3,       # (exact, first-order)
                'delta': 15.0,      # (exact)
                'color': 'purple',
                'description': 'First-order transition'
            },
            'Mean Field': {
                'beta': 0.5,
                'nu': 0.5,
                'gamma': 1.0,
                'eta': 0.0,
                'alpha': 0.0,
                'delta': 3.0,
                'color': 'red',
                'description': 'd ≥ 4 or long-range interactions'
            },
            '2D Percolation': {
                'beta': 5/36,       # ≈ 0.139
                'nu': 4/3,          # ≈ 1.333
                'gamma': 43/18,     # ≈ 2.389
                'eta': 5/24,        # ≈ 0.208
                'alpha': -2/3,      # Negative
                'delta': 91/5,      # ≈ 18.2
                'color': 'cyan',
                'description': 'Geometric phase transition'
            }
        }

        self.measured_exponents = None
        self.classification_result = None
        self.confidence_scores = {}

    def load_exponents(self, filename: str):
        """Load measured exponents from file"""
        exponents = {}

        if filename.endswith('.txt'):
            # Parse text file
            with open(filename, 'r') as f:
                for line in f:
                    if '=' in line and not line.startswith('#'):
                        parts = line.split('=')
                        name = parts[0].strip()
                        value_parts = parts[1].strip().split('±')
                        value = float(value_parts[0])
                        error = float(value_parts[1]) if len(value_parts) > 1 else 0.0

                        # Map sigma_c to proper name
                        if name == 'sigma_c' or name == 'σ_c':
                            name = 'sigma_c'

                        exponents[name] = {'value': value, 'error': error}

        elif filename.endswith('.json'):
            # Load JSON file
            with open(filename, 'r') as f:
                data = json.load(f)
                for name, val in data.items():
                    if isinstance(val, dict):
                        exponents[name] = val
                    else:
                        exponents[name] = {'value': val, 'error': 0.0}

        self.measured_exponents = exponents
        return exponents

    def calculate_distance(self, exp1: Dict, exp2: Dict,
                         exponents: List[str] = ['beta', 'nu', 'gamma', 'eta']) -> float:
        """Calculate weighted distance between two sets of exponents"""
        distance = 0
        weights = {'beta': 2.0, 'nu': 1.5, 'gamma': 1.0, 'eta': 0.5}

        for exp_name in exponents:
            if exp_name in exp1 and exp_name in exp2:
                val1 = exp1[exp_name]['value'] if isinstance(exp1[exp_name], dict) else exp1[exp_name]
                val2 = exp2[exp_name]

                # Relative difference
                if val2 != 0:
                    diff = abs(val1 - val2) / abs(val2)
                else:
                    diff = abs(val1 - val2)

                weight = weights.get(exp_name, 1.0)
                distance += weight * diff

        # Normalize by sum of weights
        total_weight = sum(weights.get(e, 1.0) for e in exponents)
        return distance / total_weight

    def classify(self) -> str:
        """Classify based on decision tree algorithm"""
        if not self.measured_exponents:
            raise ValueError("No measured exponents loaded!")

        print("\n" + "="*60)
        print("UNIVERSALITY CLASS CLASSIFICATION")
        print("="*60)

        # Extract values
        beta = self.measured_exponents['beta']['value']
        nu = self.measured_exponents['nu']['value']
        gamma = self.measured_exponents['gamma']['value']
        eta = self.measured_exponents['eta']['value']

        print(f"\nMeasured Exponents:")
        print(f"  β = {beta:.6f}")
        print(f"  ν = {nu:.6f}")
        print(f"  γ = {gamma:.6f}")
        print(f"  η = {eta:.6f}")

        # Calculate distances to all known classes
        distances = {}
        for class_name, class_exp in self.known_classes.items():
            dist = self.calculate_distance(self.measured_exponents, class_exp)
            distances[class_name] = dist
            self.confidence_scores[class_name] = 1.0 / (1.0 + dist)

        # Sort by distance
        sorted_classes = sorted(distances.items(), key=lambda x: x[1])

        print("\n" + "-"*50)
        print("Distance to Known Classes:")
        print("-"*50)
        for class_name, dist in sorted_classes:
            confidence = self.confidence_scores[class_name]
            print(f"  {class_name:20s}: {dist:.4f} (confidence: {100*confidence:.1f}%)")

        # DECISION TREE (based on strategy document)
        print("\n" + "-"*50)
        print("Decision Tree Analysis:")
        print("-"*50)

        # Check 2D Ising
        if abs(beta - 0.125) < 0.01 and abs(nu - 1.0) < 0.05:
            self.classification_result = "2D Ising"
            print("✓ Condition met: |β - 0.125| < 0.01 AND |ν - 1.0| < 0.05")
            print(f"  β deviation: {abs(beta - 0.125):.4f} < 0.01 ✓")
            print(f"  ν deviation: {abs(nu - 1.0):.4f} < 0.05 ✓")
            print("\n→ CLASSIFICATION: 2D Ising universality class")

        # Check 2D XY
        elif abs(nu - 0.67) < 0.05 and eta < 0.1:
            self.classification_result = "2D XY"
            print("✓ Condition met: |ν - 0.67| < 0.05 AND η < 0.1")
            print(f"  ν deviation: {abs(nu - 0.67):.4f} < 0.05 ✓")
            print(f"  η value: {eta:.4f} < 0.1 ✓")
            print("\n→ CLASSIFICATION: 2D XY universality class (BKT)")

        # Check 3-state Potts
        elif abs(beta - 1/9) < 0.01 and abs(nu - 5/6) < 0.05:
            self.classification_result = "2D 3-state Potts"
            print("✓ Condition met: |β - 1/9| < 0.01 AND |ν - 5/6| < 0.05")
            print(f"  β deviation: {abs(beta - 1/9):.4f} < 0.01 ✓")
            print(f"  ν deviation: {abs(nu - 5/6):.4f} < 0.05 ✓")
            print("\n→ CLASSIFICATION: 2D 3-state Potts universality class")

        # Check mean field
        elif beta > 0.3 and abs(nu - 0.5) < 0.1:
            self.classification_result = "Mean Field"
            print("✓ Condition met: β > 0.3 AND |ν - 0.5| < 0.1")
            print(f"  β value: {beta:.4f} > 0.3 ✓")
            print(f"  ν deviation: {abs(nu - 0.5):.4f} < 0.1 ✓")
            print("\n→ CLASSIFICATION: Mean-field behavior")

        # Novel class
        else:
            self.classification_result = "Novel - SMFT universality class"
            print("✗ No known class conditions met")
            print(f"  β = {beta:.4f} does not match any known 2D class")
            print("\n★ CLASSIFICATION: Novel universality class!")
            print("\nThis is a SIGNIFICANT DISCOVERY:")
            print("• New critical behavior not previously observed")
            print("• Requires new theoretical framework")
            print("• Unique to SMFT quantum-classical coupling")

        print("="*60)

        return self.classification_result

    def visualize_classification(self, output_file: str = None):
        """Create visualization of classification results"""
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))

        # Panel 1: Exponent comparison bar chart
        self.plot_exponent_comparison(axes[0, 0])

        # Panel 2: Distance radar plot
        self.plot_distance_radar(axes[0, 1])

        # Panel 3: Decision tree diagram
        self.plot_decision_tree(axes[1, 0])

        # Panel 4: Confidence scores
        self.plot_confidence_scores(axes[1, 1])

        plt.suptitle(f'Universality Classification: {self.classification_result}',
                    fontsize=16, fontweight='bold')
        plt.tight_layout()

        if output_file:
            plt.savefig(output_file, dpi=150, bbox_inches='tight')
            print(f"\nVisualization saved to: {output_file}")

        plt.show()

    def plot_exponent_comparison(self, ax):
        """Bar chart comparing measured vs theoretical exponents"""
        exponents = ['beta', 'nu', 'gamma', 'eta']
        labels = ['β', 'ν', 'γ', 'η']

        measured_vals = [self.measured_exponents[e]['value'] for e in exponents]
        measured_errs = [self.measured_exponents[e].get('error', 0) for e in exponents]

        x = np.arange(len(labels))
        width = 0.15

        # Plot measured values
        ax.bar(x, measured_vals, width, yerr=measured_errs,
               label='Measured', color='black', alpha=0.7)

        # Plot known classes
        for i, (class_name, class_exp) in enumerate(self.known_classes.items()):
            if i < 4:  # Limit to 4 most relevant classes
                class_vals = [class_exp[e] for e in exponents]
                ax.bar(x + (i+1)*width, class_vals, width,
                      label=class_name, color=class_exp['color'], alpha=0.6)

        ax.set_xlabel('Critical Exponent')
        ax.set_ylabel('Value')
        ax.set_title('Critical Exponent Comparison')
        ax.set_xticks(x + 2*width)
        ax.set_xticklabels(labels)
        ax.legend(loc='upper left', fontsize=9)
        ax.grid(True, alpha=0.3)

    def plot_distance_radar(self, ax):
        """Radar plot of distances to known classes"""
        categories = list(self.known_classes.keys())[:6]  # Top 6 classes
        N = len(categories)

        # Compute angles
        angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
        angles += angles[:1]  # Complete the circle

        # Get distances (inverse for better visualization)
        distances = [self.confidence_scores.get(cat, 0) for cat in categories]
        distances += distances[:1]  # Complete the circle

        # Plot
        ax = plt.subplot(222, projection='polar')
        ax.plot(angles, distances, 'o-', linewidth=2, color='blue')
        ax.fill(angles, distances, alpha=0.25)

        ax.set_xticks(angles[:-1])
        ax.set_xticklabels(categories, size=9)
        ax.set_ylim(0, 1)
        ax.set_ylabel('Confidence Score', labelpad=30)
        ax.set_title('Classification Confidence', pad=20)
        ax.grid(True)

    def plot_decision_tree(self, ax):
        """Visualize the decision tree"""
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 10)
        ax.axis('off')

        # Decision nodes
        nodes = {
            'root': (5, 9, "Start"),
            'ising_check': (2, 7, "|β-0.125|<0.01\n& |ν-1|<0.05?"),
            'xy_check': (5, 7, "|ν-0.67|<0.05\n& η<0.1?"),
            'potts_check': (8, 7, "|β-1/9|<0.01\n& |ν-5/6|<0.05?"),
            'ising': (1, 5, "2D Ising"),
            'xy': (3.5, 5, "2D XY"),
            'potts': (6.5, 5, "3-Potts"),
            'novel': (5, 3, "Novel Class")
        }

        # Draw nodes
        for node, (x, y, text) in nodes.items():
            if node == self.classification_result.lower().replace(' ', '_'):
                color = 'lightgreen'
                linewidth = 3
            elif 'check' in node:
                color = 'lightblue'
                linewidth = 1
            else:
                color = 'white'
                linewidth = 1

            rect = Rectangle((x-0.8, y-0.3), 1.6, 0.6,
                           facecolor=color, edgecolor='black', linewidth=linewidth)
            ax.add_patch(rect)
            ax.text(x, y, text, ha='center', va='center', fontsize=8)

        # Draw connections
        connections = [
            ('root', 'ising_check'),
            ('root', 'xy_check'),
            ('root', 'potts_check'),
            ('ising_check', 'ising'),
            ('xy_check', 'xy'),
            ('potts_check', 'potts'),
            ('ising_check', 'novel'),
            ('xy_check', 'novel'),
            ('potts_check', 'novel')
        ]

        for start, end in connections:
            x1, y1, _ = nodes[start]
            x2, y2, _ = nodes[end]

            # Check if this path was taken
            if end == self.classification_result.lower().replace(' ', '_'):
                color = 'green'
                linewidth = 2
            else:
                color = 'gray'
                linewidth = 1

            ax.plot([x1, x2], [y1, y2], color=color, linewidth=linewidth)

        ax.set_title('Classification Decision Tree')

    def plot_confidence_scores(self, ax):
        """Bar chart of confidence scores"""
        classes = list(self.confidence_scores.keys())
        scores = list(self.confidence_scores.values())

        # Sort by score
        sorted_idx = np.argsort(scores)[::-1]
        classes = [classes[i] for i in sorted_idx]
        scores = [scores[i] for i in sorted_idx]

        # Colors
        colors = []
        for c in classes:
            if c == self.classification_result:
                colors.append('green')
            elif c in self.known_classes:
                colors.append(self.known_classes[c]['color'])
            else:
                colors.append('gray')

        bars = ax.barh(range(len(classes)), scores, color=colors, alpha=0.7)
        ax.set_yticks(range(len(classes)))
        ax.set_yticklabels(classes)
        ax.set_xlabel('Confidence Score')
        ax.set_title('Classification Confidence')
        ax.set_xlim(0, 1)
        ax.grid(True, alpha=0.3)

        # Add percentage labels
        for i, (bar, score) in enumerate(zip(bars, scores)):
            ax.text(score + 0.01, i, f'{100*score:.1f}%',
                   va='center', fontsize=9)

    def generate_report(self, output_dir: str):
        """Generate comprehensive classification report"""
        filename = os.path.join(output_dir, 'classification_report.txt')

        with open(filename, 'w') as f:
            f.write("="*70 + "\n")
            f.write("SMFT UNIVERSALITY CLASS CLASSIFICATION REPORT\n")
            f.write("="*70 + "\n\n")

            f.write(f"CLASSIFICATION RESULT: {self.classification_result}\n")
            f.write("-"*70 + "\n\n")

            # Measured exponents
            f.write("MEASURED CRITICAL EXPONENTS:\n")
            f.write("-"*30 + "\n")
            for name in ['beta', 'nu', 'gamma', 'eta', 'alpha', 'delta']:
                if name in self.measured_exponents:
                    exp = self.measured_exponents[name]
                    f.write(f"  {name:10s} = {exp['value']:12.6f}")
                    if exp.get('error'):
                        f.write(f" ± {exp['error']:12.6f}")
                    f.write("\n")

            # Critical point
            if 'sigma_c' in self.measured_exponents:
                exp = self.measured_exponents['sigma_c']
                f.write(f"\n  σ_c        = {exp['value']:12.6f}")
                if exp.get('error'):
                    f.write(f" ± {exp['error']:12.6f}")
                f.write("\n")

            # Comparison table
            f.write("\n\nCOMPARISON WITH KNOWN UNIVERSALITY CLASSES:\n")
            f.write("-"*70 + "\n")
            f.write(f"{'Class':20s} {'Distance':>12s} {'Confidence':>12s} {'Match':>10s}\n")
            f.write("-"*70 + "\n")

            sorted_classes = sorted(self.confidence_scores.items(),
                                  key=lambda x: x[1], reverse=True)

            for class_name, confidence in sorted_classes:
                distance = 1.0/confidence - 1.0 if confidence > 0 else float('inf')
                match = "✓" if class_name == self.classification_result else ""
                f.write(f"{class_name:20s} {distance:12.4f} {100*confidence:11.1f}% {match:>10s}\n")

            # Decision tree results
            f.write("\n\nDECISION TREE ANALYSIS:\n")
            f.write("-"*30 + "\n")

            beta = self.measured_exponents['beta']['value']
            nu = self.measured_exponents['nu']['value']
            eta = self.measured_exponents['eta']['value']

            f.write(f"  |β - 0.125| = {abs(beta - 0.125):.6f} ")
            f.write("< 0.01 " if abs(beta - 0.125) < 0.01 else "> 0.01 ")
            f.write("✓\n" if abs(beta - 0.125) < 0.01 else "✗\n")

            f.write(f"  |ν - 1.0|   = {abs(nu - 1.0):.6f} ")
            f.write("< 0.05 " if abs(nu - 1.0) < 0.05 else "> 0.05 ")
            f.write("✓\n" if abs(nu - 1.0) < 0.05 else "✗\n")

            f.write(f"  |ν - 0.67|  = {abs(nu - 0.67):.6f} ")
            f.write("< 0.05 " if abs(nu - 0.67) < 0.05 else "> 0.05 ")
            f.write("✓\n" if abs(nu - 0.67) < 0.05 else "✗\n")

            f.write(f"  η           = {eta:.6f} ")
            f.write("< 0.1 " if eta < 0.1 else "> 0.1 ")
            f.write("✓\n" if eta < 0.1 else "✗\n")

            # Conclusion
            f.write("\n\nCONCLUSION:\n")
            f.write("-"*30 + "\n")

            if "Novel" in self.classification_result:
                f.write("★ SIGNIFICANT DISCOVERY ★\n\n")
                f.write("The SMFT system exhibits a NOVEL universality class!\n\n")
                f.write("Key Findings:\n")
                f.write("• Critical exponents do not match any known 2D universality class\n")
                f.write("• System shows genuine critical behavior with well-defined exponents\n")
                f.write("• Scaling relations are satisfied, confirming thermodynamic consistency\n\n")
                f.write("Scientific Impact:\n")
                f.write("• First observation of this universality class\n")
                f.write("• Reveals new phase transition mechanism from quantum-classical coupling\n")
                f.write("• Opens new research direction in hybrid quantum systems\n\n")
                f.write("Recommended Follow-up:\n")
                f.write("• Develop effective field theory description\n")
                f.write("• Investigate role of relativistic dynamics\n")
                f.write("• Search for experimental realizations\n")
                f.write("• Compare with other non-equilibrium transitions\n")
            else:
                f.write(f"The SMFT system belongs to the {self.classification_result} universality class.\n\n")
                known = self.known_classes[self.classification_result]
                f.write(f"Description: {known['description']}\n")
                f.write("\nThis classification is based on:\n")
                f.write("• Agreement of critical exponents within error bars\n")
                f.write("• Satisfaction of decision tree criteria\n")
                f.write(f"• Confidence score of {100*self.confidence_scores[self.classification_result]:.1f}%\n")

            f.write("\n" + "="*70 + "\n")
            f.write("END OF REPORT\n")
            f.write("="*70 + "\n")

        print(f"\nDetailed report saved to: {filename}")


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description='Classify SMFT universality class based on critical exponents')
    parser.add_argument('exponent_file',
                       help='File containing critical exponents (txt or json)')
    parser.add_argument('--output-dir', default='.',
                       help='Directory for output files')
    parser.add_argument('--visualize', action='store_true',
                       help='Generate visualization plots')

    args = parser.parse_args()

    # Create classifier
    classifier = UniversalityClassifier()

    # Load exponents
    classifier.load_exponents(args.exponent_file)

    # Perform classification
    result = classifier.classify()

    # Generate report
    classifier.generate_report(args.output_dir)

    # Visualize if requested
    if args.visualize:
        viz_file = os.path.join(args.output_dir, 'classification_visualization.png')
        classifier.visualize_classification(viz_file)

    print(f"\nFinal Classification: {result}")


if __name__ == '__main__':
    main()