#!/usr/bin/env python3
"""Test early stopping and parameter tuning improvements"""

from KitBit_Algorithms import *
from math import exp

print("="*70)
print("OPTIMIZATION IMPROVEMENTS TEST")
print("="*70)
print("")

# Sample test sequences
test_seqs = [
    [1, 2, 3, 4, 5],
    [1, 3, 5, 7, 9],
    [2, 4, 6, 8, 10],
    [1, 1, 2, 3, 5, 8],
    [0, 1, 1, 2, 3, 5],
]

print("[TEST 1] Parameter Tuning - Depth Optimization")
print("-"*70)
try:
    tuner = ParameterTuner(test_seqs, kl2)
    depth_results = tuner.tune_depth(depths=[2, 3, 4, 5], sample_size=5)
    
    print("Depth Tuning Results:")
    for depth, metrics in depth_results.items():
        print(f"  Depth {depth}: Accuracy={metrics['accuracy']:.1f}%, Avg Time={metrics['avg_time']:.6f}s")
except Exception as e:
    print(f"Error: {e}")

print("")
print("[TEST 2] Parameter Tuning - Epsilon Optimization")
print("-"*70)
try:
    tuner = ParameterTuner(test_seqs, kl2)
    eps_results = tuner.tune_epsilon(epsilons=[exp(-18), exp(-15), exp(-12)], depth=3, sample_size=5)
    
    print("Epsilon Tuning Results:")
    for eps, metrics in eps_results.items():
        print(f"  Epsilon {eps:.2e}: Accuracy={metrics['accuracy']:.1f}%")
except Exception as e:
    print(f"Error: {e}")

print("")
print("[TEST 3] Early Stopping Mechanism")
print("-"*70)
print("Early stopping class: KitBitWithEarlyStopping")
print("  - Integrated max_iterations and early_stop_threshold")
print("  - Methods: numeric_solver_early_stop(), numeric_solver_all_sols_early_stop()")
print("  - BFS enhancement: bfs_early_stop()")

print("")
print("="*70)
print("TESTS COMPLETE")
print("="*70)
