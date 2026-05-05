# KitBit

KitBit is a Python-based sequence prediction framework that reproduces and extends the algorithm from:

> **"KitBit: A New AI Model for Solving Intelligence Tests and Numerical Series"**
> V. Corsino, J. M. Gilpérez, and L. Herrera — *IEEE Transactions on Pattern Analysis and Machine Intelligence*, Vol. 45, No. 11, 2023. ([arXiv:2206.08965](https://arxiv.org/abs/2206.08965))

The original KitBit is a rule-based, training-free model that discovers patterns in numerical sequences by chaining symbolic operations called **kitas** via Breadth-First Search (BFS), guided by the Minimum Description Length (MDL) principle. This project reproduces that baseline and introduces four enhancements: composite sequence decomposition, heuristic-based search ordering, a prime number kita, and dynamic early stopping.

---

## Project Structure

```
kitbit/
├── KitBit_Algorithms.py       # Main entry point; defines and runs all solver pipelines
├── data/
│   ├── data.py                # Sequence datasets (sr0, sr1, kl2, composite_test_set)
│   └── prime_data.py          # Prime number sequence test set
├── helpers/
│   ├── generic.py             # Core KitBit and DynamicKitBit solver classes (BFS)
│   ├── composite_helper.py    # Decomposition, sanity checks, and composite prediction
│   ├── prime_helper.py        # Prime-specific prediction and fallback logic
│   └── heuristic.py          # Heuristic search runner
├── results/                   # Output files from current solver runs
├── results_old/               # Archived results from earlier experiments
└── output                     # General output log
```

---

## Running the Code

Requires Python 3 only — no external packages needed.

```bash
python KitBit_Algorithms.py
```

Results are saved to the `results/` directory and accuracy is printed to stdout:

```
[Pipeline Name] Total: N | Solved: M | Accuracy: X.XX%
```

---

## Datasets

| Dataset | Variable | Description |
|---|---|---|
| IQ series | `sr0` | IQ-test style integer sequences |
| Literature series | `sr1` | Sequences from published literature |
| OEIS | loaded from `data/` | 341,553 sequences from the Online Encyclopedia of Integer Sequences |
| Composite test set | `composite_test_set` | 158 sequences formed by interleaving two or more sub-patterns |
| Prime test set | `prime_test_set` | Sequences related to prime numbers |

---

## Reference

V. Corsino, J. M. Gilpérez, and L. Herrera, "KitBit: A New AI Model for Solving Intelligence Tests and Numerical Series," *IEEE Transactions on Pattern Analysis and Machine Intelligence*, vol. 45, no. 11, 2023. [arXiv:2206.08965](https://arxiv.org/abs/2206.08965)