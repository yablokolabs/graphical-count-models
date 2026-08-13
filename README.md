# Graphical Count Models (Mojo)

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

A Mojo library implementing the four graphical count distributions introduced in
[Graphical Models for Multivariate Count Data (arXiv:2608.11366)](https://arxiv.org/abs/2608.11366).

## Features

- **mult G** – Graphical multinomial (with‑replacement, fixed draws)
- **nm G** – Graphical negative‑multinomial (with‑replacement, fixed failures)
- **hg G** – Graphical hypergeometric (without‑replacement, fixed draws)
- **nhg G** – Graphical negative‑hypergeometric (without‑replacement, fixed failures)

All models are defined on **decomposable (chordal) graphs** and satisfy the global Markov property with respect to the graph.

The library provides:

- Exact probability mass functions (PMFs) using graph‑multinomial coefficients.
- Simple samplers for each model.
- Conjugate Bayesian hierarchies (G‑Dirichlet / G‑inverted‑Dirichlet priors) with closed‑form posteriors.
- A small example using the Rydberg‑atom excitation data from the paper.

## Installation

```bash
# Clone the repository
git clone https://github.com/your‑username/graphical-count-models.git
cd graphical-count-models

# Build the Mojo library (requires Mojo 1.0+)
mojo build src/graphical_counts.mojo
```

The library can also be installed in editable mode for development:

```bash
pip install -e .
```

## Quick start

```mojo
from graphical_counts import GraphicalMultinomial

# Define a simple decomposable graph (a chain 0‑1‑2‑3)
let G = [(0,1), (1,2), (2,3)]

# Total number of draws
let r = 10

# Create the model
let model = GraphicalMultinomial(G, r: r)

# Observed counts on the vertices (dictionary vertex → count)
let observed = {0:5, 1:2, 2:0, 3:1}

# Fit the model (maximum‑likelihood via EM)
let fitted = model.fit(observed)

# Predict a new count vector for a different total draw size
let pred = fitted.predict(r: 12)
print(pred)
```

See `examples/rydberg_demo.mojo` for a full notebook‑style walkthrough.

## Documentation

Full API documentation is generated with Sphinx and hosted on GitHub Pages.

## Contributing

Contributions are welcome! Please open an issue or submit a pull request.

## License

MIT License – see the `LICENSE` file for details.