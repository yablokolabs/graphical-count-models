# Graphical Count Models in Mojo

[![CI](https://github.com/yablokolabs/graphical-count-models/actions/workflows/ci.yml/badge.svg)](https://github.com/yablokolabs/graphical-count-models/actions/workflows/ci.yml)
[![Mojo](https://img.shields.io/badge/Mojo-1.0.0-orange)](https://mojolang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

A tested, dependency-free **small-graph reference implementation** of the four graphical count distributions in Danielewska and Kołodziejek, [*Graphical Models for Multivariate Count Data*](https://arxiv.org/abs/2608.11366) (arXiv:2608.11366v1).

The code computes the graph polynomials, the two graph-multinomial coefficient families, and exact PMFs for:

| Stopping rule | With replacement | Without replacement |
|---|---|---|
| Fixed number of draws | graphical multinomial $`\operatorname{mult}_G`$ | graphical hypergeometric $`\operatorname{hg}_G`$ |
| Fixed number of failures | graphical negative multinomial $`\operatorname{nm}_G`$ | graphical negative hypergeometric $`\operatorname{nhg}_G`$ |

> [!IMPORTANT]
> This repository is an educational and verification-oriented implementation for **small decomposable graphs**. It uses exhaustive independent-set enumeration and bounded coefficient grids. It is not a large-scale inference package and does not currently implement fitting, graph learning, Bayesian posteriors, or production samplers.

## What was checked against the paper

The implementation follows arXiv:2608.11366v1:

- graph polynomials: Definition 2.1;
- graph-multinomial coefficients: Definition 2.2;
- $`\operatorname{mult}_G`$ and $`\operatorname{nm}_G`$: Definition 3.1;
- $`\operatorname{hg}_G`$: Definition 4.1;
- $`\operatorname{nhg}_G`$: Definition 4.5.

It also enforces the paper's decomposability restriction for the four general PMFs. The exceptional graphical Bernoulli case $`\operatorname{mult}_G(1,\mathbf y)`$ is mathematically valid for arbitrary finite graphs (Remark 3.2), but the high-level PMF API deliberately keeps one uniform decomposable-graph contract in this initial release.

## Mathematics

Let $`G=(V,E)`$ be a finite simple graph and $`G^*`$ its complement. Let $`\mathcal C_{G^*}`$ be **all** cliques of $`G^*`$, including $`\varnothing`$; equivalently these are all independent sets of $`G`$.

### Graph polynomials

Definition 2.1 gives

```math
\delta_G(\mathbf y)
=\sum_{C\in\mathcal C_{G^*}}\prod_{v\in C} y_v,
\qquad
\Delta_G(\mathbf x)
=\sum_{C\in\mathcal C_{G^*}}(-1)^{|C|}\prod_{v\in C}x_v
=\delta_G(-\mathbf x).
```

The empty clique contributes the constant term $`1`$. This matters: enumerating only maximal cliques does **not** compute either polynomial.

### Two coefficient families

The first-type coefficients are defined by

```math
\delta_G(\mathbf y)^r
=\sum_{\mathbf n\in\mathbb N^V}
\binom{r}{\mathbf n}_G\mathbf y^{\mathbf n},
\qquad r\in\mathbb N_+.
```

The second-type coefficients are defined by the formal power series

```math
\Delta_G(\mathbf x)^{-s}
=\sum_{\mathbf n\in\mathbb N^V}
\left[\begin{matrix}|\mathbf n|+s-1\\\mathbf n\end{matrix}\right]_G
\mathbf x^{\mathbf n},
\qquad s>0.
```

This implementation computes both coefficients directly with truncated multivariate recurrences. For decomposable graphs, the paper also gives clique-separator product formulas (Lemma 2.3):

```math
\binom{r}{\mathbf n}_G
=
\frac{\displaystyle\prod_{C\in\mathcal C_G^+}\binom{r}{\mathbf n_C}}
{\displaystyle\prod_{S\in\mathcal S_G^-}\binom{r}{\mathbf n_S}^{\nu_S}},
```

```math
\left[\begin{matrix}|\mathbf n|+s-1\\\mathbf n\end{matrix}\right]_G
=
\frac{\displaystyle\prod_{C\in\mathcal C_G^+}
\binom{|\mathbf n_C|+s-1}{\mathbf n_C}}
{\displaystyle\prod_{S\in\mathcal S_G^-}
\binom{|\mathbf n_S|+s-1}{\mathbf n_S}^{\nu_S}}.
```

### The four PMFs

For a decomposable graph, Definition 3.1 gives

```math
\Pr(\mathbf N=\mathbf n)
=
\binom{r}{\mathbf n}_G
\mathbf y^{\mathbf n}\,
\delta_G(\mathbf y)^{-r},
\qquad
\mathbf N\sim\operatorname{mult}_G(r,\mathbf y).
```

For

```math
M_G=\left\{\mathbf x\in(0,\infty)^V:
\Delta_{G_A}(\mathbf x_A)>0\text{ for every }A\subseteq V\right\},
```

the graphical negative multinomial is

```math
\Pr(\mathbf N=\mathbf n)
=
\left[\begin{matrix}|\mathbf n|+r-1\\\mathbf n\end{matrix}\right]_G
\mathbf x^{\mathbf n}\,
\Delta_G(\mathbf x)^r,
\qquad
\mathbf N\sim\operatorname{nm}_G(r,\mathbf x).
```

The implementation checks the full induced-subgraph condition defining $`M_G`$.

Definition 4.1 gives the graphical hypergeometric PMF

```math
\Pr(\mathbf N=\mathbf n)
=
\frac{
\binom{r}{\mathbf n}_G
\binom{M-r}{\mathbf K-\mathbf n}_G
}{
\binom{M}{\mathbf K}_G
},
\qquad
\mathbf N\sim\operatorname{hg}_G(M,\mathbf K,r).
```

Its exact support is

```math
S^{\operatorname{hg}_G}_{M,\mathbf K,r}
=
\left\{\mathbf n\in\mathbb N^V:
 n_v\le K_v\ \forall v,
\ 0\le r-|\mathbf n_C|\le M-|\mathbf K_C|\ \forall C\in\mathcal C_G
\right\}.
```

Definition 4.5 gives the graphical negative hypergeometric PMF

```math
\Pr(\mathbf N=\mathbf n)
=
\frac{
\left[\begin{matrix}|\mathbf n|+r-1\\\mathbf n\end{matrix}\right]_G
\left[\begin{matrix}M-r+|\mathbf K|-|\mathbf n|\\\mathbf K-\mathbf n\end{matrix}\right]_G
}{
\left[\begin{matrix}M+|\mathbf K|\\\mathbf K\end{matrix}\right]_G
},
\qquad
\mathbf N\sim\operatorname{nhg}_G(M,\mathbf K,r),
```

with

```math
S^{\operatorname{nhg}_G}_{M,\mathbf K,r}
=\{\mathbf n\in\mathbb N^V:0\le n_v\le K_v\ \forall v\}.
```

### Conditioning relationships

The without-replacement models arise by conditioning independent with-replacement models on their sum:

```math
\mathbf N_1\mid(\mathbf N_1+\mathbf N_2=\mathbf K)
\sim\operatorname{hg}_G(M,\mathbf K,r),
```

when $`\mathbf N_1\sim\operatorname{mult}_G(r,\mathbf y)`$ and $`\mathbf N_2\sim\operatorname{mult}_G(M-r,\mathbf y)`$; and

```math
\mathbf N_1\mid(\mathbf N_1+\mathbf N_2=\mathbf K)
\sim\operatorname{nhg}_G(M,\mathbf K,r),
```

when $`\mathbf N_1\sim\operatorname{nm}_G(r,\mathbf x)`$ and $`\mathbf N_2\sim\operatorname{nm}_G(M-r+1,\mathbf x)`$ (Propositions 4.2 and 4.6).

### Classical limits

These provide useful correctness checks:

- **complete $`G`$:** the four families reduce to the classical multivariate multinomial, negative multinomial, hypergeometric, and negative hypergeometric distributions;
- **empty $`G`$:** they factor into products of the corresponding univariate laws.

The test suite exercises both limits and checks that the finite-support PMFs normalize to one.

## Install and run

[Pixi](https://pixi.sh/) is the recommended reproducible environment manager.

```bash
git clone https://github.com/yablokolabs/graphical-count-models.git
cd graphical-count-models
pixi run check
```

The tasks are:

```bash
pixi run format       # format Mojo sources
pixi run format-check # verify formatting
pixi run build        # precompile package with warnings as errors
pixi run test         # run six numerical tests
pixi run example      # run the complete-graph example
```

If Mojo 1.0 is already installed, the direct commands are:

```bash
mojo precompile src/graphical_counts -o graphical_counts.mojoc --Werror
mojo run --Werror -I src tests/test_graphical_counts.mojo
mojo run --Werror -I src examples/basic.mojo
```

## API example

Vertices must be contiguous integers `0..vertex_count-1`; graphs use a flattened symmetric adjacency matrix.

```mojo
from graphical_counts import (
    add_edge,
    empty_adjacency,
    graphical_multinomial_pmf,
)


def main() raises:
    # Chain 0--1--2, a decomposable graph.
    var adjacency = empty_adjacency(3)
    add_edge(adjacency, 3, 0, 1)
    add_edge(adjacency, 3, 1, 2)

    var y: List[Float64] = [0.5, 1.0, 1.5]
    var counts: List[Int] = [1, 0, 1]
    var probability = graphical_multinomial_pmf(
        adjacency, 3, 2, y, counts
    )
    print(probability)
```

See [`examples/basic.mojo`](examples/basic.mojo) for a runnable example of $`\operatorname{mult}_G`$ and $`\operatorname{hg}_G`$.

## Practical use cases

The paper's core modeling pattern is: vertices represent candidate events or objects; an edge in $`G`$ means the pair is incompatible; one admissible configuration is an independent set of $`G`$ (a clique of $`G^*`$); the observed vector aggregates vertex occurrences across configurations.

| Domain | Graph construction | Useful question | Candidate family |
|---|---|---|---|
| **Rydberg-atom blockade experiments** | atoms are vertices; blockade interactions are edges | Distribution of blockade-consistent excitation patterns | $`\operatorname{mult}_G(1,\mathbf y)`$; $`\operatorname{hg}_G`$ for subsets of a fixed admissible collection |
| **Wireless/CSMA scheduling** | links are vertices; interference creates edges | Link activity counts across approximately independent stationary schedule snapshots | $`\operatorname{mult}_G`$; $`\operatorname{hg}_G`$ for a subset of fixed schedule history |
| **Cancer genomics** | genes/pathways are vertices; hypothesized mutual exclusions are edges | Alteration prevalence under an idealized exclusion null, or a finite-cohort subset with fixed totals | $`\operatorname{mult}_G`$ / $`\operatorname{hg}_G`$ |
| **Spatial hard-core systems** | locations/objects are vertices; overlap or minimum-distance conflicts are edges | Occupancy counts across repeated admissible spatial configurations | $`\operatorname{mult}_G`$; $`\operatorname{hg}_G`$ for conditioned finite collections |
| **Adsorption and packing** | candidate placements are vertices; geometric conflicts are edges | Counts accumulated before a rejection/jamming event under an idealized stopped mechanism | $`\operatorname{nm}_G`$ / $`\operatorname{nhg}_G`$ |
| **Loss networks and resource sharing** | calls/routes/jobs are vertices; shared-capacity conflicts are edges | Occupancy counts over independent stationary snapshots | $`\operatorname{mult}_G`$ / $`\operatorname{hg}_G`$ |
| **Abstract polymer models** | polymers are vertices; incompatibility creates edges | Weighted admissible configurations and repeated polymer occupancy counts | $`\operatorname{mult}_G(1,\mathbf y)`$ and $`\operatorname{mult}_G`$ |
| **Conflict-constrained portfolios** *(extension)* | candidate positions are vertices; prohibited co-holdings are edges | Exact probability calculations for small rule-constrained allocation snapshots | $`\operatorname{mult}_G`$ / $`\operatorname{hg}_G`$ |
| **Maintenance or deployment windows** *(extension)* | tasks are vertices; unsafe concurrent tasks are edges | Counts across feasible windows or subsets of a fixed rollout plan | $`\operatorname{mult}_G`$ / $`\operatorname{hg}_G`$ |
| **Experiment design with incompatible treatments** *(extension)* | treatments are vertices; co-assignment exclusions are edges | Conditional allocation probabilities under fixed aggregate exposure | $`\operatorname{hg}_G`$ |

The first seven rows are discussed in Sections 7–8 of the paper. The last three are plausible extensions of the same independent-set mechanism, not claims validated by the paper.

### Important modeling caveats

- The paper treats graph-constrained support as an **idealization**. Out-of-support observations can indicate noise, imperfect exclusion, omitted interactions, or a misspecified graph.
- Cancer mutual exclusivity is often a tendency rather than a hard law; use the model as a conditional law or idealized null where appropriate.
- Consecutive CSMA schedules are generally dependent. The exact $`\operatorname{mult}_G`$ interpretation requires independent stationary samples; widely spaced observations may only approximate it.
- Negative-model sampling requires the specific Cartier–Foata admissibility and transition mechanism from the paper. A generic “count until any failure” process does not automatically have an $`\operatorname{nm}_G`$ or $`\operatorname{nhg}_G`$ law.
- The general four-family construction requires a decomposable graph. Do not silently apply the coefficient PMFs to a non-chordal graph and infer the global Markov property.

## Implementation scope and complexity

- `vertex_count <= 16`;
- graph polynomial/domain checks enumerate up to $`2^{|V|}`$ subsets;
- coefficient calls allocate $`\prod_v(n_v+1)`$ states, capped at 2,000,000;
- calculations use `Float64` and are not log-space stabilized;
- PMF calls prioritize traceability to Definitions 2.1–4.5 over asymptotic speed;
- no random sampler, estimator, graph learner, data loader, posterior/predictive layer, or Rydberg dataset is bundled.

For larger systems, implement clique-separator/DAG factorization, sparse polynomial representations, log-gamma arithmetic, and dedicated sampling/inference algorithms before treating this as production software.

## Citation

```bibtex
@misc{danielewska2026graphical,
  title         = {Graphical Models for Multivariate Count Data},
  author        = {Iza Danielewska and Bartosz Ko{\l}odziejek},
  year          = {2026},
  eprint        = {2608.11366},
  archivePrefix = {arXiv},
  primaryClass  = {stat.ME},
  doi           = {10.48550/arXiv.2608.11366}
}
```

This repository is an independent implementation and is not affiliated with the paper's authors.

## License

[MIT](LICENSE)
