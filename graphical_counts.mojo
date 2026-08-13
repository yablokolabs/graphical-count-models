// graphical_counts.mojo
// Implementation of the four graphical count distributions from
// Danielewska & Kołodziejek, arXiv:2608.11366.
//
// The code is deliberately concise; full numerical stability checks,
// extensive comments, and edge‑case handling can be added later.

import Python

# ----------------------------------------------------------------------
# Graph utilities
# ----------------------------------------------------------------------
struct Graph:
    # Simple adjacency list representation.
    var adj: Dictionary[Int, List[Int]]

    fn __init__(inout self, vertices: List[Int], edges: List[(Int, Int)]) -> None:
        self.adj = {}
        for v in vertices:
            self.adj[v] = []
        for (u, v) in edges:
            self.adj[u].append(v)
            self.adj[v].append(u)

    # Return list of maximal cliques (Bron–Kerbosch with pivot)
    fn maximal_cliques(self) -> List[List[Int]]:
        # Placeholder: use Python networkx for maximal cliques.
        let nx = Python.import_module("networkx")
        let G = nx.Graph()
        for v in self.adj.keys():
            G.add_node(v)
        for (u, nbrs) in self.adj.items():
            for w in nbrs:
                if u < w:
                    G.add_edge(u, w)
        let cliques = nx.find_cliques(G)
        return List[List[Int]](cliques)

    # Check if graph is chordal (decomposable)
    fn is_chordal(self) -> Bool:
        let nx = Python.import_module("networkx")
        let G = nx.Graph()
        for v in self.adj.keys():
            G.add_node(v)
        for (u, nbrs) in self.adj.items():
            for w in nbrs:
                if u < w:
                    G.add_edge(u, w)
        return nx.is_chordal(G)

# ----------------------------------------------------------------------
# Independence polynomial δ_G(y) and signed polynomial Δ_G(x)
# ----------------------------------------------------------------------
fn independence_polynomial(g: Graph, y: List[Float]) -> Float:
    # Compute sum over all cliques of complement graph G*.
    # We'll use Python networkx to get complement cliques.
    let nx = Python.import_module("networkx")
    let G = nx.Graph()
    for v in g.adj.keys():
        G.add_node(v)
    for (u, nbrs) in g.adj.items():
        for w in nbrs:
            if u < w:
                G.add_edge(u, w)
    let Gc = nx.complement(G)
    var sum_val = 0.0
    for clique in Gc.find_cliques():
        var term = 1.0
        for v in clique:
            term *= y[v]
        sum_val += term
    return sum_val

fn signed_polynomial(g: Graph, x: List[Float]) -> Float:
    # Δ_G(x) = Σ_{C∈C(G*)} (-1)^{|C|} ∏_{v∈C} x_v
    let nx = Python.import_module("networkx")
    let G = nx.Graph()
    for v in g.adj.keys():
        G.add_node(v)
    for (u, nbrs) in g.adj.items():
        for w in nbrs:
            if u < w:
                G.add_edge(u, w)
    let Gc = nx.complement(G)
    var sum_val = 0.0
    for clique in Gc.find_cliques():
        var term = 1.0
        for v in clique:
            term *= x[v]
        sum_val += ((-1) ** len(clique)) * term
    return sum_val

# ----------------------------------------------------------------------
# Graph‑multinomial coefficients (first and second type)
# ----------------------------------------------------------------------
fn graph_multinomial_coeff_first(g: Graph, r: Int, n: List[Int]) -> Float:
    # δ_G(y)^r coefficient of ∏ y_v^{n_v}
    # We'll expand using Python sympy for brevity.
    let sympy = Python.import_module("sympy")
    let y = sympy.symbols(' '.join([f'y{i}' for i in range(len(n))]))
    let poly = sympy.expand((independence_polynomial(g, y)) ** r)
    var coeff = 0.0
    for i in range(len(n)):
        let mon = sympy.Symbol(f'y{i}') ** n[i]
        coeff = sympy.expand(poly).coeff(mon)
    return Float(coeff)

fn graph_multinomial_coeff_second(g: Graph, s: Float, n: List[Int]) -> Float:
    # Coefficient of ∏ x_v^{n_v} in (Δ_G(x) - s)^{-1}
    # Use series expansion via sympy.
    let sympy = Python.import_module("sympy")
    let x = sympy.symbols(' '.join([f'x{i}' for i in range(len(n))]))
    let poly = sympy.series((signed_polynomial(g, x) - s), x, 0, sum(n) + 1).removeO()
    var coeff = 0.0
    for i in range(len(n)):
        let mon = sympy.Symbol(f'x{i}') ** n[i]
        coeff = sympy.expand(poly).coeff(mon)
    return Float(coeff)

# ----------------------------------------------------------------------
# Distribution classes
# ----------------------------------------------------------------------
struct GraphicalMultinomial:
    let g: Graph
    let r: Int
    let y: List[Float]

    fn pmf(self, n: List[Int]) -> Float:
        if not self.g.is_chordal():
            raise "Graph must be decomposable for GraphicalMultinomial."
        let coeff = graph_multinomial_coeff_first(self.g, self.r, n)
        let denom = independence_polynomial(self.g, self.y) ** self.r
        return coeff / denom

    fn sample(self) -> List[Int]:
        # Simple rejection sampling: draw a clique from complement according to weights, repeat r times.
        let nx = Python.import_module("networkx")
        let G = nx.Graph()
        for v in self.adj.keys():
            G.add_node(v)
        for (u, nbrs) in self.adj.items():
            for w in nbrs:
                if u < w:
                    G.add_edge(u, w)
        let Gc = nx.complement(G)
        # Compute weights proportional to ∏ y_v^{1_{v∈C}}
        var weights = {}
        var total = 0.0
        for clique in Gc.find_cliques():
            var w = 1.0
            for v in clique:
                w *= self.y[v]
            weights[clique] = w
            total += w
        # Normalize
        for k in weights.keys():
            weights[k] = weights[k] / total
        # Sample r cliques
        var counts = List[Int](repeating: 0, count: len(self.g.adj.keys()))
        for _ in range(self.r):
            # Choose a clique according to weights
            let rand = Python.random.random()
            var cum = 0.0
            var chosen = None[List[Int]]
            for (clique, prob) in weights.items():
                cum += prob
                if cum >= rand:
                    chosen = clique
                    break
            # Increment counts
            for v in chosen!:
                counts[v] += 1
        return counts

# ----------------------------------------------------------------------
# Graphical Hypergeometric (conditional on sum of two mult G draws)
# ----------------------------------------------------------------------
struct GraphicalHypergeometric:
    let g: Graph
    let M: Int
    let K: List[Int]   # total content vector
    let r: Int

    fn pmf(self, n: List[Int]) -> Float:
        if not self.g.is_chordal():
            raise "Graph must be decomposable for GraphicalHypergeometric."
        let coeff_num = graph_multinomial_coeff_first(self.g, self.r, n)
        let coeff_den = graph_multinomial_coeff_first(self.g, self.M - self.r, List[Int](K = self.K - n))
        let coeff_total = graph_multinomial_coeff_first(self.g, self.M, self.K)
        return coeff_num * coeff_den / coeff_total

    fn sample(self) -> List[Int]:
        # Sample two independent mult G draws and condition on their sum.
        let y = List[Float](repeating: 1.0, count: len(self.g.adj.keys()))
        let mult = GraphicalMultinomial(g: self.g, r: self.r, y: y)
        let mult2 = GraphicalMultinomial(g: self.g, r: self.M - self.r, y: y)
        let n1 = mult.sample()
        let n2 = mult2.sample()
        # If sum matches K, accept; otherwise repeat (simple rejection)
        while List[Int](K = n1 + n2) != self.K:
            n1 = mult.sample()
            n2 = mult2.sample()
        return n1

# ----------------------------------------------------------------------
# Graphical Negative Multinomial (failure‑stopped)
# ----------------------------------------------------------------------
struct GraphicalNegativeMultinomial:
    let g: Graph
    let r: Float
    let x: List[Float]

    fn pmf(self, n: List[Int]) -> Float:
        if not self.g.is_chordal():
            raise "Graph must be decomposable for GraphicalNegativeMultinomial."
        let coeff = graph_multinomial_coeff_second(self.g, self.r, n)
        let denom = signed_polynomial(self.g, self.x) ** self.r
        return coeff / denom

    fn sample(self) -> List[Int]:
        # Use the Markov‑chain construction described in the paper.
        # For brevity we fall back to rejection sampling using the pmf.
        var n = List[Int](repeating: 0, count: len(self.g.adj.keys()))
        # Simple bound: iterate over plausible support (very naive)
        let max_iter = 10000
        var attempts = 0
        while attempts < max_iter:
            # Propose a count vector by drawing from a Poisson with mean r * x_i
            var proposal = List[Int](repeating: 0, count: len(self.x))
            for i in range(len(self.x)):
                let lam = self.r * self.x[i]
                let k = Python.random.poisson(lam)
                proposal[i] = k
            let prob = self.pmf(proposal)
            if Python.random.random() < prob / (self.pmf(List[Int](repeating: 0, count: len(self.x)))) :
                return proposal
            attempts += 1
        raise "Sampling failed after many attempts."

# ----------------------------------------------------------------------
# Graphical Negative Hypergeometric (conditional version)
# ----------------------------------------------------------------------
struct GraphicalNegativeHypergeometric:
    let g: Graph
    let M: Int
    let K: List[Int]
    let r: Float

    fn pmf(self, n: List[Int]) -> Float:
        if not self.g.is_chordal():
            raise "Graph must be decomposable for GraphicalNegativeHypergeometric."
        let coeff_num = graph_multinomial_coeff_second(self.g, self.r, n)
        let coeff_den = graph_multinomial_coeff_second(self.g, self.M - self.r + 1, List[Int](K = self.K - n))
        let coeff_total = graph_multinomial_coeff_second(self.g, self.M + 1, self.K)
        return coeff_num * coeff_den / coeff_total

    fn sample(self) -> List[Int]:
        # Sample two independent negative multinomial draws and condition.
        let nm = GraphicalNegativeMultinomial(g: self.g, r: self.r, x: List[Float](repeating: 1.0, count: len(self.g.adj.keys())))
        let nm2 = GraphicalNegativeMultinomial(g: self.g, r: self.M - self.r + 1, x: List[Float](repeating: 1.0, count: len(self.g.adj.keys())))
        let n1 = nm.sample()
        let n2 = nm2.sample()
        while List[Int](K = n1 + n2) != self.K:
            n1 = nm.sample()
            n2 = nm2.sample()
        return n1

# ----------------------------------------------------------------------
# Bayesian priors (placeholder)
# ----------------------------------------------------------------------
struct GDirichlet:
    let g: Graph
    let alpha: List[Float]
    let beta: Float

    # Placeholder: return prior density (not implemented)
    fn density(self, theta: List[Float]) -> Float:
        # In a full implementation we would compute Δ_G(theta) etc.
        return 1.0

# ----------------------------------------------------------------------
# End of library
# ----------------------------------------------------------------------