"""Exact small-graph reference routines for arXiv:2608.11366v1.

Vertices are numbered 0, ..., vertex_count - 1. Graphs are represented by a
flattened symmetric adjacency matrix produced by `empty_adjacency()` and
`add_edge()`. The algorithms enumerate independent sets and truncated
multivariate series, so they are intentionally bounded reference routines—not
large-scale inference software.
"""

from std.math import abs


comptime MAX_VERTICES = 16
comptime MAX_COEFFICIENT_STATES = 2_000_000


def empty_adjacency(vertex_count: Int) raises -> List[Int]:
    """Create a flattened adjacency matrix for an edgeless graph."""
    if vertex_count <= 0 or vertex_count > MAX_VERTICES:
        raise Error("vertex_count must be in 1..16")
    var adjacency = List[Int]()
    for _ in range(vertex_count * vertex_count):
        adjacency.append(0)
    return adjacency^


def add_edge(
    mut adjacency: List[Int], vertex_count: Int, u: Int, v: Int
) raises:
    """Add an undirected edge to a flattened adjacency matrix."""
    _validate_matrix_shape(adjacency, vertex_count)
    if u < 0 or u >= vertex_count or v < 0 or v >= vertex_count or u == v:
        raise Error("edge endpoints must be distinct valid vertex indices")
    adjacency[u * vertex_count + v] = 1
    adjacency[v * vertex_count + u] = 1


def is_chordal(adjacency: List[Int], vertex_count: Int) raises -> Bool:
    """Check decomposability with maximum-cardinality search."""
    _validate_graph(adjacency, vertex_count)
    var weight = List[Int]()
    var number = List[Int]()
    for _ in range(vertex_count):
        weight.append(0)
        number.append(-1)

    for selected_count in range(vertex_count):
        var selected = -1
        var best_weight = -1
        for v in range(vertex_count):
            if number[v] == -1 and weight[v] > best_weight:
                selected = v
                best_weight = weight[v]
        number[selected] = vertex_count - selected_count
        for u in range(vertex_count):
            if number[u] == -1 and adjacency[selected * vertex_count + u] == 1:
                weight[u] += 1

    for v in range(vertex_count):
        for u in range(vertex_count):
            if adjacency[v * vertex_count + u] == 1 and number[u] > number[v]:
                for w in range(u + 1, vertex_count):
                    if (
                        adjacency[v * vertex_count + w] == 1
                        and number[w] > number[v]
                        and adjacency[u * vertex_count + w] == 0
                    ):
                        return False
    return True


def independence_polynomial(
    adjacency: List[Int], vertex_count: Int, y: List[Float64]
) raises -> Float64:
    r"""Evaluate δ_G(y), summing over every independent set, including ∅."""
    _validate_graph(adjacency, vertex_count)
    _validate_vector_length(y, vertex_count, "y")
    var total = 0.0
    for mask in range(1 << vertex_count):
        if _is_independent(adjacency, vertex_count, mask):
            var term = 1.0
            for v in range(vertex_count):
                if _contains(mask, v):
                    term *= y[v]
            total += term
    return total


def signed_graph_polynomial(
    adjacency: List[Int], vertex_count: Int, x: List[Float64]
) raises -> Float64:
    r"""Evaluate Δ_G(x) = δ_G(-x)."""
    _validate_graph(adjacency, vertex_count)
    _validate_vector_length(x, vertex_count, "x")
    return _signed_graph_polynomial_on(
        adjacency, vertex_count, x, (1 << vertex_count) - 1
    )


def first_graph_multinomial_coefficient(
    adjacency: List[Int], vertex_count: Int, r: Int, n: List[Int]
) raises -> Float64:
    r"""Return [y^n] δ_G(y)^r (Definition 2.2, first type)."""
    _validate_graph(adjacency, vertex_count)
    _validate_count_vector(n, vertex_count, "n")
    if r <= 0:
        raise Error("r must be a positive integer")
    var state_count = _coefficient_state_count(n)
    var current = _zeros(state_count)
    current[0] = 1.0

    for _ in range(r):
        var next_values = _zeros(state_count)
        for state in range(state_count):
            if current[state] == 0.0:
                continue
            for mask in range(1 << vertex_count):
                if not _is_independent(adjacency, vertex_count, mask):
                    continue
                var next_state = _shift_state(state, n, mask)
                if next_state >= 0:
                    next_values[next_state] += current[state]
        current = next_values^
    return current[state_count - 1]


def second_graph_multinomial_coefficient(
    adjacency: List[Int], vertex_count: Int, s: Float64, n: List[Int]
) raises -> Float64:
    r"""Return [x^n] Δ_G(x)^(-s) (Definition 2.2, second type)."""
    _validate_graph(adjacency, vertex_count)
    _validate_count_vector(n, vertex_count, "n")
    if s <= 0.0:
        raise Error("s must be positive")

    var state_count = _coefficient_state_count(n)
    var target_degree = _sum_int(n)
    var q_power = _zeros(state_count)
    q_power[0] = 1.0
    var result = q_power[state_count - 1]
    var generalized_binomial = 1.0

    for m in range(1, target_degree + 1):
        var next_values = _zeros(state_count)
        for state in range(state_count):
            if q_power[state] == 0.0:
                continue
            for mask in range(1, 1 << vertex_count):
                if not _is_independent(adjacency, vertex_count, mask):
                    continue
                var next_state = _shift_state(state, n, mask)
                if next_state >= 0:
                    var sign = -1.0 if _popcount(mask) % 2 == 1 else 1.0
                    next_values[next_state] += q_power[state] * sign
        q_power = next_values^
        generalized_binomial *= -(s + Float64(m - 1)) / Float64(m)
        result += generalized_binomial * q_power[state_count - 1]
    return result


def graphical_multinomial_pmf(
    adjacency: List[Int],
    vertex_count: Int,
    r: Int,
    y: List[Float64],
    n: List[Int],
) raises -> Float64:
    r"""Evaluate mult_G(r, y) from Definition 3.1(1)."""
    _require_decomposable(adjacency, vertex_count)
    _validate_positive_vector(y, vertex_count, "y")
    var coefficient = first_graph_multinomial_coefficient(
        adjacency, vertex_count, r, n
    )
    if coefficient == 0.0:
        return 0.0
    var weight = _monomial(y, n)
    var normalizer = independence_polynomial(
        adjacency, vertex_count, y
    ) ** Float64(r)
    return coefficient * weight / normalizer


def graphical_negative_multinomial_pmf(
    adjacency: List[Int],
    vertex_count: Int,
    r: Float64,
    x: List[Float64],
    n: List[Int],
) raises -> Float64:
    r"""Evaluate nm_G(r, x) from Definition 3.1(2)."""
    _require_decomposable(adjacency, vertex_count)
    if r <= 0.0:
        raise Error("r must be positive")
    _validate_positive_vector(x, vertex_count, "x")
    _require_negative_multinomial_domain(adjacency, vertex_count, x)
    var coefficient = second_graph_multinomial_coefficient(
        adjacency, vertex_count, r, n
    )
    return (
        coefficient
        * _monomial(x, n)
        * signed_graph_polynomial(adjacency, vertex_count, x) ** r
    )


def graphical_hypergeometric_pmf(
    adjacency: List[Int],
    vertex_count: Int,
    M: Int,
    K: List[Int],
    r: Int,
    n: List[Int],
) raises -> Float64:
    r"""Evaluate hg_G(M, K, r) from Definition 4.1."""
    _require_decomposable(adjacency, vertex_count)
    _validate_count_vector(K, vertex_count, "K")
    _validate_count_vector(n, vertex_count, "n")
    if M <= 0 or r <= 0 or M <= r:
        raise Error("hg_G requires positive integers M, r with M > r")
    if (
        first_graph_multinomial_coefficient(adjacency, vertex_count, M, K)
        == 0.0
    ):
        raise Error("K is outside N_{G,M}")
    var remainder = _subtract_counts(K, n)
    if len(remainder) == 0:
        return 0.0
    var denominator = first_graph_multinomial_coefficient(
        adjacency, vertex_count, M, K
    )
    var left = first_graph_multinomial_coefficient(
        adjacency, vertex_count, r, n
    )
    var right = first_graph_multinomial_coefficient(
        adjacency, vertex_count, M - r, remainder
    )
    return left * right / denominator


def graphical_negative_hypergeometric_pmf(
    adjacency: List[Int],
    vertex_count: Int,
    M: Int,
    K: List[Int],
    r: Float64,
    n: List[Int],
) raises -> Float64:
    r"""Evaluate nhg_G(M, K, r) from Definition 4.5."""
    _require_decomposable(adjacency, vertex_count)
    _validate_count_vector(K, vertex_count, "K")
    _validate_count_vector(n, vertex_count, "n")
    if M <= 0 or r <= 0.0 or r > Float64(M):
        raise Error("nhg_G requires M > 0 and 0 < r <= M")
    var remainder = _subtract_counts(K, n)
    if len(remainder) == 0:
        return 0.0
    var left = second_graph_multinomial_coefficient(
        adjacency, vertex_count, r, n
    )
    var right = second_graph_multinomial_coefficient(
        adjacency, vertex_count, Float64(M) - r + 1.0, remainder
    )
    var denominator = second_graph_multinomial_coefficient(
        adjacency, vertex_count, Float64(M) + 1.0, K
    )
    return left * right / denominator


def _validate_matrix_shape(adjacency: List[Int], vertex_count: Int) raises:
    if vertex_count <= 0 or vertex_count > MAX_VERTICES:
        raise Error("vertex_count must be in 1..16")
    if len(adjacency) != vertex_count * vertex_count:
        raise Error("adjacency must contain vertex_count^2 entries")


def _validate_graph(adjacency: List[Int], vertex_count: Int) raises:
    _validate_matrix_shape(adjacency, vertex_count)
    for u in range(vertex_count):
        if adjacency[u * vertex_count + u] != 0:
            raise Error("self-loops are not supported")
        for v in range(vertex_count):
            var uv = adjacency[u * vertex_count + v]
            if (uv != 0 and uv != 1) or uv != adjacency[v * vertex_count + u]:
                raise Error("adjacency must be symmetric and binary")


def _require_decomposable(adjacency: List[Int], vertex_count: Int) raises:
    if not is_chordal(adjacency, vertex_count):
        raise Error(
            "the four general PMFs require a decomposable (chordal) graph"
        )


def _validate_vector_length[
    T: Movable
](values: List[T], vertex_count: Int, name: String) raises:
    if len(values) != vertex_count:
        raise Error(name + " must have one value per vertex")


def _validate_count_vector(
    values: List[Int], vertex_count: Int, name: String
) raises:
    _validate_vector_length(values, vertex_count, name)
    for value in values:
        if value < 0:
            raise Error(name + " must contain nonnegative counts")


def _validate_positive_vector(
    values: List[Float64], vertex_count: Int, name: String
) raises:
    _validate_vector_length(values, vertex_count, name)
    for value in values:
        if value <= 0.0:
            raise Error(name + " must contain strictly positive values")


def _require_negative_multinomial_domain(
    adjacency: List[Int], vertex_count: Int, x: List[Float64]
) raises:
    for induced_mask in range(1 << vertex_count):
        if (
            _signed_graph_polynomial_on(
                adjacency, vertex_count, x, induced_mask
            )
            <= 0.0
        ):
            raise Error(
                "x is outside M_G: an induced signed polynomial is nonpositive"
            )


def _signed_graph_polynomial_on(
    adjacency: List[Int], vertex_count: Int, x: List[Float64], induced_mask: Int
) -> Float64:
    var total = 0.0
    for mask in range(1 << vertex_count):
        if (mask & induced_mask) != mask:
            continue
        if _is_independent(adjacency, vertex_count, mask):
            var term = 1.0
            for v in range(vertex_count):
                if _contains(mask, v):
                    term *= x[v]
            total += -term if _popcount(mask) % 2 == 1 else term
    return total


def _is_independent(adjacency: List[Int], vertex_count: Int, mask: Int) -> Bool:
    for u in range(vertex_count):
        if not _contains(mask, u):
            continue
        for v in range(u + 1, vertex_count):
            if _contains(mask, v) and adjacency[u * vertex_count + v] == 1:
                return False
    return True


def _contains(mask: Int, vertex: Int) -> Bool:
    return (mask & (1 << vertex)) != 0


def _popcount(mask: Int) -> Int:
    var value = mask
    var count = 0
    while value != 0:
        count += value & 1
        value = value >> 1
    return count


def _zeros(size: Int) -> List[Float64]:
    var values = List[Float64]()
    for _ in range(size):
        values.append(0.0)
    return values^


def _coefficient_state_count(target: List[Int]) raises -> Int:
    var count = 1
    for value in target:
        if value + 1 > MAX_COEFFICIENT_STATES // count:
            raise Error("coefficient grid exceeds 2,000,000 states")
        count *= value + 1
    return count


def _shift_state(state: Int, target: List[Int], mask: Int) -> Int:
    var stride = 1
    var next_state = state
    for v in range(len(target)):
        if _contains(mask, v):
            var exponent = (state // stride) % (target[v] + 1)
            if exponent == target[v]:
                return -1
            next_state += stride
        stride *= target[v] + 1
    return next_state


def _sum_int(values: List[Int]) -> Int:
    var total = 0
    for value in values:
        total += value
    return total


def _monomial(base: List[Float64], exponent: List[Int]) -> Float64:
    var result = 1.0
    for i in range(len(base)):
        for _ in range(exponent[i]):
            result *= base[i]
    return result


def _subtract_counts(total: List[Int], part: List[Int]) -> List[Int]:
    if len(total) != len(part):
        return List[Int]()
    var result = List[Int]()
    for i in range(len(total)):
        if part[i] > total[i]:
            return List[Int]()
        result.append(total[i] - part[i])
    return result^
