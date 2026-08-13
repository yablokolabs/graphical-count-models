from std.math import abs
from std.testing import TestSuite, assert_equal, assert_false, assert_true

from graphical_counts import (
    add_edge,
    empty_adjacency,
    first_graph_multinomial_coefficient,
    graphical_hypergeometric_pmf,
    graphical_multinomial_pmf,
    graphical_negative_hypergeometric_pmf,
    graphical_negative_multinomial_pmf,
    independence_polynomial,
    is_chordal,
    second_graph_multinomial_coefficient,
    signed_graph_polynomial,
)


def assert_close(
    actual: Float64, expected: Float64, tolerance: Float64 = 1e-12
) raises:
    assert_true(abs(actual - expected) <= tolerance)


def complete_graph(vertex_count: Int) raises -> List[Int]:
    var adjacency = empty_adjacency(vertex_count)
    for u in range(vertex_count):
        for v in range(u + 1, vertex_count):
            add_edge(adjacency, vertex_count, u, v)
    return adjacency^


def test_graph_polynomials_include_all_independent_sets() raises:
    # Path 0--1--2 has independent sets ∅, {0}, {1}, {2}, {0,2}.
    var adjacency = empty_adjacency(3)
    add_edge(adjacency, 3, 0, 1)
    add_edge(adjacency, 3, 1, 2)
    var y: List[Float64] = [1.0, 2.0, 3.0]
    var x: List[Float64] = [0.1, 0.2, 0.3]
    assert_close(independence_polynomial(adjacency, 3, y), 10.0)
    assert_close(signed_graph_polynomial(adjacency, 3, x), 0.43)


def test_chordality_check() raises:
    var path = empty_adjacency(4)
    add_edge(path, 4, 0, 1)
    add_edge(path, 4, 1, 2)
    add_edge(path, 4, 2, 3)
    assert_true(is_chordal(path, 4))

    var cycle = path.copy()
    add_edge(cycle, 4, 3, 0)
    assert_false(is_chordal(cycle, 4))


def test_coefficient_families_match_classical_limits() raises:
    var complete = complete_graph(2)
    var empty = empty_adjacency(2)
    var n: List[Int] = [1, 1]

    # Complete G: δ_G = 1 + y0 + y1 and Δ_G = 1 - x0 - x1.
    assert_close(first_graph_multinomial_coefficient(complete, 2, 2, n), 2.0)
    assert_close(second_graph_multinomial_coefficient(complete, 2, 1.0, n), 2.0)

    # Empty G: δ_G = (1+y0)(1+y1), Δ_G = (1-x0)(1-x1).
    assert_close(first_graph_multinomial_coefficient(empty, 2, 2, n), 4.0)
    assert_close(second_graph_multinomial_coefficient(empty, 2, 1.0, n), 1.0)


def test_multinomial_and_negative_multinomial_pmf() raises:
    var complete = complete_graph(2)
    var counts: List[Int] = [1, 1]
    var y: List[Float64] = [1.0, 1.0]
    var x: List[Float64] = [0.2, 0.3]

    assert_close(
        graphical_multinomial_pmf(complete, 2, 2, y, counts), 2.0 / 9.0
    )
    assert_close(
        graphical_negative_multinomial_pmf(complete, 2, 1.0, x, counts),
        0.06,
    )


def test_graphical_hypergeometric_normalizes() raises:
    var complete = complete_graph(2)
    var K: List[Int] = [2, 1]
    var total = 0.0
    for n0 in range(3):
        for n1 in range(2):
            var n: List[Int] = [n0, n1]
            total += graphical_hypergeometric_pmf(complete, 2, 5, K, 2, n)
    assert_close(total, 1.0)
    var target: List[Int] = [1, 0]
    assert_close(
        graphical_hypergeometric_pmf(complete, 2, 5, K, 2, target), 0.4
    )


def test_graphical_negative_hypergeometric_normalizes() raises:
    var complete = complete_graph(2)
    var K: List[Int] = [2, 1]
    var total = 0.0
    for n0 in range(3):
        for n1 in range(2):
            var n: List[Int] = [n0, n1]
            total += graphical_negative_hypergeometric_pmf(
                complete, 2, 4, K, 2.0, n
            )
    assert_close(total, 1.0)
    var target: List[Int] = [1, 0]
    assert_close(
        graphical_negative_hypergeometric_pmf(complete, 2, 4, K, 2.0, target),
        8.0 / 35.0,
    )


def main() raises:
    TestSuite.discover_tests[__functions_in_module()]().run()
