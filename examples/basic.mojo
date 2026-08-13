from graphical_counts import (
    add_edge,
    empty_adjacency,
    graphical_hypergeometric_pmf,
    graphical_multinomial_pmf,
)


def main() raises:
    # The complete graph recovers the classical multivariate distributions.
    var adjacency = empty_adjacency(3)
    add_edge(adjacency, 3, 0, 1)
    add_edge(adjacency, 3, 0, 2)
    add_edge(adjacency, 3, 1, 2)

    var y: List[Float64] = [0.5, 1.0, 1.5]
    var counts: List[Int] = [1, 1, 0]
    var mult_probability = graphical_multinomial_pmf(adjacency, 3, 2, y, counts)

    var population_counts: List[Int] = [2, 1, 1]
    var sample_counts: List[Int] = [1, 0, 0]
    var hg_probability = graphical_hypergeometric_pmf(
        adjacency, 3, 6, population_counts, 2, sample_counts
    )

    print("mult_G probability:", mult_probability)
    print("hg_G probability:", hg_probability)
