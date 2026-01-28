# SizeStructuredPopulations

Code for the paper [Size-structured populations with growth
fluctuations: Feynman–Kac formula and decoupling](https://arxiv.org/abs/2508.14680). The figures themselves are generated in the notebooks in paper-figures. example-notebooks contains some additional examples, showing how to setup the models.  

To speed up population simulations, `grow_tree!` can spawn subtrees in parallel by setting `parallel_depth` (e.g. `grow_tree!(cell, model, terminate; parallel_depth=2)`).
