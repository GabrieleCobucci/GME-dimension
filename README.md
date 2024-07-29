# Convex programming for GME-dimension

Detecting, quantifying and characterising relevant forms of entanglement is of central interest in quantum information science. In recent years, high-dimensional systems have emerged as a major frontier for entanglement. In bipartite systems, entanglement has been detected far beyond qubit dimensions but the high-dimensional multipartite regime is more challenging. We propose to quantify the genuine multipartite entanglement dimension via the worst-case Schmidt number needed to generate the state by classically mixing entangled states over different cuts of the system.
Here we detect this form of entanglement by convex programming methods as proposed in [1].
In particular, in this repository we provide the functions needed to evaluate the linear programs when the pure target state is:

n-partite high-dimensional GHZ state;
n-partite high-dimensional linear cluster state.


# List of the functions

## Cluster_state-LP

LPGMENmixercl.m
This function evaluates the LP (YALMIP) with a pure linear cluster state and finds the maximum visibility to have a specific GME-dimension

kpositivecluster.m
This function applies the k-positive generalised reduction map on a density matrix diagonal in the graph basis

Stirling2nd.m [2]

SetPartition.m [2]

## GHZ_state-LP

LPGMENmixerghz.m
This function evaluates the LP (YALMIP) with a pure n-partite high-dimensional GHZ state and finds the maximum visibility to have a specific GME-dimension

kpositiveghz.m
This function applies the k-positive generalised reduction map on a density matrix diagonal in the GHZ basis

Stirling2dn.m [2]

SetPartition.m [2]

# References
[1]: Gabriele Cobucci, Armin Tavakoli (2024). Characterising and detecting genuinely high-dimensional genuine multipartite entanglement (https://arxiv.org/abs/2402.06234)
[2]: Bruno Luong (2024). Set partition (https://www.mathworks.com/matlabcentral/fileexchange/24133-set-partition), MATLAB Central File Exchange. Retrieved February 5, 2024.
