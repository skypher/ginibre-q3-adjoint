"""Root system data for exceptional compact simple Lie groups.

Cartan matrix convention: A[j][i] = <alpha_j, alpha_i^vee>, i.e., rows are coroots
and columns are roots, so that alpha_i = sum_k A[k][i] omega_k (i-th column).
Bourbaki numbering throughout.

Positive roots are generated algorithmically from the Cartan matrix via
height-stratified enumeration using the alpha-string structure:
  beta + alpha_i is a root  iff  p - <beta, alpha_i^vee> > 0
where p is the max integer with beta - p*alpha_i still a root.
"""

CARTAN_MATRICES = {
    'G2': [
        [2, -1],
        [-3, 2],
    ],
    'F4': [
        [2, -1, 0, 0],
        [-1, 2, -2, 0],
        [0, -1, 2, -1],
        [0, 0, -1, 2],
    ],
    'E6': [
        [2, 0, -1, 0, 0, 0],
        [0, 2, 0, -1, 0, 0],
        [-1, 0, 2, -1, 0, 0],
        [0, -1, -1, 2, -1, 0],
        [0, 0, 0, -1, 2, -1],
        [0, 0, 0, 0, -1, 2],
    ],
    'E7': [
        [2, 0, -1, 0, 0, 0, 0],
        [0, 2, 0, -1, 0, 0, 0],
        [-1, 0, 2, -1, 0, 0, 0],
        [0, -1, -1, 2, -1, 0, 0],
        [0, 0, 0, -1, 2, -1, 0],
        [0, 0, 0, 0, -1, 2, -1],
        [0, 0, 0, 0, 0, -1, 2],
    ],
    'E8': [
        [2, 0, -1, 0, 0, 0, 0, 0],
        [0, 2, 0, -1, 0, 0, 0, 0],
        [-1, 0, 2, -1, 0, 0, 0, 0],
        [0, -1, -1, 2, -1, 0, 0, 0],
        [0, 0, 0, -1, 2, -1, 0, 0],
        [0, 0, 0, 0, -1, 2, -1, 0],
        [0, 0, 0, 0, 0, -1, 2, -1],
        [0, 0, 0, 0, 0, 0, -1, 2],
    ],
    # Classical families for testing: A_1, B_2, etc.
    'A1': [[2]],
    'A2': [[2, -1], [-1, 2]],
    'B2': [[2, -2], [-1, 2]],
    'B3': [[2, -1, 0], [-1, 2, -2], [0, -1, 2]],
    'D4': [[2, 0, -1, 0], [0, 2, -1, 0], [-1, -1, 2, -1], [0, 0, -1, 2]],
}


def positive_roots_in_simple_basis(cartan):
    """Enumerate positive roots (as tuples of simple-root coefficients).

    Uses the alpha-string trick: beta + alpha_i is a root iff there is room
    in the string.  Return positive roots sorted by height.
    """
    r = len(cartan)
    # Start: simple roots = standard basis vectors.
    simple = [tuple(1 if j == i else 0 for j in range(r)) for i in range(r)]
    roots = set(simple)
    # height-stratified BFS
    frontier = list(simple)
    while frontier:
        new_frontier = []
        for beta in frontier:
            for i in range(r):
                # <beta, alpha_i^vee> = sum_j beta_j A[j][i]
                pairing = sum(beta[j] * cartan[j][i] for j in range(r))
                # max p with beta - p*alpha_i still a root
                p = 0
                while True:
                    cand = tuple(beta[j] - (p + 1) * (1 if j == i else 0) for j in range(r))
                    if cand in roots:
                        p += 1
                    else:
                        break
                # room upward: q = p - pairing; if > 0, beta + alpha_i is a root
                q_max = p - pairing
                if q_max > 0:
                    cand_up = tuple(beta[j] + (1 if j == i else 0) for j in range(r))
                    if cand_up not in roots:
                        roots.add(cand_up)
                        new_frontier.append(cand_up)
        frontier = new_frontier
    return sorted(roots, key=lambda b: sum(b))


def simple_to_fundamental(beta, cartan):
    """Convert root in simple-root basis to fundamental weight basis.

    omega_k is defined by <omega_k, alpha_j^vee> = delta_{kj}.  So if
    alpha_i = sum_k b_{ki} omega_k, then b_{ji} = <alpha_i, alpha_j^vee> =
    cartan[i][j].  Hence b_{ki} = cartan[i][k], and
      beta = sum_i beta_i alpha_i = sum_k (sum_i beta_i * cartan[i][k]) omega_k.
    """
    r = len(cartan)
    return tuple(sum(beta[i] * cartan[i][k] for i in range(r)) for k in range(r))


def adjoint_weights_fw(group):
    """Return list of (weight_in_fw, multiplicity) for the adjoint representation."""
    cartan = CARTAN_MATRICES[group]
    r = len(cartan)
    pos = positive_roots_in_simple_basis(cartan)
    weights = []
    for beta in pos:
        fw = simple_to_fundamental(beta, cartan)
        weights.append((fw, 1))
        weights.append((tuple(-x for x in fw), 1))
    # zero weight with multiplicity = rank
    weights.append((tuple([0] * r), r))
    return weights


def highest_root_fw(group):
    """Highest root in fundamental weight basis (= highest weight of adj)."""
    cartan = CARTAN_MATRICES[group]
    pos = positive_roots_in_simple_basis(cartan)
    return simple_to_fundamental(pos[-1], cartan)


def rank(group):
    return len(CARTAN_MATRICES[group])


def dim_adj(group):
    cartan = CARTAN_MATRICES[group]
    return len(cartan) + 2 * len(positive_roots_in_simple_basis(cartan))


if __name__ == '__main__':
    for g in ['A1', 'G2', 'F4', 'E6', 'E7', 'E8']:
        pr = positive_roots_in_simple_basis(CARTAN_MATRICES[g])
        print(f'{g}: rank={rank(g)}, |R+|={len(pr)}, dim(adj)={dim_adj(g)}')
        print(f'  highest root (fw basis): {highest_root_fw(g)}')
