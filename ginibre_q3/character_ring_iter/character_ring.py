"""Klimyk's formula for iterated tensor with the adjoint representation.

Given V = sum_lambda c_lambda V_lambda as a dict {lambda_fw: c_lambda} in
fundamental-weight coordinates, compute V' = V ⊗ adj using

  V_lambda ⊗ V_theta = sum_{eta in weights(V_theta)} m_theta(eta) · [|lambda + eta|]_dom

where [|sigma|]_dom is the Racah–Speiser dominant reflection with rho-shift:

  sigma -> sigma + rho.  Reflect to the strict-dominant chamber (all coords > 0)
  tracking sign.  If any coord is 0 at any point, the contribution is 0.
  Finally subtract rho.

For computing adjoint tensor moments m_k = mult(triv, adj^⊗k), we start from
V = triv (V_0, coords all zero) and iterate k times.  The multiplicity of
the trivial irreducible at each step is V[(0,...,0)].
"""

from collections import defaultdict
from lie_data import CARTAN_MATRICES, adjoint_weights_fw, rank


class LieGroup:
    def __init__(self, name):
        self.name = name
        self.cartan = CARTAN_MATRICES[name]
        self.rank = len(self.cartan)
        self.adj_weights = adjoint_weights_fw(name)
        self._dominant_cache = {}

    def dominant_reflect(self, mu_fw):
        """Apply Weyl dot-action reflections to bring mu + rho to strict dominant.

        mu_fw: tuple of ints, fundamental-weight coordinates.
        Returns (dom_mu_fw, sign) where dom_mu_fw is dominant and sign ∈ {+1, -1}.
        Returns None if mu + rho hits a wall (any coord = 0 at any reflection point).
        """
        if mu_fw in self._dominant_cache:
            return self._dominant_cache[mu_fw]

        r = self.rank
        # Shift by rho (= (1,1,...,1) in fw basis).
        n = [x + 1 for x in mu_fw]
        sign = 1
        # Bounded iteration (Weyl group is finite; guard loop)
        for _step in range(10000):
            # Find coord with n[i] <= 0
            bad = -1
            for i in range(r):
                if n[i] < 0:
                    bad = i
                    break
                if n[i] == 0:
                    self._dominant_cache[mu_fw] = None
                    return None  # wall
            if bad < 0:
                break  # all n[i] > 0, done
            # Reflect by s_{alpha_bad}.  alpha_i in FW basis has k-th coord
            # C[i][k], so s_{alpha_i}(mu)_k = n_k - n_i * C[i][k].
            # For k=bad: n[bad] - n[bad]*2 = -n[bad].
            ni = n[bad]
            for k in range(r):
                n[k] -= ni * self.cartan[bad][k]
            sign = -sign
        else:
            raise RuntimeError(f'dominant_reflect did not terminate on {mu_fw}')
        # Subtract rho back
        dom = tuple(x - 1 for x in n)
        res = (dom, sign)
        self._dominant_cache[mu_fw] = res
        return res

    def tensor_with_adj(self, V):
        """V: dict {weight_tuple: int_multiplicity}.  Return V ⊗ adj as dict."""
        W = defaultdict(int)
        for mu, c_mu in V.items():
            if c_mu == 0:
                continue
            for eta, m_eta in self.adj_weights:
                # target = mu + eta in fundamental weight basis
                target = tuple(mu[i] + eta[i] for i in range(self.rank))
                res = self.dominant_reflect(target)
                if res is None:
                    continue
                dom, sign = res
                W[dom] += sign * c_mu * m_eta
        # Prune zeros
        return {k: v for k, v in W.items() if v != 0}


def compute_moments(group_name, max_k, verbose=False):
    """Compute m_k = mult(triv, adj^⊗k) for k = 0, 1, ..., max_k."""
    G = LieGroup(group_name)
    zero = tuple([0] * G.rank)
    V = {zero: 1}  # V_triv
    moments = [V.get(zero, 0)]  # m_0 = 1
    for k in range(1, max_k + 1):
        V = G.tensor_with_adj(V)
        m_k = V.get(zero, 0)
        moments.append(m_k)
        if verbose:
            print(f'  m_{k} = {m_k}  (|support| = {len(V)})', flush=True)
    return moments


if __name__ == '__main__':
    import sys
    group = sys.argv[1] if len(sys.argv) > 1 else 'G2'
    max_k = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    print(f'Computing adjoint moments for {group}, k = 0..{max_k}:')
    moments = compute_moments(group, max_k, verbose=True)
    print()
    print(f'Full moment list: {moments}')
