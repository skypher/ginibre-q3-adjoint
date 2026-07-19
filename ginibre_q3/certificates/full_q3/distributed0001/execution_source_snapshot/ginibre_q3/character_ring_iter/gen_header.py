#!/usr/bin/env python3
"""Generate lie_data.h from lie_data.py."""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lie_data import CARTAN_MATRICES, adjoint_weights_fw

lines = [
    '// Auto-generated from lie_data.py',
    '#pragma once',
    '#include <array>',
    '#include <vector>',
    '',
]
# Keep this list synchronized with the namespaces consumed by the checked-in
# header.  The classical verifier constructs B/C/D data algorithmically and
# does not use generated lie_B2, lie_B3, or lie_D4 namespaces.
for g in ['G2', 'F4', 'E6', 'E7', 'E8', 'A1', 'A2']:
    cartan = CARTAN_MATRICES[g]
    r = len(cartan)
    adj = adjoint_weights_fw(g)
    lines.append(f'namespace lie_{g} {{')
    lines.append(f'    constexpr int rank = {r};')
    lines.append(f'    constexpr int cartan[{r}][{r}] = {{')
    for row in cartan:
        lines.append('        {' + ', '.join(str(x) for x in row) + '},')
    lines.append('    };')
    lines.append(f'    constexpr int num_adj_weights = {len(adj)};')
    lines.append(f'    constexpr int adj_weight_coord[{len(adj)}][{r}] = {{')
    for w, m in adj:
        lines.append('        {' + ', '.join(str(x) for x in w) + '},')
    lines.append('    };')
    lines.append(f'    constexpr int adj_weight_mult[{len(adj)}] = {{')
    lines.append('        ' + ', '.join(str(m) for w, m in adj))
    lines.append('    };')
    lines.append('}')
    lines.append('')

dst = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'lie_data.h')
with open(dst, 'w') as f:
    f.write('\n'.join(lines))
print(f'wrote {dst}')
