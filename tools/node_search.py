import numpy as np
import re
from collections import defaultdict

# --------------------- USER-ADJUSTABLE PARAMETERS ---------------------
PSI_THRESHOLD_RATIO = 0.05   # Minimum |psi| fraction for a valid node
Q_TOLERANCE_POINTS = 2       # Node may wander across this many grid steps
# ----------------------------------------------------------------------

def find_zero_crossings(x, y, threshold):
    """
    Identify zero crossings in y(x) where sign changes and both adjacent values
    are above a given threshold. Returns list of interpolated x values.
    """
    zero_crossings = []
    for i in range(1, len(y)):
        if y[i-1] * y[i] < 0:
            if abs(y[i-1]) >= threshold and abs(y[i]) >= threshold:
                # Linear interpolation to estimate zero crossing
                x0, x1 = x[i-1], x[i]
                y0, y1 = y[i-1], y[i]
                x_zero = x0 - y0 * (x1 - x0) / (y1 - y0)
                zero_crossings.append(x_zero)
    return zero_crossings

def analyze_wavefunction_nodes(filepath):
    """
    Read a FEMvib-format wavefunction file and return node locations along
    each vibrational coordinate, accounting for digitization noise.
    """
    # Load and parse the file
    data = []
    with open(filepath, 'r') as f:
        for line in f:
            if re.match(r'^\s*#', line):
                continue  # Skip comment lines
            tokens = line.strip().split()
            if len(tokens) < 3:
                continue
            data.append([float(tok) for tok in tokens])

    data = np.array(data)
    coords = data[:, :-1]
    psi = data[:, -1]
    ndim = coords.shape[1]

    # Get coordinate axes and reshape the ψ array
    axes = [np.unique(coords[:, i]) for i in range(ndim)]
    grid_shape = tuple(len(ax) for ax in axes)
    psi_grid = psi.reshape(grid_shape)

    psi_max = np.max(np.abs(psi_grid))
    threshold = PSI_THRESHOLD_RATIO * psi_max
    node_locations = defaultdict(list)

    # Loop over each axis and scan 1D slices
    for axis in range(ndim):
        other_axes = [i for i in range(ndim) if i != axis]
        other_grids = [axes[i] for i in other_axes]

        # Loop over the indices of other axes
        it = np.ndindex(*[len(g) for g in other_grids])
        for idx in it:
            slicer = [slice(None) if i == axis else idx[other_axes.index(i)]
                      for i in range(ndim)]
            psi_line = psi_grid[tuple(slicer)]
            q_line = axes[axis]

            zeros = find_zero_crossings(q_line, psi_line, threshold)
            for z in zeros:
                coord_vals = [axes[i][idx[other_axes.index(i)]]
                              if i != axis else z for i in range(ndim)]
                node_locations[axis].append(tuple(coord_vals))

    return node_locations

# ------------------ Example execution (uncomment to use) ------------------
 if __name__ == "__main__":
     filename = "xx.txt"
     nodes = analyze_wavefunction_nodes(filename)
     for axis, node_list in nodes.items():
         print(f"\nNodes along Q{axis+1}:")
         for pt in node_list:
             print("  ", pt)

