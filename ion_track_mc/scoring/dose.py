"""Dose scoring: energy deposited per unit mass, D = E / m  [Gy = J/kg].

Two geometries (scoring/README_scoring.md §2):
  - `VoxelDoseScorer`: a box grid, for depth-dose (Bragg) curves and 2D dose
    maps;
  - `RadialDoseScorer`: cylindrical shells around the track axis, for the
    radial dose profile D(r).

Each scorer adds every deposit to the bin it falls in; dose is computed at
the end by dividing each bin's energy by the mass of water it contains.
"""
import numpy as np

from ion_track_mc import constants
from ion_track_mc.physics.medium import Medium
from ion_track_mc.scoring.scorer import Scorer


class VoxelDoseScorer(Scorer):
    """Cartesian voxel grid dose scorer.

    Bin edges are given explicitly per axis (`x_edges`, `y_edges`,
    `z_edges`), in nm -- so voxel shape/aspect ratio is fully under the
    caller's control (e.g. one wide transverse bin + many thin z bins for a
    depth-dose/Bragg curve, see `make_bragg_curve_scorer`).

    Deposits outside every edge range are silently dropped -- this scorer
    only accounts for energy landing inside its grid. Size the grid
    generously (or check `total_deposited_eV()` against the run's known
    total loss) if completeness matters.
    """

    def __init__(self, medium: Medium, x_edges: np.ndarray, y_edges: np.ndarray, z_edges: np.ndarray):
        self.medium = medium
        self.x_edges = np.asarray(x_edges, dtype=float)
        self.y_edges = np.asarray(y_edges, dtype=float)
        self.z_edges = np.asarray(z_edges, dtype=float)
        shape = (len(self.x_edges) - 1, len(self.y_edges) - 1, len(self.z_edges) - 1)
        self.energy_grid_eV = np.zeros(shape)  # [eV] accumulated per voxel

    def score(self, position: np.ndarray, energy: float) -> None:
        # Index of the bin each coordinate falls in (-1 or len if outside).
        ix = int(np.searchsorted(self.x_edges, position[0], side="right") - 1)
        iy = int(np.searchsorted(self.y_edges, position[1], side="right") - 1)
        iz = int(np.searchsorted(self.z_edges, position[2], side="right") - 1)
        shape = self.energy_grid_eV.shape
        if 0 <= ix < shape[0] and 0 <= iy < shape[1] and 0 <= iz < shape[2]:
            self.energy_grid_eV[ix, iy, iz] += energy

    def voxel_volumes_nm3(self) -> np.ndarray:
        dx = np.diff(self.x_edges)
        dy = np.diff(self.y_edges)
        dz = np.diff(self.z_edges)
        return dx[:, None, None] * dy[None, :, None] * dz[None, None, :]

    def dose_Gy(self) -> np.ndarray:
        """Dose [Gy] per voxel: energy[J] / mass[kg], mass from `medium.density_g_cm3`."""
        energy_J = self.energy_grid_eV * constants.EV_TO_JOULE
        volume_cm3 = self.voxel_volumes_nm3() / constants.CM3_TO_NM3
        mass_kg = self.medium.density_g_cm3 * volume_cm3 * 1.0e-3   # g -> kg
        return energy_J / mass_kg

    def total_deposited_eV(self) -> float:
        return float(self.energy_grid_eV.sum())


class RadialDoseScorer(Scorer):
    """Concentric cylindrical-shell dose scorer around a track axis.

    Bins energy deposits by perpendicular distance from an infinite line
    (`axis_point` + `axis_direction`) -- the natural geometry for a track's
    radial dose profile D(r): how dose falls off with distance from the
    path of the ion.

    `length_nm` is the axial extent the scored deposits are expected to
    span (e.g. the track's total path length): shell volumes are annular
    cylinders, `pi * (r_out^2 - r_in^2) * length_nm`. This scorer does NOT
    itself restrict deposits by axial position -- `length_nm` is purely for
    volume/dose normalization, so pick it to match the axial extent the
    scored deposits actually cover.
    """

    def __init__(
        self,
        medium: Medium,
        r_edges: np.ndarray,
        length_nm: float,
        axis_point: np.ndarray | None = None,
        axis_direction: np.ndarray | None = None,
    ):
        self.medium = medium
        self.r_edges = np.asarray(r_edges, dtype=float)
        self.length_nm = length_nm
        self.axis_point = np.zeros(3) if axis_point is None else np.asarray(axis_point, dtype=float)
        direction = np.array([0.0, 0.0, 1.0]) if axis_direction is None else np.asarray(axis_direction, dtype=float)
        self.axis_direction = direction / np.linalg.norm(direction)
        self.energy_shell_eV = np.zeros(len(self.r_edges) - 1)  # [eV] per shell

    def _radial_distance(self, position: np.ndarray) -> float:
        # Remove the component along the axis; what is left is perpendicular.
        rel = position - self.axis_point
        along = np.dot(rel, self.axis_direction)
        perp = rel - along * self.axis_direction
        return float(np.linalg.norm(perp))

    def score(self, position: np.ndarray, energy: float) -> None:
        r = self._radial_distance(position)
        idx = int(np.searchsorted(self.r_edges, r, side="right") - 1)
        if 0 <= idx < len(self.energy_shell_eV):
            self.energy_shell_eV[idx] += energy

    def shell_volumes_nm3(self) -> np.ndarray:
        r_in, r_out = self.r_edges[:-1], self.r_edges[1:]
        # Annular cylinder: pi (r_out^2 - r_in^2) times the axial length.
        return np.pi * (r_out ** 2 - r_in ** 2) * self.length_nm

    def dose_Gy(self) -> np.ndarray:
        """Dose [Gy] per shell: energy[J] / mass[kg]."""
        energy_J = self.energy_shell_eV * constants.EV_TO_JOULE
        volume_cm3 = self.shell_volumes_nm3() / constants.CM3_TO_NM3
        mass_kg = self.medium.density_g_cm3 * volume_cm3 * 1.0e-3
        return energy_J / mass_kg

    def total_deposited_eV(self) -> float:
        return float(self.energy_shell_eV.sum())


def make_bragg_curve_scorer(
    medium: Medium, z_edges: np.ndarray, transverse_half_width_nm: float
) -> VoxelDoseScorer:
    """Depth-dose (Bragg curve) scorer: a single transverse bin spanning
    [-w, w] x [-w, w] nm, many z bins. `scorer.dose_Gy()[0, 0, :]` is the
    Bragg curve; `scorer.z_edges` gives the corresponding depth bin edges.

    The beam travels along +z. Choose w wide enough to catch the whole
    lateral spread of the track (in this code, caused by elastic
    scattering only: delta-electrons deposit on the proton's path), or the dose
    near the end of the range is underestimated.
    """
    w = transverse_half_width_nm
    return VoxelDoseScorer(
        medium, x_edges=np.array([-w, w]), y_edges=np.array([-w, w]), z_edges=z_edges
    )
