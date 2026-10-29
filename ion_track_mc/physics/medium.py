"""Medium: the target material, characterised by density and molar mass.

Water is treated as a single "effective target": the cross sections are
molecular (per H2O molecule), not built up from separate H and O atoms.
Only liquid water is defined; composite media such as tissue or PMMA are not
supported. See physics/README_physics.md §3.
"""
from dataclasses import dataclass
from pathlib import Path

from ion_track_mc import constants


@dataclass(frozen=True)
class Medium:
    """A target material, treated as a single effective molecular species.

    `target_label` is the {Target} token used in cross-section filenames
    (data/README_data_format.md naming convention, e.g. "H2O"); it is what
    `check_medium_matches_filename` cross-checks against.

    `target_nuclear_charge`/`target_nucleon_count` (Z_T, A_T) are used by
    the elastic-scattering formula (physics/README_physics.md §4.5). They
    default to `None`; WATER sets them (Z_T = 1 + 1 + 8 = 10, the total
    nuclear charge of H2O; A_T = 18 nucleons), and any other medium must set them before elastic
    scattering can be used with it.
    """
    name: str
    density_g_cm3: float
    molar_mass_g_mol: float
    target_label: str
    target_nuclear_charge: float | None = None
    target_nucleon_count: float | None = None

    @property
    def N(self) -> float:
        """Molecular number density [nm^-3] (internal unit policy).

        Multiplying N by a per-molecule cross section sigma [nm^2] gives the
        macroscopic cross section Sigma = N*sigma [nm^-1]: the expected number
        of collisions per nm of path. 1/Sigma is the mean free path.
        """
        # (g/cm^3) / (g/mol) = mol/cm^3; times Avogadro = molecules/cm^3.
        N_cm3 = (self.density_g_cm3 / self.molar_mass_g_mol) * constants.AVOGADRO
        return N_cm3 / constants.CM3_TO_NM3


WATER = Medium(name="water", density_g_cm3=1.0, molar_mass_g_mol=18.015, target_label="H2O",
               target_nuclear_charge=10.0, target_nucleon_count=18.0)


def target_label_from_filename(filename: str) -> str:
    """Extract the {Target} token from a cross-section filename.

    Per data/README_data_format.md:
        {Projectile}_+_{Target}_{channel}_TCS_{method}.dat
    """
    stem = Path(filename).stem
    parts = stem.split("_+_", 1)
    if len(parts) != 2:
        raise ValueError(
            f"filename does not follow the {{Projectile}}_+_{{Target}}_... "
            f"convention: {filename}"
        )
    target = parts[1].split("_", 1)[0]
    return target


def check_medium_matches_filename(medium: Medium, filename: str) -> None:
    """Consistency check: does a cross-section file's {Target} token match `medium`?

    Guards against wiring, e.g., a PMMA cross-section table to a Medium named
    water. Raises ValueError on mismatch.
    """
    target = target_label_from_filename(filename)
    if target != medium.target_label:
        raise ValueError(
            f"cross-section file target '{target}' does not match medium "
            f"'{medium.name}' (expected target label '{medium.target_label}')"
        )
