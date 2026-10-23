import warnings

from yt._typing import KnownFieldsT
from yt.fields.field_detector import FieldDetector
from yt.fields.field_info_container import FieldInfoContainer
from yt.utilities.physical_constants import (
    boltzmann_constant_cgs,
    mass_hydrogen_cgs,
)

b_units = "code_magnetic"
ra_units = "code_length / code_time**2"
rho_units = "code_density"
vel_units = "code_velocity"
pressure_units = "code_pressure"
specific_ener_units = "code_velocity**2"
ang_mom_units = "code_mass * code_velocity * code_length"


class RAMSES2FieldInfo(FieldInfoContainer):
    known_other_fields: KnownFieldsT = (
        ("Density", (rho_units, ["density"], None)),
        ("x-velocity", (vel_units, ["velocity_x"], None)),
        ("y-velocity", (vel_units, ["velocity_y"], None)),
        ("z-velocity", (vel_units, ["velocity_z"], None)),
        ("Pressure", (pressure_units, ["pressure"], None)),
        ("Metallicity", ("", ["metallicity"], None)),
        ("B_x", (b_units, ["magnetic_field_x"], None)),
        ("B_y", (b_units, ["magnetic_field_y"], None)),
        ("B_z", (b_units, ["magnetic_field_z"], None)),
        ("x-acceleration", (ra_units, ["acceleration_x"], None)),
        ("y-acceleration", (ra_units, ["acceleration_y"], None)),
        ("z-acceleration", (ra_units, ["acceleration_z"], None)),
        ("Potential", (specific_ener_units, ["potential"], None)),
    )
    known_particle_fields: KnownFieldsT = (
        ("particle_position_x", ("code_length", [], None)),
        ("particle_position_y", ("code_length", [], None)),
        ("particle_position_z", ("code_length", [], None)),
        ("particle_velocity_x", (vel_units, [], None)),
        ("particle_velocity_y", (vel_units, [], None)),
        ("particle_velocity_z", (vel_units, [], None)),
        ("particle_mass", ("code_mass", [], None)),
        ("particle_potential", (specific_ener_units, [], None)),
        ("particle_metallicity", ("", [], None)),
        ("particle_acceleration_x", (ra_units, [], None)),
        ("particle_acceleration_y", (ra_units, [], None)),
        ("particle_acceleration_z", (ra_units, [], None)),
        ("particle_angular_momentum_x", (ang_mom_units, [], None)),
        ("particle_angular_momentum_y", (ang_mom_units, [], None)),
        ("particle_angular_momentum_z", (ang_mom_units, [], None)),
        ("particle_birth_time", ("code_time", ["creation_time"], None)),
        ("particle_merging_time", ("code_time", [], None)),
        ("particle_refinement_level", ("", [], None)),
        ("particle_identity", ("", ["particle_index"], None)),
        ("particle_merging_identity", ("", [], None)),
        ("particle_tracking_identity", ("", [], None)),
    )

    def setup_particle_fields(self, ptype):
        super().setup_particle_fields(ptype)

        if (ptype, "particle_birth_time") not in self.field_list:
            return

        def _age(data):
            return data.ds.current_time - data[ptype, "particle_birth_time"]

        self.add_field(
            (ptype, "age"),
            sampling_type="particle",
            function=_age,
            units=self.ds.unit_system["time"],
        )

    def setup_fluid_fields(self):
        def _temperature_over_mu(data):
            rv = data["gas", "pressure"] / data["gas", "density"]
            rv *= mass_hydrogen_cgs / boltzmann_constant_cgs
            return rv

        self.add_field(
            ("gas", "temperature_over_mu"),
            sampling_type="cell",
            function=_temperature_over_mu,
            units=self.ds.unit_system["temperature"],
        )

        def _temperature(data):
            if not isinstance(data, FieldDetector):
                warnings.warn(
                    "The mean molecular weight is not stored in RAMSES II outputs. "
                    "yt will return T/µ instead of T — this is equivalent to "
                    "assuming µ=1.0. To suppress this, derive the temperature "
                    "from temperature_over_mu with some values for mu.",
                    category=RuntimeWarning,
                    stacklevel=1,
                )
            return data["gas", "temperature_over_mu"]

        self.add_field(
            ("gas", "temperature"),
            sampling_type="cell",
            function=_temperature,
            units=self.ds.unit_system["temperature"],
        )
