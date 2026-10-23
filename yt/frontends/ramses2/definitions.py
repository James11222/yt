import re

# RAMSES II (and mini-ramses) outputs live in a folder named output_XXXXX
# (or backup_XXXXX for restart dumps) containing
# - info.txt, params.bin, namelist.txt, compilation.txt, timer.txt
# - amr.NNNNN                          the oct structure
# - hydro.NNNNN, grav.NNNNN, rt.NNNNN  cell-based fields (one block per oct)
# - part.NNNNN, star.NNNNN, sink.NNNNN, tree.NNNNN, trac.NNNNN
#                                      one file per particle type
# - <prefix>_header.txt                text descriptors of the binary files
# where NNNNN is the file index, running from 1 to nfile.
OUTPUT_DIR_RE = re.compile(r"output_(\d{5})")

RAMSES2_FILE_RE = re.compile(
    r"((amr|hydro|grav|rt|cr|part|star|sink|tree|trac)\.\d{5}"
    r"|(hydro|grav|rt|cr|part|star|sink|tree|trac)_header\.txt"
    r"|info\.txt|params\.bin)"
)

# Map of the fluid file prefixes to the yt field type they are loaded into
FLUID_FILES = {
    "hydro": "ramses2",
    "grav": "gravity",
    "rt": "rt",
    "cr": "cr",
}

# Map of the particle file prefixes to the yt particle type they are loaded into
PARTICLE_FILES = {
    "part": "io",
    "star": "star",
    "sink": "sink",
    "tree": "tree",
    "trac": "tracer",
}

# Names of the variables in <prefix>_header.txt mapped onto the names yt
# uses on disk for the fluid fields
FLUID_FIELD_NAMES = {
    "density": "Density",
    "velocity_x": "x-velocity",
    "velocity_y": "y-velocity",
    "velocity_z": "z-velocity",
    "thermal_pressure": "Pressure",
    "magnetic_field_x": "B_x",
    "magnetic_field_y": "B_y",
    "magnetic_field_z": "B_z",
    "metal_mass_fraction": "Metallicity",
    "potential": "Potential",
    "accel_x": "x-acceleration",
    "accel_y": "y-acceleration",
    "accel_z": "z-acceleration",
}

# Content of the particle files, in the order they are written to disk.
# For each entry of the "Particle fields" line of <prefix>_header.txt, give
# the yt field name(s) and the on-disk dtype. "id" fields are either 4 or 8
# bytes long depending on the compilation flags of the code: their size is
# inferred from the file size at read time.
PARTICLE_FIELD_DESCRIPTORS = {
    "pos": (
        ("particle_position_x", "particle_position_y", "particle_position_z"),
        "f4",
    ),
    "vel": (
        ("particle_velocity_x", "particle_velocity_y", "particle_velocity_z"),
        "f4",
    ),
    "mass": (("particle_mass",), "f4"),
    "potential": (("particle_potential",), "f4"),
    "metallicity": (("particle_metallicity",), "f4"),
    "accel": (
        (
            "particle_acceleration_x",
            "particle_acceleration_y",
            "particle_acceleration_z",
        ),
        "f4",
    ),
    "angmom": (
        (
            "particle_angular_momentum_x",
            "particle_angular_momentum_y",
            "particle_angular_momentum_z",
        ),
        "f4",
    ),
    "birth_date": (("particle_birth_time",), "f4"),
    "merging_date": (("particle_merging_time",), "f4"),
    "level": (("particle_refinement_level",), "i4"),
    "birth_id": (("particle_identity",), "id"),
    "merging_id": (("particle_merging_identity",), "id"),
    "tracking_id": (("particle_tracking_identity",), "id"),
    "size": (("particle_size",), "f4"),
    "charge": (("particle_charge",), "f4"),
}

field_aliases = {
    "standard_five": (
        "Density",
        "x-velocity",
        "y-velocity",
        "z-velocity",
        "Pressure",
    ),
}
