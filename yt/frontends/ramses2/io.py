import os
from collections import defaultdict

import numpy as np

from yt.utilities.io_handler import BaseIOHandler
from yt.utilities.logger import ytLogger as mylog

from .definitions import FLUID_FILES, PARTICLE_FIELD_DESCRIPTORS, PARTICLE_FILES

_PTYPE_TO_PREFIX = {ptype: prefix for prefix, ptype in PARTICLE_FILES.items()}
_FTYPE_TO_PREFIX = {ftype: prefix for prefix, ftype in FLUID_FILES.items()}


class IOHandlerRAMSES2(BaseIOHandler):
    _dataset_type = "ramses2"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        tr = defaultdict(list)
        ftypes = {ft for ft, _ in fields}
        for chunk in chunks:
            for ft in ftypes:
                field_subs = [f for f in fields if f[0] == ft]
                field_list = [f for f in self.ds._fields_in_file if f[0] == ft]
                for subset in chunk.objs:
                    fn = subset.domain.fluid_fn(_FTYPE_TO_PREFIX[ft])
                    if not os.path.exists(fn):
                        continue
                    rv = subset.fill(fn, field_subs, selector, field_list)
                    for ftype, fname in field_subs:
                        d = rv.pop(fname)
                        if d.size == 0:
                            continue
                        mylog.debug(
                            "Filling %s with %s (%0.3e %0.3e) (%s zones)",
                            fname,
                            d.size,
                            d.min(),
                            d.max(),
                            d.size,
                        )
                        tr[ftype, fname].append(d)
        d = {}
        for field in fields:
            tmp = tr.pop(field, None)
            d[field] = np.concatenate(tmp) if tmp else np.empty(0, dtype="d")
        return d

    @staticmethod
    def _particle_count(fn):
        if not os.path.exists(fn):
            return 0
        return int(np.fromfile(fn, dtype="<i4", count=2)[1])

    def _particle_layout(self, fn, ptype):
        """Return the number of particles and the offset and dtype of each
        field in a particle file.

        Particle files contain a header (ndim, npart) followed by one
        contiguous array per field, in the order given by <prefix>_header.txt.
        """
        npart = self._particle_count(fn)
        names = self.ds.parameters[f"{_PTYPE_TO_PREFIX[ptype]}_fields"]
        if npart == 0:
            return 0, {}

        nbytes = os.path.getsize(fn) - 8
        n_id = sum(PARTICLE_FIELD_DESCRIPTORS[n][1] == "id" for n in names)
        n_other = sum(
            len(PARTICLE_FIELD_DESCRIPTORS[n][0])
            for n in names
            if PARTICLE_FIELD_DESCRIPTORS[n][1] != "id"
        )
        id_size = (nbytes // npart - 4 * n_other) // max(n_id, 1)
        if id_size not in (4, 8) or nbytes != npart * (4 * n_other + id_size * n_id):
            raise RuntimeError(f"Could not infer the layout of particle file {fn}")

        layout = {}
        offset = 8
        for name in names:
            fnames, kind = PARTICLE_FIELD_DESCRIPTORS[name]
            dtype = f"<i{id_size}" if kind == "id" else f"<{kind}"
            for fname in fnames:
                layout[fname] = (offset, dtype)
                offset += npart * np.dtype(dtype).itemsize
        return npart, layout

    def _read_particle_file(self, fn, ptype, fields):
        npart, layout = self._particle_layout(fn, ptype)
        ds = self.ds
        tr = {}
        for fname in fields:
            offset, dtype = layout[fname]
            data = np.fromfile(fn, dtype=dtype, count=npart, offset=offset)
            data = data.astype("float64")
            if fname.startswith("particle_position_"):
                # Positions are stored in [0, boxlen]
                data /= ds.parameters["boxlen"]
            elif fname in ("particle_birth_time", "particle_merging_time"):
                if ds.cosmological_simulation:
                    # Times are stored as the proper time since today in
                    # units of 1/H0 = unit_t / aexp**2
                    aexp = ds.parameters["aexp"]
                    data = ds.current_time.d + (data - ds.parameters["texp"]) / aexp**2
            tr[fname] = data
        return tr

    def _read_particle_coords(self, chunks, ptf):
        pn = "particle_position_%s"
        for chunk in chunks:
            for subset in chunk.objs:
                for ptype in sorted(ptf):
                    fn = subset.domain.particle_fn(_PTYPE_TO_PREFIX[ptype])
                    if self._particle_count(fn) == 0:
                        continue
                    rv = self._read_particle_file(fn, ptype, [pn % ax for ax in "xyz"])
                    yield ptype, (rv[pn % "x"], rv[pn % "y"], rv[pn % "z"]), 0.0

    def _read_particle_fields(self, chunks, ptf, selector):
        pn = "particle_position_%s"
        for chunk in chunks:
            for subset in chunk.objs:
                for ptype, field_list in sorted(ptf.items()):
                    fn = subset.domain.particle_fn(_PTYPE_TO_PREFIX[ptype])
                    if self._particle_count(fn) == 0:
                        continue
                    pos_fields = [pn % ax for ax in "xyz"]
                    rv = self._read_particle_file(
                        fn, ptype, set(field_list) | set(pos_fields)
                    )
                    mask = selector.select_points(
                        rv[pn % "x"], rv[pn % "y"], rv[pn % "z"], 0.0
                    )
                    if mask is None:
                        continue
                    for fname in field_list:
                        yield (ptype, fname), rv[fname][mask]
