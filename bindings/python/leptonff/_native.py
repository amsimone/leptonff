"""Private ctypes declarations for the versioned, double-based C interface."""
import ctypes as C
from ctypes.util import find_library
import os
from pathlib import Path


class FFConfig(C.Structure):
    _fields_ = [("struct_size", C.c_size_t), ("abi_version", C.c_uint)] + [
        (name, C.c_int) for name in ("model", "alpha_scheme", "running", "flavours")
    ] + [(name, C.c_double) for name in (
        "alpha", "alpha_reference_scale", "muon_mass", "initial_scale", "final_scale"
    )] + [("grid_bins", C.c_size_t)] + [(name, C.c_double) for name in (
        "log_span", "absolute_tolerance", "relative_tolerance"
    )]


class DecayConfig(C.Structure):
    _fields_ = [("struct_size", C.c_size_t), ("abi_version", C.c_uint),
                ("process", C.c_int), ("model", C.c_int)] + [
        (name, C.c_double) for name in ("alpha", "muon_mass", "mu_frag")
    ] + [("grid_bins", C.c_size_t), ("log_span", C.c_double), ("hard_samples", C.c_uint)] + [
        (name, C.c_double) for name in (
            "soft_cutoff", "collinear_cutoff", "minimum_z", "muon_pt_min", "muon_eta_max",
            "neutrino_pt_min", "transverse_mass_min", "mass_min", "mass_max", "hadronic_matching_scale"
        )
    ]


class Observable(C.Structure):
    _fields_ = [("observable", C.c_char_p), ("edges", C.POINTER(C.c_double)), ("edge_count", C.c_size_t)]


class DecayTheory(C.Structure):
    _fields_ = [("struct_size", C.c_size_t), ("abi_version", C.c_uint),
                ("alpha_scheme", C.c_int), ("running", C.c_int)] + [
        (name, C.c_double) for name in ("alpha", "alpha_reference_scale", "initial_scale")]


class Prediction(C.Structure):
    _fields_ = [(name, C.c_double) for name in ("born", "order1", "matched", "order2", "higher")]


class Diagnostics(C.Structure):
    _fields_ = [("abi_version", C.c_uint), ("photonic_unitary_conditioning", C.c_int)] + [
        (name, C.c_double) for name in (
            "a_initial", "source_number_moment", "raw_inclusive_factor", "wrap_bound", "tail_bound"
        )
    ]


Function = C.CFUNCTYPE(C.c_double, C.c_double, C.c_void_p)


class Library:
    def __init__(self, path=None, *, decay=False):
        name = path or os.environ.get("LEPTONFF_LIBRARY") or find_library("leptonff")
        if not name and not decay:
            name = find_library("leptonff_core")
        if not name:
            raise OSError("Install a shared leptonff library and set LEPTONFF_LIBRARY to its path")
        expanded = Path(name).expanduser()
        self.lib = C.CDLL(str(expanded.resolve()) if expanded.exists() else str(name))
        self._declare("leptonff_last_error", C.c_char_p)
        self._declare("leptonff_version", C.c_char_p)
        self.version = self.lib.leptonff_version().decode("ascii")
        self._declare("leptonff_fragmentation_config_default", None, C.POINTER(FFConfig))
        self._declare("leptonff_fragmentation_create", C.c_int, C.POINTER(FFConfig), C.POINTER(C.c_void_p))
        self._declare("leptonff_fragmentation_destroy", None, C.c_void_p)
        self._declare("leptonff_fragmentation_moment", C.c_int, C.c_void_p, C.c_int, C.c_int, C.c_double, C.POINTER(C.c_double))
        self._declare("leptonff_fragmentation_integral", C.c_int, C.c_void_p, C.c_int, C.c_int, C.c_double, C.c_double, C.POINTER(C.c_double))
        self._declare("leptonff_fragmentation_convolve", C.c_int, C.c_void_p, C.c_int, C.c_int, C.c_double, Function, C.c_void_p, C.POINTER(C.c_double))
        if decay:
            if not hasattr(self.lib, "leptonff_decay_create"):
                raise OSError("Decay calculations require the full shared leptonff library")
            self._declare("leptonff_decay_config_default", None, C.POINTER(DecayConfig))
            self._declare("leptonff_decay_create", C.c_int, C.POINTER(DecayConfig), C.POINTER(Observable), C.c_size_t, C.POINTER(C.c_void_p))
            self._declare("leptonff_decay_destroy", None, C.c_void_p)
            self._declare("leptonff_decay_output_count", C.c_size_t, C.c_void_p)
            self._declare("leptonff_decay_evaluate", C.c_int, C.c_void_p, C.POINTER(C.c_double), C.POINTER(C.c_double), C.c_uint64, C.POINTER(Prediction), C.c_size_t)
            self._declare("leptonff_decay_get_diagnostics", C.c_int, C.c_void_p, C.POINTER(Diagnostics))
            self._declare("leptonff_decay_response", C.c_int, C.c_void_p, C.POINTER(C.c_double), C.c_uint64, C.c_uint, C.c_uint, C.POINTER(C.c_double), C.c_size_t)
            self._declare("leptonff_decay_theory_default", None, C.POINTER(DecayTheory))
            self._declare("leptonff_decay_create_with_theory", C.c_int, C.POINTER(DecayConfig), C.POINTER(DecayTheory), C.POINTER(Observable), C.c_size_t, C.POINTER(C.c_void_p))
            self._declare("leptonff_decay_with_observables", C.c_int, C.c_void_p, C.POINTER(DecayConfig), C.POINTER(Observable), C.c_size_t, C.POINTER(C.c_void_p))
            self._declare("leptonff_decay_evaluate_batch", C.c_int, C.c_void_p, C.POINTER(C.c_double), C.POINTER(C.c_double), C.POINTER(C.c_uint64), C.c_size_t, C.c_uint, C.POINTER(Prediction), C.c_size_t)
            self._declare("leptonff_decay_configuration_json", C.c_int, C.c_void_p, C.c_void_p, C.POINTER(C.c_size_t))
            self._declare("leptonff_decay_output_label", C.c_int, C.c_void_p, C.c_size_t, C.c_void_p, C.POINTER(C.c_size_t))

    def _declare(self, name, result, *args):
        function = getattr(self.lib, name)
        function.restype = result
        function.argtypes = list(args)

    def check(self, status):
        if status:
            raise ValueError(self.lib.leptonff_last_error().decode("utf-8", errors="replace"))


def configuration(library, native_type, defaults, values):
    native = native_type()
    defaults(C.byref(native))
    if native.struct_size != C.sizeof(native) or native.abi_version != 1:
        raise RuntimeError("Incompatible leptonff C ABI layout/version")
    for name, value in vars(values).items():
        field_type = dict(native_type._fields_)[name]
        if field_type in (C.c_size_t, C.c_uint) and (
            not isinstance(value, int) or not 0 <= value < 2 ** (8 * C.sizeof(field_type))
        ):
            raise ValueError(f"{name} must be a non-negative integer fitting the native ABI")
        setattr(native, name, value)
    return native
