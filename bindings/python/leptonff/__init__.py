"""Numerical bindings to leptonff's C ABI; a shared native library is required.

Fragmentation preserves delta/plus contributions through moments, integrals
and convolutions. Decay exposes the separately matched W/Z prescription.
"""
import ctypes as C
from dataclasses import dataclass, replace
from enum import IntEnum
import math
import json
from . import _native as native

__version__ = "0.27.0"
__all__ = ["Species", "FFModel", "AlphaScheme", "Running", "Flavours", "Process",
           "DecayModel", "FragmentationConfig", "Fragmentation", "DecayConfig", "DecayTheory", "Decay", "Prediction", "Spectrum"]


class Species(IntEnum):
    MUON_MINUS = 0
    MUON_PLUS = 1
    PHOTON = 2


class FFModel(IntEnum):
    ORDER1 = 1
    COUPLED_NLL = 2


class AlphaScheme(IntEnum):
    ON_SHELL = 0
    MSBAR = 1


class Running(IntEnum):
    FIXED = 0
    ONE_LOOP = 1
    TWO_LOOP = 2


class Flavours(IntEnum):
    ONE_FAMILY = 0
    LEPTON_THRESHOLDS = 1


class Process(IntEnum):
    W_PLUS = 24
    W_MINUS = -24
    Z = 23


class DecayModel(IntEnum):
    PHOTONIC_LL = 1
    PHOTONIC_NLL = 2
    LEPTONIC_NLL = 3
    QED_HVP = 4


@dataclass(frozen=True)
class FragmentationConfig:
    model: FFModel = FFModel.ORDER1
    alpha_scheme: AlphaScheme = AlphaScheme.MSBAR
    running: Running = Running.TWO_LOOP
    flavours: Flavours = Flavours.ONE_FAMILY
    alpha: float = 1 / 137.035999084
    alpha_reference_scale: float = .1056583755
    muon_mass: float = .1056583755
    initial_scale: float = .1056583755
    final_scale: float = 91.1876
    grid_bins: int = 0
    log_span: float = 40
    absolute_tolerance: float = 1e-10
    relative_tolerance: float = 1e-10


class _Context:
    @property
    def native_version(self):
        """Version of the linked native library (separate from the binding)."""
        return self._library.version

    def _require_open(self):
        if not self._handle.value:
            raise RuntimeError("leptonff context is closed")

    def close(self):
        if getattr(self, "_handle", None) and self._handle.value:
            self._destroy(self._handle)
            self._handle.value = None

    def __enter__(self):
        self._require_open()
        return self

    def __exit__(self, *_):
        self.close()

    def __del__(self):
        self.close()


class Fragmentation(_Context):
    """Prepare once and reuse. Coupled x-space bins require grid_bins>0.

    Moments use the continuum solver. Coupled bin integrals use the signed
    energy-weighted lattice and need convergence checks for sharp bins.
    """
    def __init__(self, config=None, *, library=None):
        self.config = config or FragmentationConfig()
        self._library = native.Library(library)
        lib = self._library.lib
        self._destroy = lib.leptonff_fragmentation_destroy
        self._handle = C.c_void_p()
        cfg = native.configuration(self._library, native.FFConfig, lib.leptonff_fragmentation_config_default, self.config)
        self._library.check(lib.leptonff_fragmentation_create(C.byref(cfg), C.byref(self._handle)))

    def moment(self, N, parent=Species.MUON_MINUS, observed=Species.MUON_MINUS):
        self._require_open()
        result = C.c_double()
        self._library.check(self._library.lib.leptonff_fragmentation_moment(self._handle, parent, observed, N, C.byref(result)))
        return result.value

    def integral(self, lower, upper, parent=Species.MUON_MINUS, observed=Species.MUON_MINUS):
        self._require_open()
        result = C.c_double()
        self._library.check(self._library.lib.leptonff_fragmentation_integral(self._handle, parent, observed, lower, upper, C.byref(result)))
        return result.value

    def convolve(self, function, x, parent=Species.MUON_MINUS, observed=Species.MUON_MINUS):
        self._require_open()
        failure = []

        @native.Function
        def callback(z, _):
            try:
                return float(function(z))
            except Exception as error:
                failure.append(error)
                return math.nan

        result = C.c_double()
        status = self._library.lib.leptonff_fragmentation_convolve(self._handle, parent, observed, x, callback, None, C.byref(result))
        if failure:
            raise failure[0]
        self._library.check(status)
        return result.value


@dataclass(frozen=True)
class DecayConfig:
    """ABI-v1 decay presets: alpha is fixed for photonic modes and on-shell
    for leptonic/HVP modes. Photonic NLL is a unitary-conditioned shape.
    """
    process: Process = Process.W_PLUS
    model: DecayModel = DecayModel.PHOTONIC_NLL
    alpha: float = 1 / 137.035999084
    muon_mass: float = .1056583755
    mu_frag: float = 80.385
    grid_bins: int = 4096
    log_span: float = 40
    hard_samples: int = 16
    soft_cutoff: float = 1e-5
    collinear_cutoff: float = 1e-5
    minimum_z: float = 1e-4
    muon_pt_min: float = 25
    muon_eta_max: float = 2.4
    neutrino_pt_min: float = 0
    transverse_mass_min: float = 0
    mass_min: float = 60
    mass_max: float = 120
    hadronic_matching_scale: float = 2


@dataclass(frozen=True)
class Prediction:
    """Related predictions: order1 includes Born; higher=matched-order1.

    order2 is the retained second-order diagnostic, not a complete NNLO result.
    """
    born: float
    order1: float
    matched: float
    order2: float
    higher: float


@dataclass(frozen=True)
class DecayTheory:
    """Explicit alpha input; MSbar reference scale must equal initial_scale.

    Photonic modes require FIXED running and initial_scale=muon_mass.
    The leptonic preset supports independent running order and alpha scheme.
    """
    alpha_scheme: AlphaScheme = AlphaScheme.MSBAR
    running: Running = Running.FIXED
    alpha: float = 1 / 137.035999084
    alpha_reference_scale: float = .1056583755
    initial_scale: float = .1056583755


@dataclass(frozen=True)
class Spectrum:
    total: Prediction
    bins: tuple[Prediction, ...]


def _momentum(values):
    values = tuple(values)
    if len(values) != 4 or not all(math.isfinite(v) for v in values):
        raise ValueError("Momentum must contain four finite values (energy, px, py, pz) in GeV")
    return (C.c_double * 4)(*values)


def _seed(value):
    if not isinstance(value, int) or not 0 <= value < 2**64:
        raise ValueError("Radiation seed must be an unsigned 64-bit integer")
    return value


class Decay(_Context):
    """Immutable matched W/Z context. Host weights and covariance are external.

    Observables is a mapping of names to edges, e.g. {'pt': [25,30,40,math.inf]}.
    Use the same radiation seed for every mapped point in a subtraction group.
    """
    def __init__(self, config=None, observables=None, *, theory=None, library=None):
        self.config = config or DecayConfig()
        self.theory = theory
        self._library = native.Library(library, decay=True)
        lib = self._library.lib
        self._destroy = lib.leptonff_decay_destroy
        self._handle = C.c_void_p()
        cfg = native.configuration(self._library, native.DecayConfig, lib.leptonff_decay_config_default, self.config)
        array, definitions = self._definitions(observables)
        if theory is None:
            self._library.check(lib.leptonff_decay_create(C.byref(cfg), array, len(definitions), C.byref(self._handle)))
        else:
            t = native.configuration(self._library, native.DecayTheory, lib.leptonff_decay_theory_default, theory)
            self._library.check(lib.leptonff_decay_create_with_theory(C.byref(cfg), C.byref(t), array, len(definitions), C.byref(self._handle)))
        self.output_count = lib.leptonff_decay_output_count(self._handle)

    @staticmethod
    def _definitions(observables):
        definitions = []
        for name, values in (observables or {}).items():
            edges = (C.c_double * len(values))(*values)
            definitions.append(native.Observable(name.encode(), edges, len(values)))
        array = (native.Observable * len(definitions))(*definitions) if definitions else None
        return array, definitions

    def with_observables(self, observables=None, **cuts):
        """Create an independent measurement view sharing the prepared physics."""
        self._require_open()
        allowed = {"muon_pt_min", "muon_eta_max", "neutrino_pt_min", "transverse_mass_min", "mass_min", "mass_max"}
        if cuts.keys() - allowed:
            raise ValueError("A measurement view may only change cuts and binning")
        child = object.__new__(type(self))
        child.config = replace(self.config, **cuts)
        child.theory = self.theory
        child._library, child._destroy, child._handle = self._library, self._destroy, C.c_void_p()
        lib = self._library.lib
        cfg = native.configuration(self._library, native.DecayConfig, lib.leptonff_decay_config_default, child.config)
        array, definitions = self._definitions(observables)
        self._library.check(lib.leptonff_decay_with_observables(self._handle, C.byref(cfg), array, len(definitions), C.byref(child._handle)))
        child.output_count = lib.leptonff_decay_output_count(child._handle)
        return child

    def _string(self, function, *arguments):
        self._require_open()
        capacity = C.c_size_t()
        self._library.check(function(self._handle, *arguments, None, C.byref(capacity)))
        result = C.create_string_buffer(capacity.value)
        self._library.check(function(self._handle, *arguments, result, C.byref(capacity)))
        return result.value.decode("utf-8")

    @property
    def output_labels(self):
        return tuple(self._string(self._library.lib.leptonff_decay_output_label, i) for i in range(self.output_count))

    def resolved_configuration(self):
        return json.loads(self._string(self._library.lib.leptonff_decay_configuration_json))

    def evaluate_batch(self, first, second, *, seeds=None, threads=1):
        """Evaluate complete momentum arrays in one native call; no NumPy needed."""
        self._require_open()
        first, second = list(first), list(second)
        if len(first) != len(second):
            raise ValueError("First and second momentum arrays must have the same length")
        count = len(first)
        seeds = list(seeds) if seeds is not None else [0] * count
        if len(seeds) != count:
            raise ValueError("One radiation seed is required per point")
        if not isinstance(threads, int) or not 1 <= threads <= 1024:
            raise ValueError("threads must be in [1,1024]")
        a = (C.c_double * (4 * count))(*(v for p in first for v in _momentum(p)))
        b = (C.c_double * (4 * count))(*(v for p in second for v in _momentum(p)))
        stream = (C.c_uint64 * count)(*(_seed(s) for s in seeds))
        outputs = (native.Prediction * (count * self.output_count))()
        self._library.check(self._library.lib.leptonff_decay_evaluate_batch(
            self._handle, a, b, stream, count, threads, outputs, len(outputs)))
        results = []
        for i in range(count):
            row = tuple(Prediction(*(getattr(p, n) for n, _ in p._fields_))
                        for p in outputs[i*self.output_count:(i+1)*self.output_count])
            results.append(Spectrum(row[0], row[1:]))
        return tuple(results)

    def evaluate_array(self, first, second, *, seeds=None, threads=1):
        """NumPy float64 output (points, outputs, 5), in Prediction field order.

        Contiguous inputs pass directly to C; other layouts are copied once.
        output_labels names axis 1; arrays use (E,px,py,pz) in GeV.
        """
        import numpy as np
        self._require_open()
        a = np.ascontiguousarray(first, dtype=np.float64)
        b = np.ascontiguousarray(second, dtype=np.float64)
        if a.ndim != 2 or a.shape[1] != 4 or b.shape != a.shape:
            raise ValueError("Momenta must have matching shape (points,4)")
        if not isinstance(threads, int) or not 1 <= threads <= 1024:
            raise ValueError("threads must be in [1,1024]")
        if seeds is None:
            stream = np.zeros(len(a), dtype=np.uint64)
        else:
            # Validate before uint conversion, preventing negative/wrapped seeds.
            values = list(seeds)
            if len(values) != len(a):
                raise ValueError("One radiation seed is required per point")
            import operator
            stream = np.array([_seed(operator.index(s)) for s in values], dtype=np.uint64)
        result = np.empty((len(a), self.output_count, 5), dtype=np.float64)
        pointer = lambda x: x.ctypes.data_as(C.POINTER(C.c_double))
        self._library.check(self._library.lib.leptonff_decay_evaluate_batch(
            self._handle, pointer(a), pointer(b), stream.ctypes.data_as(C.POINTER(C.c_uint64)),
            len(a), threads, result.ctypes.data_as(C.POINTER(native.Prediction)), result.size // 5))
        return result

    def evaluate(self, first, second, *, seed=0):
        self._require_open()
        outputs = (native.Prediction * self.output_count)()
        self._library.check(self._library.lib.leptonff_decay_evaluate(
            self._handle, _momentum(first), _momentum(second), _seed(seed), outputs, len(outputs)))
        predictions = tuple(Prediction(*(getattr(p, name) for name, _ in p._fields_)) for p in outputs)
        return Spectrum(predictions[0], predictions[1:])

    def diagnostics(self):
        self._require_open()
        result = native.Diagnostics()
        self._library.check(self._library.lib.leptonff_decay_get_diagnostics(self._handle, C.byref(result)))
        return {name: getattr(result, name) for name, _ in result._fields_}

    def response(self, parent, *, seed=0, cosine_nodes=16, azimuth_nodes=32):
        """Return one tuple per output: five 3x3 complex Cartesian tensors.

        Output zero is the fiducial total, followed by physical bins. Component
        order is born/order1/matched/order2/higher. Refine angular nodes and
        hard_samples independently before assigning numerical precision.
        """
        self._require_open()
        if not isinstance(cosine_nodes, int) or not isinstance(azimuth_nodes, int) or not (
            0 < cosine_nodes < 2**32 and 0 < azimuth_nodes < 2**32
        ):
            raise ValueError("Angular node counts must be positive unsigned integers")
        outputs = (C.c_double * (90 * self.output_count))()
        self._library.check(self._library.lib.leptonff_decay_response(
            self._handle, _momentum(parent), _seed(seed), cosine_nodes, azimuth_nodes, outputs, len(outputs)))
        values = iter(outputs)
        return tuple(tuple(tuple(tuple(complex(next(values), next(values)) for _ in range(3))
                                 for _ in range(3)) for _ in range(5)) for _ in range(self.output_count))
