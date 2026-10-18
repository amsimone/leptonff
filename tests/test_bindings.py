"""Exercise Python against the actual native C ABI, including error paths."""
from dataclasses import replace
import json
import math
from pathlib import Path
import sys

root = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(root / "bindings/python"))
from leptonff import (__version__, Fragmentation, FragmentationConfig, FFModel, Species,
                      Decay, DecayConfig, DecayTheory, DecayModel, Running, Process)

library = sys.argv[1]
reference = json.loads((root / "benchmarks/canonical/order1.json").read_text())
with Fragmentation(library=library) as ff:
    assert ff.native_version == __version__
    for N in (1, 2, 3):
        assert abs(ff.moment(N) - reference["values"][f"ff_N{N}"]) < 1e-10
    assert math.isfinite(ff.integral(.9, 1))
    assert abs(sum(ff.integral(lo, hi) for lo, hi in ((0, .3), (.3, .9), (.9, 1))) - 1) < 1e-10
    assert math.isfinite(ff.convolve(lambda z: z*(1-z), .4))
    try:
        ff.moment(1, observed=Species.PHOTON)
    except ValueError:
        pass
    else:
        raise AssertionError("Soft photon number moment accepted")
    try:
        ff.convolve(lambda z: 1 / 0, .4)
    except ZeroDivisionError:
        pass
    else:
        raise AssertionError("Python callback exception lost")
try:
    ff.moment(2)
except RuntimeError:
    pass
else:
    raise AssertionError("Closed native handle accepted")
with Fragmentation(FragmentationConfig(alpha=0), library=library) as ff:
    assert abs(ff.convolve(lambda z: z*(1-z), .4) - .24) < 1e-14
with Fragmentation(FragmentationConfig(model=FFModel.COUPLED_NLL), library=library) as ff:
    assert math.isfinite(ff.moment(2))
try:
    Fragmentation(FragmentationConfig(alpha_reference_scale=2), library=library)
except ValueError:
    pass
else:
    raise AssertionError("Implicit MSbar input-scale change accepted")

if "--decay" in sys.argv:
    # Every generator-independent process ID must reach the corresponding
    # native preset, and zero alpha must preserve Born predictions bin by bin.
    config = DecayConfig(alpha=0, grid_bins=128, hard_samples=2)
    for process in Process:
        observable = "mass" if process == Process.Z else "pt"
        edges = [60, 70, 80, 90, 100, 110, 120] if process == Process.Z else [25, 30, 35, 40, 45, 60, math.inf]
        with Decay(replace(config, process=process), {observable: edges}, library=library) as decay:
            result = decay.evaluate([50, 40, 0, 30], [50, -40, 0, -30], seed=42)
            assert result.total.born == result.total.matched == 1
            assert decay.output_labels[0] == "total" and len(decay.output_labels) == decay.output_count
            assert decay.resolved_configuration()["alpha_msbar_mu0"] == 0
            batch = decay.evaluate_batch([[50,40,0,30]]*3, [[50,-40,0,-30]]*3, seeds=[42]*3, threads=2)
            assert all(row == result for row in batch)
            assert decay.evaluate_batch([], [], threads=2) == ()
            try:
                import numpy as np
            except ImportError:
                pass
            else:
                array = decay.evaluate_array([[50,40,0,30]]*3, [[50,-40,0,-30]]*3, seeds=[42]*3, threads=2)
                assert array.shape == (3, decay.output_count, 5)
                assert np.array_equal(array[:,0,2], np.ones(3))
            assert abs(sum(p.matched for p in result.bins) - 1) < 1e-14
            diagnostics = decay.diagnostics()
            assert diagnostics["photonic_unitary_conditioning"] == 1
            assert abs(diagnostics["source_number_moment"] + .719311605774509) < 1e-12
            tensors = decay.response([100, 0, 0, 0], cosine_nodes=4, azimuth_nodes=8)
            assert len(tensors) == 7 and len(tensors[0]) == 5
            for output in tensors:
                for tensor in output:
                    for i in range(3):
                        for j in range(3):
                            assert abs(tensor[i][j] - tensor[j][i].conjugate()) < 1e-12
            try:
                decay.evaluate([50], [50, -40, 0, -30])
            except ValueError:
                pass
            else:
                raise AssertionError("Malformed momentum accepted")
            view = decay.with_observables({observable: [edges[0], edges[-1]]})
        with view:
            assert view.output_count == 2
            assert view.evaluate([50,40,0,30], [50,-40,0,-30]).total == result.total
    with Decay(replace(config, model=DecayModel.LEPTONIC_NLL), {"pt": [25,math.inf]},
               theory=DecayTheory(running=Running.TWO_LOOP, initial_scale=.2, alpha_reference_scale=.2), library=library) as explicit:
        assert explicit.resolved_configuration()["mu0"] == .2
print("Python/native fragmentation, matching, lifecycle and exception checks passed")
