"""Use a full shared library, with LEPTONFF_LIBRARY set when needed."""
import math
from leptonff import Decay, DecayConfig

first = [[50, 40, 0, 30], [50, 30, 40, 0]]
second = [[50, -40, 0, -30], [50, -30, -40, 0]]
seeds = [42, 137]
with Decay(DecayConfig(grid_bins=1024)) as physics:
    with physics.with_observables({"pt": [25, 30, 35, 40, math.inf]}) as study:
        results = study.evaluate_batch(first, second, seeds=seeds, threads=2)
        print(study.output_labels)
        print([row.total.matched for row in results])
        print(study.resolved_configuration())
