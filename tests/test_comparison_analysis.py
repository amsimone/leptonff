"""Correlated tool ratios and cancellation of identical normalised shapes."""
import math
import csv
import tempfile
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'adapters/mcfm'))
from analyse_comparison import aggregate, chunk_difference, quadratic_diagnostic, ratio, rows_for
from analyse_spectra import COMPONENTS, RATES, definitions

# A common fluctuating sample has an exactly known ratio despite nonzero
# marginal errors. Losing its cross covariance would give a spurious error.
assert ratio([10., 20.], [[4., 8.], [8., 16.]], 1, 0) == [2., 0.]
assert ratio([1., -.1], [[.1, 0.], [0., 1.]], 0, 1) is None
assert ratio([10., -20.], [[.01, 0.], [0., .01]], 0, 1)[0] == -.5

edges = definitions('wplus')
labels = RATES+[f'sp_{obs}{k+1:02d}_{c}' for obs, values in edges.items()
                for k in range(len(values)-1) for c in COMPONENTS]
vectors = []
for event, weight in enumerate([1.2, -.3, .9, 2.]):
    row = dict.fromkeys(labels, 0.)
    for obs, values in edges.items():
        n = len(values)-1
        for k in range(n):
            # Preserve a migration fluctuation that cannot be reduced to a
            # common overall normalisation. PHOTOS is a constant multiple
            # of the analytic density, so their normalised shapes coincide.
            m = weight*(1+.2*((event+k) % n))
            b = weight*(1+.1*k)
            components = [b, m-.1*b, m, .01*b, .1*b, b, 1.03*m, .03*m]
            for c, x in zip(COMPONENTS, components):
                row[f'sp_{obs}{k+1:02d}_{c}'] = x
        for c, rate in zip(COMPONENTS, RATES):
            # Total rate must be identical for every observable: rescale each
            # bin vector to the pT total before constructing event products.
            total = sum(row[f'sp_{obs}{k+1:02d}_{c}'] for k in range(n))
            target = row[rate] if obs != 'pt' else total
            for k in range(n):
                row[f'sp_{obs}{k+1:02d}_{c}'] *= target/total
            row[rate] = target
    vectors.append([row[label] for label in labels])
vectors *= 20
mean = [sum(event[i] for event in vectors)/len(vectors) for i in range(len(labels))]
cov = [[sum((event[i]-mean[i])*(event[j]-mean[j]) for event in vectors)/(len(vectors)*(len(vectors)-1))
        for j in range(len(labels))] for i in range(len(labels))]
rows = rows_for('wplus', labels, mean, cov)
for data in rows.values():
    assert abs(sum(r['normalised_fractions']['match'][0] for r in data['bins'])-1) < 1e-14
    for row in data['bins']:
        assert abs(row['photos_over_match'][0]-1.03) < 1e-14
        assert row['photos_over_match'][1] < 1e-7
        delta, sigma = row['normalised_photos_minus_match']
        assert abs(delta) < 1e-14 and sigma < 1e-7

reports = [{'labels': ['x'], 'mean_pb': [10.], 'covariance_pb2': [[4.]]},
           {'labels': ['x'], 'mean_pb': [12.], 'covariance_pb2': [[9.]]}]
pooled_mean, pooled_cov = aggregate(reports)
assert pooled_mean == [11.] and pooled_cov == [[3.25]]
try:
    import numpy
except ImportError:
    pass  # NumPy is optional for the core suite and required for campaign diagnostics.
else:
    for scale in [1e-120, 1., 1e120]:
        q = quadratic_diagnostic([.5*math.sqrt(scale), -.5*math.sqrt(scale)],
                                 [[scale, -scale], [-scale, scale]])
        assert q['covariance_rank'] == 1 and abs(q['chi2']-.25) < 1e-12
        assert q['null_projection_norm'] < 1e-12*math.sqrt(scale)
with tempfile.TemporaryDirectory() as directory:
    left, right = Path(directory)/'left', Path(directory)/'right'
    for path, photons in [(left, [2000., 4000.]), (right, [3000., 7000.])]:
        path.mkdir()
        with (path/'iteration_rates.csv').open('w') as stream:
            writer = csv.writer(stream)
            for it, value in enumerate(photons):
                for label, x in [('sp_pt01_match', 10000.), ('dressed_matched', 10000.),
                                 ('sp_pt01_photos', value), ('photos_post', value)]:
                    writer.writerow([1, 1, 1, it, 100, label, x, 1])
    # A PHOTOS-cutoff test must measure PHOTOS, not the invariant analytic column.
    assert chunk_difference(left, right)['dressed_matched'] == [0., 0.]
    assert chunk_difference(left, right, 'photos')['photos_post'] == [2., 1.]
print('Paired ratios, signed/weak denominators, shape covariance and independent pooling passed')
