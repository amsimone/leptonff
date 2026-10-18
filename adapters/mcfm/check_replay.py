#!/usr/bin/env python3
"""Compare offline production replay with the online MCFM decay covariance."""
import argparse
import json
from pathlib import Path
import sys

from analyse_spectra import analyse, definitions


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('run', type=Path)
    p.add_argument('replay', type=Path)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    host = analyse(a.run)
    replay = json.loads(a.replay.read_text())
    meta = json.loads((a.run/'manifest.json').read_text())
    if meta['environment'].get('LEPTONFF_PHOTONIC_NLL') != '1':
        raise RuntimeError('This reference check expects photonic NLL on both paths')
    if (replay['physics_id'] != meta['production_archive']['physics_id'] or
        replay['configuration']['hard_samples'] != int(meta['environment']['LEPTONFF_HARD_SAMPLES']) or
        replay['configuration']['grid_bins'] != int(meta['environment']['LEPTONFF_GRID_BINS'])):
        raise RuntimeError('Replay inputs do not match the online validation run')
    pairs = []
    components = ['born', 'order1', 'matched', 'order2', 'higher']
    rates = ['ff_born', 'decay_order1', 'dressed_matched', 'ff_order2', 'dressed_minus_order1']
    suffixes = ['born', 'nlo', 'match', 'ff2', 'higher']
    for component, rate in zip(components, rates):
        pairs.append((replay['labels'].index('central/photonic_nll/total/'+component), host['labels'].index(rate)))
    global_bin = 0
    for obs, edges in definitions(meta['process']).items():
        for k in range(len(edges)-1):
            for component, suffix in zip(components, suffixes):
                label = f'central/photonic_nll/{obs}:{global_bin}/{component}'
                pairs.append((replay['labels'].index(label), host['labels'].index(f'sp_{obs}{k+1:02d}_{suffix}')))
            global_bin += 1
    differences = [abs(replay['mean_pb'][i]-host['mean_pb'][j]) for i,j in pairs]
    cov_errors, cov_pass = [], True
    for i,j in pairs:
        for k,l in pairs:
            expected = host['covariance_pb2'][j][l]
            error = abs(replay['covariance_pb2'][i][k]-expected)
            cov_errors.append(error)
            cov_pass &= error <= 1e-7+1e-11*abs(expected)
    report = dict(process=meta['process'], channels=len(pairs),
                  maximum_mean_difference_pb=max(differences),
                  maximum_covariance_difference_pb2=max(cov_errors),
                  passed=max(differences)<1e-7 and cov_pass,
                  scope='all analytic rates/bins and full covariance; common-iteration estimator, independent implementations')
    a.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))
    if not report['passed']:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
