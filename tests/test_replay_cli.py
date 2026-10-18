"""Synthetic grouped archive: all model routes and paired covariance."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

engine=Path(sys.argv[1]).resolve()
with tempfile.TemporaryDirectory(prefix='leptonff-replay-') as directory:
    base=Path(directory);archive=base/'fixture.lff'
    lines=['LEPTONFF_PRODUCTION 1','META 1 "synthetic-fixture" "test-physics" "test-run" "fb" 1045 1 0 1','WEIGHTS "central"']
    for iteration in [1,2]:
        scale=1+.1*iteration
        lines += ['GROUP 1 2',f'POINT 123 50 40 0 30 50 -40 0 -30 {scale}',
                  f'POINT 123 50 40 0 30 50 -40 0 -30 {-scale/4}',
                  'GROUP 2 1',f'POINT 456 50 30 40 0 50 -30 -40 0 {scale/2}',f'BLOCK 1 1 {iteration} 4 2']
    lines += ['END'];archive.write_text('\n'.join(lines)+'\n')
    (base/'bins.conf').write_text('pt 25 35 inf\n')
    models=['born','photonic_ll','photonic_nll','leptonic_nll','qed_hvp']
    command=[str(engine),str(archive),'--output',str(base/'result.json'),'--binning',str(base/'bins.conf'),
             '--grid-bins','128','--hard-samples','2']
    for model in models:command += ['--model',model]
    subprocess.run(command,check=True)
    data=json.loads((base/'result.json').read_text())
    parallel=command.copy();parallel[parallel.index(str(base/'result.json'))]=str(base/'parallel.json')
    subprocess.run(parallel+['--threads','2','--batch-groups','1'],check=True)
    threaded=json.loads((base/'parallel.json').read_text())
    assert data['mean_pb']==threaded['mean_pb'] and data['covariance_pb2']==threaded['covariance_pb2']
    parallel[parallel.index(str(base/'parallel.json'))]=str(base/'threaded-batch.json')
    subprocess.run(parallel+['--threads','3','--batch-groups','8'],check=True)
    threaded=json.loads((base/'threaded-batch.json').read_text())
    assert data['mean_pb']==threaded['mean_pb'] and data['covariance_pb2']==threaded['covariance_pb2']
    for layout in ('diagonal','per-weight'):
        compact=command.copy();compact[compact.index(str(base/'result.json'))]=str(base/f'{layout}.json')
        subprocess.run(compact+['--covariance',layout],check=True)
        sparse=json.loads((base/f'{layout}.json').read_text())
        assert sparse['schema']==2 and 'covariance_pb2' not in sparse
        assert sparse['mean_pb']==data['mean_pb']
        for (i,j),value in zip(sparse['covariance_pairs'],sparse['covariance_values_pb2']):
            assert value==data['covariance_pb2'][i][j]
    assert data['models']==models and len(data['bins'])==2
    assert data['groups']==4 and data['terms']==6 and len(data['labels'])==75
    means=data['mean_pb'];cov=data['covariance_pb2']
    for model in models[1:]:
        i=data['labels'].index('central/born/total/born');j=data['labels'].index(f'central/{model}/total/born')
        assert means[i]==means[j]
        assert abs(cov[i][i]+cov[j][j]-2*cov[i][j])<1e-15
    for component in ['born','order1']:
        i=data['labels'].index(f'central/photonic_ll/total/{component}');j=data['labels'].index(f'central/photonic_nll/total/{component}')
        assert means[i]==means[j]
        assert abs(cov[i][i]+cov[j][j]-2*cov[i][j])<1e-15
    # Larger channel vectors can be replayed without allocating a dense matrix.
    large=base/'large.lff'
    weights=['central']+[f'w{i}' for i in range(6)]
    large.write_text('LEPTONFF_PRODUCTION 1\nMETA 31 "fixture" "physics" "run" "fb" 42 1 0 7\nWEIGHTS '
                     +' '.join(f'"{w}"' for w in weights)+'\nBLOCK 1 1 1 2 0\nEND\n')
    command=[str(engine),str(large),'--output',str(base/'large.json'),'--subdivide','4',
             '--grid-bins','128','--hard-samples','2','--covariance','diagonal']
    for model in models:command += ['--model',model]
    subprocess.run(command,check=True)
    large_result=json.loads((base/'large.json').read_text())
    assert len(large_result['mean_pb'])==11550 and all(x==0 for x in large_result['variance_pb2'])
print('Five models replay common grouped records with correct paired covariance; synthetic fixture only')
