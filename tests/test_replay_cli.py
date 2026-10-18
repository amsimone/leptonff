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
print('Five models replay common grouped records with correct paired covariance; synthetic fixture only')
