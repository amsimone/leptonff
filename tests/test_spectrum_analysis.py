"""Independent event-vector fixture: signs, zeros, migration and covariance."""
import csv,json,math,sys,tempfile
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'adapters/mcfm'))
from analyse_spectra import analyse,definitions,COMPONENTS,RATES

with tempfile.TemporaryDirectory() as directory:
    path=Path(directory);edges=definitions('wplus')
    labels=RATES+[f'sp_{o}{k+1:02d}_{c}' for o,e in edges.items() for k in range(len(e)-1) for c in COMPONENTS]
    groups=[];covariances=[];calls=[4,8]
    for iteration in range(2):
        events=[]
        for event,weight in enumerate([1.2,-.35,.8,0]):
            row=dict.fromkeys(labels,0.)
            for observable,values in edges.items():
                count=len(values)-1;born=(event+iteration)%count;post=(born+1)%count
                for k in range(count):
                    b=weight*(k==born);nlo=b+.1*weight*((k==post)-(k==born))
                    match=nlo+.02*weight*(k==post);psborn=b;photos=b+.08*weight*((k==post)-(k==born))
                    vals=[b,nlo,match,.03*weight*(k==post),match-nlo,psborn,photos,photos-match]
                    for c,x in zip(COMPONENTS,vals):row[f'sp_{observable}{k+1:02d}_{c}']=x
                for c,rate in zip(COMPONENTS,RATES):
                    row[rate]=sum(row[f'sp_{observable}{k+1:02d}_{c}'] for k in range(count))
            events.append(row)
        sums={label:sum(e[label] for e in events) for label in labels}
        covariance={(l,r):(calls[iteration]*sum(e[l]*e[r] for e in events)-sums[l]*sums[r])/(calls[iteration]-1) for i,l in enumerate(labels) for r in labels[i:]}
        groups.append(sums);covariances.append(covariance)
    with (path/'iteration_rates.csv').open('w') as file:
        writer=csv.writer(file)
        for it,g in enumerate(groups):
            for label,x in g.items():writer.writerow([1+it,1,1,it+1,calls[it],label,x,covariances[it][label,label]])
    with (path/'event_covariance.csv').open('w') as file:
        writer=csv.writer(file)
        for it,cov in enumerate(covariances):
            for (left,right),v in cov.items():writer.writerow([1+it,1,1,it+1,calls[it],left,right,v])
    native={}
    for label in labels:
        nonzero=[(g[label],covariances[it][label,label]) for it,g in enumerate(groups) if g[label]!=0]
        inv=sum(1/v for x,v in nonzero)
        native[label]=[[0,1,sum(x/v for x,v in nonzero)/inv if inv else 0,math.sqrt(1/inv) if inv else 0]]
    meta={'status':'completed','process':'wplus','environment':{'LEPTONFF_SPECTRA':'1'},
          'histograms_fb':{k:v for k,v in native.items() if not k.startswith('sp_')},
          'spectra_histograms_fb':{k:v for k,v in native.items() if k.startswith('sp_')}}
    (path/'manifest.json').write_text(json.dumps(meta))
    result=analyse(path)
    for i,label in enumerate(labels):
        expected=sum(g[label] for g in groups)/2*1e-3
        assert abs(result['mean_pb'][i]-expected)<1e-15
        variance=sum(c[label,label] for c in covariances)/4*1e-6
        assert abs(result['covariance_pb2'][i][i]-variance)<1e-18
    for obs,rows in result['spectra'].items():
        assert abs(sum(r['normalised_fractions']['match'][0] for r in rows)-1)<1e-14
    assert result['maximum_covariance_bin_rate_residual_pb2']<1e-18
    # A paired difference must use the correlation retained by the event
    # products, rather than quadrature of independently weighted histograms.
    i=labels.index('photos_minus_coupled');p=labels.index('photos_post');m=labels.index('coupled_matched');cov=result['covariance_pb2']
    assert abs(cov[i][i]-(cov[p][p]+cov[m][m]-2*cov[p][m]))<1e-18
    assert cov[i][i] < cov[p][p]+cov[m][m]
    content=(path/'event_covariance.csv').read_text().splitlines()
    (path/'event_covariance.csv').write_text('\n'.join(content[:-1])+'\n')
    try:analyse(path)
    except RuntimeError:pass
    else:raise AssertionError('Incomplete covariance was accepted')
print('Signed events, resumed stages, variable calls, migration and paired covariance passed')
