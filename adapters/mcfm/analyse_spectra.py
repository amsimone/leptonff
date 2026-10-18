#!/usr/bin/env python3
"""Observed spectra with common iteration weights and grouped-event covariance."""
import argparse,csv,hashlib,json,math
from collections import defaultdict
from pathlib import Path

COMPONENTS=['born','nlo','match','order2','higher','psborn','photos','psdiff']
RATES=['ff_born','decay_order1','coupled_matched','coupled_order2','coupled_minus_order1',
       'photos_born','photos_post','photos_minus_coupled']

def definitions(process):
    return ({'mass':[60,70,80,85,90,95,100,110,120], 'minpt':[25,30,35,40,50,None],
             'maxeta':[0,.6,1.2,1.8,2.4]} if process=='z' else
            {'pt':[25,30,35,40,45,60,None], 'mt':[0,40,60,70,80,90,120,None],
             'abseta':[0,.6,1.2,1.8,2.4]})

def analyse(path):
    meta=json.loads((path/'manifest.json').read_text())
    if meta['status']!='completed' or meta['environment'].get('LEPTONFF_SPECTRA')!='1':
        raise RuntimeError('A completed physical-spectrum run is required')
    photonic=meta['environment'].get('LEPTONFF_PHOTONIC_LL')=='1'
    photonic_nll=meta['environment'].get('LEPTONFF_PHOTONIC_NLL')=='1'
    flavours=meta.get('flavour_matching',{}).get('profile','one-family')
    components=['ff2' if c=='order2' else c for c in COMPONENTS] if photonic else COMPONENTS
    rates=['ff_born','decay_order1','dressed_matched','ff_order2','dressed_minus_order1',
           'photos_born','photos_post','photos_minus_dressed'] if photonic else RATES
    if 'config_sha256' in meta and hashlib.sha256((path/'input.ini').read_bytes()).hexdigest()!=meta['config_sha256']:
        raise RuntimeError('Run configuration has changed')
    if meta.get('adapter_build') and meta['executable_sha256']!=meta['adapter_build']['adapter_executable_sha256']:
        raise RuntimeError('Executable and build manifest disagree')
    edges=definitions(meta['process'])
    bins=[(obs,k,values[k],values[k+1]) for obs,values in edges.items() for k in range(len(values)-1)]
    labels=rates+[f'sp_{obs}{k+1:02d}_{c}' for obs,k,_,_ in bins for c in components]
    index={key:i for i,key in enumerate(labels)};n=len(labels)
    groups=defaultdict(dict)
    with (path/'iteration_rates.csv').open() as file:
        for row in csv.reader(file):
            if row[5] not in index:continue
            stage,part,ips,iteration,calls=map(int,row[:5]);key=(part,ips,iteration)
            if stage<1:raise RuntimeError('Warm-up iterations are not production data')
            if row[5] in groups[key]:raise RuntimeError('Duplicate raw iteration')
            x,v=map(float,row[6:]);
            if not all(map(math.isfinite,[x,v])) or v<0:raise RuntimeError('Invalid raw moments')
            groups[key][row[5]]=(x,v,calls)
    if not groups or any(set(g)!=set(labels) for g in groups.values()):raise RuntimeError('Incomplete raw spectrum')
    counts=defaultdict(int);callcounts=defaultdict(set)
    for (part,ips,it),g in groups.items():
        counts[(part,ips)]+=1;callcounts[(part,ips)].add(g[labels[0]][2])
        if len({r[2] for r in g.values()})!=1:raise RuntimeError('Inconsistent calls within an iteration')
    if any(c<2 for c in counts.values()):
        raise RuntimeError('Common-weight audit requires at least two iterations per part')
    means=[0.]*n;covariance=[[0.]*n for _ in labels];chunk_cov=[[0.]*n for _ in labels]
    maxclosure=0.
    for key,g in groups.items():
        for obs,values in edges.items():
            for c,rate in zip(components,rates):
                total=sum(g[f'sp_{obs}{k+1:02d}_{c}'][0] for k in range(len(values)-1))
                maxclosure=max(maxclosure,abs(total-g[rate][0])*1e-3)
        for obs,k,_,_ in bins:
            v=lambda c:g[f'sp_{obs}{k+1:02d}_{c}'][0]
            maxclosure=max(maxclosure,abs(v('higher')-v('match')+v('nlo'))*1e-3,
                           abs(v('psdiff')-v('photos')+v('match'))*1e-3)
        weight=1/counts[key[:2]]
        for label,i in index.items():means[i]+=g[label][0]*weight*1e-3
    if maxclosure>1e-7:raise RuntimeError('Physical spectrum algebra/coverage failure')
    # Cross-products were formed AFTER summing every real/dipole counterevent.
    seen=set();maxdiag=0.
    with (path/'event_covariance.csv').open() as file:
        for row in csv.reader(file):
            stage,part,ips,it,calls=map(int,row[:5]);left,right=row[5:7];v=float(row[7])
            key=(part,ips,it);i,j=index[left],index[right];pair=(key,min(i,j),max(i,j))
            if stage<1 or key not in groups or calls!=groups[key][left][2]:raise RuntimeError('Unmatched covariance iteration')
            if pair in seen or not math.isfinite(v):raise RuntimeError('Duplicate/nonfinite event covariance')
            seen.add(pair);scale=1e-6/counts[key[:2]]**2
            covariance[i][j]+=v*scale
            if i!=j:covariance[j][i]+=v*scale
            else:maxdiag=max(maxdiag,abs(v-groups[key][left][1])*1e-6)
    if len(seen)!=len(groups)*n*(n+1)//2 or maxdiag>1e-5:
        raise RuntimeError('Incomplete covariance or disagreement with host variances')
    for part,count in counts.items():
        data=[g for key,g in groups.items() if key[:2]==part]
        avg=[sum(g[label][0] for g in data)/count for label in labels]
        for i,left in enumerate(labels):
            for j,right in enumerate(labels):
                chunk_cov[i][j]+=sum((g[left][0]-avg[i])*(g[right][0]-avg[j]) for g in data)/(count*(count-1))*1e-6
    # The full event covariance must reproduce total-rate covariance when
    # summing ANY spectrum, including cross-channel PHOTOS/analytic products.
    maxcovclosure=0.
    for obs,values in edges.items():
        for c,rate in zip(components,rates):
            indices=[index[f'sp_{obs}{k+1:02d}_{c}'] for k in range(len(values)-1)]
            for j in range(n):
                summed=sum(covariance[i][j] for i in indices)
                maxcovclosure=max(maxcovclosure,abs(summed-covariance[index[rate]][j]))
    scale=max(1.,max(abs(x) for row in covariance for x in row))
    if maxcovclosure>2e-8*scale:raise RuntimeError('Covariance bin/rate closure failure')
    # Native histogram estimates have independent inverse-variance weights;
    # reconstruct them but use common weights for every physical comparison.
    maxnative=0.
    for label in labels:
        native=[0.,0.]
        for part,count in counts.items():
            nonzero=[g[label] for key,g in groups.items() if key[:2]==part and g[label][0]!=0]
            if nonzero:
                if any(v<=0 for x,v,calls in nonzero):raise RuntimeError('Nonzero native estimate without variance')
                inv=sum(1/v for x,v,calls in nonzero)
                native[0]+=sum(x/v for x,v,calls in nonzero)/inv*1e-3;native[1]+=1/inv*1e-6
        actual=(meta['spectra_histograms_fb'] if label.startswith('sp_') else meta['histograms_fb'])[label][0]
        maxnative=max(maxnative,abs(native[0]-actual[2]*1e-3),abs(math.sqrt(native[1])-actual[3]*1e-3))
    if maxnative>1e-3:raise RuntimeError('Native histogram reconstruction failed')
    def linear_error(coefficients):
        var=sum(x*y*covariance[i][j] for i,x in coefficients for j,y in coefficients)
        if var < -1e-8*scale:raise RuntimeError('Negative propagated variance')
        return math.sqrt(max(0.,var))
    results={}
    for obs,values in edges.items():
        rows=[]
        for k in range(len(values)-1):
            row={'lower':values[k],'upper':values[k+1],'integrals_pb':{},'normalised_fractions':{},'densities_pb':{}}
            width=None if values[k+1] is None else values[k+1]-values[k]
            for c,rate in zip(components,rates):
                i=index[f'sp_{obs}{k+1:02d}_{c}'];error=linear_error([(i,1)])
                row['integrals_pb'][c]=[means[i],error]
                row['densities_pb'][c]=[means[i]/width,error/width] if width else None
                d=index[rate]
                if c in ['born','nlo','match','psborn','photos'] and means[d]!=0:
                    row['normalised_fractions'][c]=[means[i]/means[d],
                        linear_error([(i,1/means[d]),(d,-means[i]/means[d]**2)])]
            rows.append(row)
        results[obs]=rows
    normalised_covariance={}
    for obs,values in edges.items():
        normalised_covariance[obs]={}
        for component,rate in zip(components,rates):
            if component not in ['born','nlo','match','psborn','photos']:continue
            d=index[rate];denominator=means[d]
            if denominator==0:continue
            names=[f'sp_{obs}{k+1:02d}_{component}' for k in range(len(values)-1)]
            indices=[index[name] for name in names]
            matrix=[[(covariance[i][j]-means[j]/denominator*covariance[i][d]
                       -means[i]/denominator*covariance[d][j]
                       +means[i]*means[j]/denominator**2*covariance[d][d])/denominator**2
                     for j in indices] for i in indices]
            magnitude=max(1e-20,max(abs(v) for row in matrix for v in row))
            if max(abs(sum(row)) for row in matrix)>1e-7*magnitude:
                raise RuntimeError('Normalised covariance violates the sum constraint')
            for k,row in enumerate(matrix):
                error=results[obs][k]['normalised_fractions'][component][1]
                if abs(row[k]-error**2)>1e-7*magnitude:
                    raise RuntimeError('Normalised covariance diagonal disagrees with propagated errors')
            normalised_covariance[obs][component]={'labels':names,'covariance':matrix}
    return {'analysis_schema_version':3,'analysis_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'process':meta['process'],'integration_warnings':meta.get('integration_warnings',[]),'model':('photonic_nll' if photonic_nll else 'photonic' if photonic else 'coupled'),'spectra':results,'labels':labels,'mean_pb':means,
            'covariance_pb2':covariance,'normalised_shape_covariance':normalised_covariance,'chunk_covariance_diagnostic_pb2':chunk_cov,
            'iterations_per_contribution':{str(k):v for k,v in counts.items()},
            'calls_per_contribution':{str(k):sorted(v) for k,v in callcounts.items()},
            'maximum_bin_rate_or_overlap_residual_pb':maxclosure,
            'maximum_event_covariance_diagonal_residual_pb2':maxdiag,
            'maximum_covariance_bin_rate_residual_pb2':maxcovclosure,
            'maximum_native_histogram_reconstruction_residual_pb':maxnative,
            'provenance':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                          [path/'manifest.json',path/'iteration_rates.csv',path/'event_covariance.csv']},
            'flavour_matching':meta.get('flavour_matching',{'profile':'one-family'}),
            'resolved_flavour_coupling':meta.get('resolved_flavour_coupling'),
            'hadronic_input':meta.get('hadronic_input'),
            'coefficient2_scope':'FF-only photonic NLL boundary/kernel coefficient' if photonic_nll else 'FF-only LL coefficient' if photonic else 'Full retained coupled matching coefficient',
            'comparison_scope':'Aligned fixed-alpha/no-pair content; unitary-conditioned photonic NLL' if photonic_nll else 'Aligned fixed-alpha/no-pair content; LL with one finite insertion, no P1' if photonic else ('One-family MSbar NLL vs fixed-alpha/no-pair PHOTOS: diagnostic only' if flavours=='one-family' else (('Pure-QED with HVP anchor and pointlike boundary; strong FF evolution off' if (meta.get('hadronic_input') or {}).get('qcd_evolution_off') else 'HVP anchor and LO QCD with pointlike hadronic FF boundary')+' vs fixed-alpha/no-pair PHOTOS: diagnostic only' if flavours=='hadronic-pointlike' else 'Variable-flavour pure-QED MSbar NLL vs fixed-alpha/no-pair PHOTOS: diagnostic only')),
            'scope':'Equal weights per production iteration within each contribution; each iteration retains its own call count; full cross-bin/channel grouped-event covariance; chunk covariance is a tail diagnostic, not a precision bound'}

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('run',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    report=analyse(a.run);a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k.startswith('maximum_')},indent=2))
if __name__=='__main__':main()
