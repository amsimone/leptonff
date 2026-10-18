#!/usr/bin/env python3
"""Independent LO-QCD analytic moments, ordered exponential and HVP anchor."""
import argparse,hashlib,json,math,subprocess
from pathlib import Path
from check_flavour_reference import inputs,running,mul,distance,harmonic

def action(k,d):
    out=[r[:] for r in d];term=[r[:] for r in d]
    for order in range(1,100):
        term=[[x/order for x in row] for row in mul(k,term)]
        out=[[x+y for x,y in zip(a,b)] for a,b in zip(out,term)]
        if max(abs(x) for row in term for x in row)<1e-18:return out
    raise RuntimeError('Independent exponential action did not converge')

def qcd(case,active):
    N=case['N'];nf=sum(on and f['colours']==3 for on,f in zip(active,case['fermions']))
    n=2*len(active)+2;g=n-1;p=[[0.]*n for _ in range(n)]
    qq=4/3*(1.5-2*harmonic(N)+1/N-1/(N+1))
    qg=4/3*(2/(N-1)-2/N+1/(N+1))
    gq=.5*(1/N-2/(N+1)+2/(N+2))
    p[g][g]=5.5-nf/3+6*(-harmonic(N)+1/(N-1)-1/N+1/(N+1)-1/(N+2))
    for f,(on,x) in enumerate(zip(active,case['fermions'])):
        if not on or x['colours']!=3:continue
        for i in [2*f,2*f+1]:p[i][i]=qq;p[i][g]=qg;p[g][i]=gq
    return p,5.5-nf/3

def strong(case,scale):
    if scale<case['start_scale']:return 0.
    ref=case['reference_scale'];b=case['alpha_s_reference']/(2*math.pi)
    edges=sorted({f['threshold'] for f in case['fermions'] if f['colours']==3 and min(scale,ref)<f['threshold']<max(scale,ref)},reverse=scale<ref)+[scale]
    for edge in edges:
        nf=sum(f['colours']==3 and f['threshold']<=math.sqrt(ref*edge) for f in case['fermions'])
        b/=1+b*(5.5-nf/3)*2*math.log(edge/ref);ref=edge
    return b

def reference(case,steps):
    base={**case,'fermions':[{**f,'mass':f['threshold'] if f['pointlike'] else f['mass']} for f in case['fermions']]}
    f=case['fermions'];n=2*len(f)+2;gamma=n-2;selected=[2*case['muon_index'],2*case['muon_index']+1,gamma]
    low=case['mu0'];a=case['alpha_msbar_at_mu0']/(2*math.pi)
    _,_,boundary,_,_=inputs(base,[x['threshold']<=low for x in f]);boundary.append([0.,0.,0.])
    e=[[float(i==j) for j in selected] for i in range(n)]
    d=[[x+a*y for x,y in zip(row,other)] for row,other in zip(e,boundary)]
    for high in sorted({x['threshold'] for x in f if low<x['threshold']<=case['mu_frag']}|{case['mu_frag'],case['hvp_scale']}):
        active=[x['threshold']<=math.sqrt(low*high) for x in f]
        p0,p1,_,b0,b1=inputs(base,active)
        p0=[r+[0.] for r in p0]+[[0.]*n];p1=[r+[0.] for r in p1]+[[0.]*n]
        ps,bs=qcd(case,active);L=2*math.log(high/low);dt=L/steps
        b=strong(case,math.sqrt(low*high))
        # Transport the analytic midpoint coupling back to the left edge.
        b=b/(1-b*bs*L/2) if b else 0.
        for k in range(steps):
            am=running(a,(k+.5)*dt,b0,b1);bm=b/(1+b*bs*(k+.5)*dt)
            generator=[[dt*(am*p0[i][j]+am*am*p1[i][j]+bm*ps[i][j]) for j in range(n)] for i in range(n)]
            d=action(generator,d);e=action(generator,e)
        a=running(a,L,b0,b1);oldd=d[gamma][:];olde=e[gamma][:];match=case['hvp_coefficient'] if high==case['hvp_scale'] else 0.
        for h,x in enumerate(f):
            if x['threshold']!=high or high<=case['mu0']:continue
            if x['pointlike']:kernel=0.
            else:
                log=2*math.log(high/x['mass']);N=case['N'];q=x['charge']**2
                kernel=q*((N*N+N+2)/(N*(N*N-1))*(log-1)+4/(N-1)**2-4/N**2+2/(N+1)**2)
                match+=2*x['colours']*q*log/3
            for i in [2*h,2*h+1]:d[i]=[a*kernel*v for v in oldd];e[i]=[a*kernel*v for v in olde]
        d[gamma]=[(1-a*match)*v for v in oldd];e[gamma]=[(1-a*match)*v for v in olde]
        a/=1-a*match;low=high
    return d,e,a

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('executable',type=Path);p.add_argument('--output',type=Path);p.add_argument('--bins',type=int,default=128);a=p.parse_args()
    result=json.loads(subprocess.run([str(a.executable.resolve()),str(a.bins)],text=True,capture_output=True,check=True).stdout)
    rows=[];maximum=0.;refine=0.;lattice=0.;anchor_error=0.
    for case in result['cases']:
        coarse=reference(case,32);fine=reference(case,64)
        error=max(distance([[x+(x-y)/3 for x,y in zip(r,s)] for r,s in zip(fine[k],coarse[k])],case['fragmentation' if k==0 else 'evolution']) for k in [0,1])
        maximum=max(maximum,error);refine=max(refine,max(distance(fine[k],coarse[k])/3 for k in [0,1]))
        if 'lattice_moment' in case:lattice=max(lattice,distance(case['lattice_moment'],case['fragmentation']))
        alpha=1/137.035999084;pi=math.pi
        am=alpha/(1-alpha*2*math.log(.1056583755/.00051099895069)/(3*pi))/(2*pi)
        am=running(am,2*math.log(1.77693/.1056583755),4/3,2)
        am=running(am,2*math.log(91.1876/1.77693),2,3)
        shift=.02760+1-alpha/(2*pi*am)
        shift+=alpha/pi*(55/27+.1193/pi*(605/108-44*1.2020569031595942854/9))
        anchor_error=max(anchor_error,abs(case['a_final']*2*pi-alpha/(1-shift)))
        rows.append({'N':case['N'],'alpha_s_reference':case['alpha_s_reference'],'maximum_difference':error})
    if maximum>3e-8 or lattice>3e-8 or anchor_error>1e-13:raise RuntimeError(f'{maximum=}, {lattice=}, {anchor_error=}')
    report={'passed':True,'cases':rows,'maximum_evolution_difference':maximum,'maximum_midpoint_refinement':refine,
        'maximum_lattice_moment_difference':lattice,'anchor_error':anchor_error,'executable_sha256':hashlib.sha256(a.executable.read_bytes()).hexdigest(),
        'scope':'Independent analytic LO QCD moments and exponential action; sourced QED moments; pointlike hadronic boundary with fixed HVP anchor; not mixed-kernel accuracy', 'cpp':result}
    if a.output:a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='cpp'},indent=2))
if __name__=='__main__':main()
