#!/usr/bin/env python3
"""Independent flavour/threshold check: source moments and exponential midpoint.

No compiled kernels are used to form the reference. Scalar running uses a
bracketed implicit solution; evolution uses matrix exponentials and Richardson
extrapolation, unlike the production joint RK4 solver.
"""
import argparse,hashlib,json,math,subprocess
from pathlib import Path
from functools import lru_cache
from check_nlo_kernel_reference import integrate_rows,harmonic,ns_moment,singlet_moment
from check_moment_reference import analytic_inputs

def mul(a,b):
    return [[sum(x*y for x,y in zip(row,column)) for column in zip(*b)] for row in a]
def identity(n):return [[float(i==j) for j in range(n)] for i in range(n)]
def norm(m):return max(abs(x) for row in m for x in row)
def distance(a,b):return max(abs(x-y) for row,other in zip(a,b) for x,y in zip(row,other))
def expm(a):
    n=len(a);result=identity(n);term=identity(n)
    for k in range(1,100):
        term=[[x/k for x in row] for row in mul(a,term)]
        result=[[x+y for x,y in zip(row,other)] for row,other in zip(result,term)]
        if norm(term)<1e-18:return result
    raise RuntimeError('Independent exponential did not converge')

@lru_cache(None)
def sources(N):
    p0,p1,d1=analytic_inputs(N);r=integrate_rows(80,0.,N)
    p=1/N-2/(N+1)+2/(N+2)
    logz=-1/N**2+2/(N+1)**2-2/(N+2)**2
    logs=-harmonic(N)/N+2*harmonic(N+1)/(N+1)-2*harmonic(N+2)/(N+2)
    vacuum=-4/(3*N)-8*p/9-4*(logz+logs)/3
    return p0,d1,r,vacuum

def inputs(case,active):
    N=case['N'];fermions=case['fermions'];n=2*len(fermions)+1;g=n-1
    lo,boundary,raw,vacuum=sources(N)
    q2=[f['charge']**2 for f in fermions];col=[f['colours'] for f in fermions]
    s2=sum(c*q for c,q,on in zip(col,q2,active) if on);s4=sum(c*q*q for c,q,on in zip(col,q2,active) if on)
    p0=[[0.]*n for _ in range(n)];p1=[[0.]*n for _ in range(n)]
    p0[g][g]=-2*s2/3;p1[g][g]=s4*raw['photon_to_photon']
    for f,on in enumerate(active):
        if not on:continue
        for i in [2*f,2*f+1]:
            p0[i][i]=q2[f]*lo[0][0];p0[i][g]=q2[f]*lo[0][2];p0[g][i]=col[f]*q2[f]*lo[2][0]
            p1[i][g]=q2[f]**2*raw['lepton_to_photon']
            p1[g][i]=col[f]*(q2[f]**2*(raw['photon_to_lepton']-vacuum)+q2[f]*s2*vacuum)
            for h,other in enumerate(active):
                if not other:continue
                for j in [2*h,2*h+1]:
                    p1[i][j]=col[h]*q2[f]*q2[h]*singlet_moment(N)
                    if i==j:
                        p1[i][j]+=q2[f]**2*(raw['photonic']+raw['same_charge_interference'])+q2[f]*s2*ns_moment(N)
                    elif f==h:p1[i][j]+=q2[f]**2*raw['opposite_charge_interference']
    selected=[2*case['muon_index'],2*case['muon_index']+1,g]
    b=[[0.]*3 for _ in range(n)];mu0=case['mu0']
    for f,on in enumerate(active):
        if not on:continue
        log=2*math.log(mu0/fermions[f]['mass'])
        for i in [2*f,2*f+1]:
            b[i][2]=q2[f]*(boundary[0][2]+log*lo[0][2])
            for c,j in enumerate(selected[:2]):
                if i==j:b[i][c]=q2[f]*(boundary[0][0]+log*lo[0][0])
        b[g][2]-=2*col[f]*q2[f]*log/3
        if f==case['muon_index']:
            b[g][0]=b[g][1]=col[f]*q2[f]*log*lo[2][0]
    return p0,p1,b,2*s2/3,s4

def running(a,L,b,c):
    if a==0:return 0.
    if c==0:return a/(1-b*a*L)
    # Solve the exact implicit equation using a bracket, not Newton in 1/a.
    def elapsed(ratio):
        return (ratio-1)/(b*a*ratio)+c/b**2*(math.log1p(c*a*(ratio-1)/(b+c*a))-math.log(ratio))
    lo,hi=1.,1/(2*math.pi*a)
    for _ in range(70):
        mid=(lo+hi)/2
        if elapsed(mid)<L:lo=mid
        else:hi=mid
    return a*(lo+hi)/2

def reference(case,steps,alpha=None):
    f=case['fermions'];n=2*len(f)+1;g=n-1;selected=[2*case['muon_index'],2*case['muon_index']+1,g]
    low=case['mu0'];active=[x['threshold']<=low for x in f]
    p0,p1,b,_,_=inputs(case,active)
    a=(case['alpha_msbar_at_mu0'] if alpha is None else alpha)/(2*math.pi)
    e=[[float(i==j) for j in selected] for i in range(n)]
    d=[[x+a*y for x,y in zip(row,other)] for row,other in zip(e,b)]
    edges=sorted({x['threshold'] for x in f if low<x['threshold']<=case['mu_frag']}|{case['mu_frag']})
    for high in edges:
        active=[x['threshold']<=math.sqrt(low*high) for x in f]
        p0,p1,_,b0,b1=inputs(case,active);b1=b1 if case['coupling']=='two_loop' else 0
        L=2*math.log(high/low);dt=L/steps
        for k in range(steps):
            am=running(a,(k+.5)*dt,b0,b1)
            u=expm([[dt*(am*x+am*am*y) for x,y in zip(row,other)] for row,other in zip(p0,p1)])
            d=mul(u,d);e=mul(u,e)
        a=running(a,L,b0,b1)
        new=[i for i,x in enumerate(f) if x['threshold']==high and high>case['mu0']]
        oldd=d[g][:];olde=e[g][:];log_sum=0
        for h in new:
            q=f[h]['charge']**2;log=2*math.log(high/f[h]['mass']);N=case['N']
            kernel=q*((N*N+N+2)/(N*(N*N-1))*(log-1)+4/(N-1)**2-4/N**2+2/(N+1)**2)
            for i in [2*h,2*h+1]:d[i]=[a*kernel*x for x in oldd];e[i]=[a*kernel*x for x in olde]
            log_sum+=2*f[h]['colours']*q*log/3
        d[g]=[(1-a*log_sum)*x for x in oldd];e[g]=[(1-a*log_sum)*x for x in olde]
        a=a/(1-a*log_sum);low=high
    return d,e,a

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('executable',type=Path);p.add_argument('--output',type=Path);p.add_argument('--bins',type=int,default=256);a=p.parse_args()
    result=json.loads(subprocess.run([str(a.executable.resolve()),str(a.bins)],check=True,text=True,capture_output=True).stdout)
    maximum=0.;kernel_error=0.;refinement=0.;rows=[]
    for case in result['cases']:
        p0,p1,_,_,_=inputs(case,[True]*len(case['fermions']))
        kernel_error=max(kernel_error,distance(p0,case['final_kernel0']),distance(p1,case['final_kernel1']))
        coarse=reference(case,32);fine=reference(case,64)
        extrapolated=[[[x+(x-y)/3 for x,y in zip(row,other)] for row,other in zip(fine[k],coarse[k])] for k in [0,1]]
        shift=max(distance(fine[k],coarse[k])/3 for k in [0,1]);refinement=max(refinement,shift)
        error=max(distance(extrapolated[0],case['fragmentation']),distance(extrapolated[1],case['evolution']),abs(fine[2]-case['a_final']))
        maximum=max(maximum,error)
        if 'lattice_moment2' in case:
            if distance(case['lattice_moment2'],case['continuum_moment2'])>5e-9:raise RuntimeError('Lattice N2 failure')
            if distance(case['lattice_moment3'],case['fragmentation'])>5e-9:raise RuntimeError('Lattice N3 failure')
        rows.append({'profile':case['profile'],'thresholds':[f['threshold'] for f in case['fermions']],
                     'coupling':case['coupling'],'error':error,'midpoint_refinement':shift})
    if maximum>3e-9 or kernel_error>3e-6:raise RuntimeError(f'Independent flavour reference mismatch: {maximum=}, {kernel_error=}')
    report={'passed':True,'cases':rows,'maximum_evolution_difference':maximum,'maximum_source_kernel_difference':kernel_error,
            'maximum_midpoint_refinement':refinement,'cpp':result,
            'executable_sha256':hashlib.sha256(a.executable.read_bytes()).hexdigest(),
            'scope':'Independent analytic charge factors, source/continuation moments, implicit coupling, finite time-like maps and exponential-midpoint evolution; test quarks are pure-QED synthetic content'}
    if a.output:a.output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['cpp','cases']},indent=2))
if __name__=='__main__':main()
