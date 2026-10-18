#!/usr/bin/env python3
"""Export observed bare-muon spectra and paired PHOTOS differences."""
import argparse,json,os
from pathlib import Path
os.environ.setdefault('XDG_CACHE_HOME','/private/tmp/leptonff-cache')
os.environ.setdefault('MPLCONFIGDIR','/private/tmp/leptonff-matplotlib')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('report',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    report=json.loads(a.report.read_text());fig,axes=plt.subplots(2,3,figsize=(13,6.5),gridspec_kw={'height_ratios':[2.3,1]},sharex='col')
    hadronic=report.get('flavour_matching',{}).get('profile')=='hadronic-pointlike'
    qcd_off=(report.get('hadronic_input') or {}).get('qcd_evolution_off',False)
    match_label=('QED + HVP (pointlike)' if qcd_off else 'QED + LO QCD (pointlike)') if hadronic else 'photonic NLO + LL' if report.get('model')=='photonic' else 'coupled NLL match'
    titles={'pt':r'$p_T(\mu)$ [GeV]','mt':r'$m_T(\mu,\nu)$ [GeV]','abseta':r'$|\eta(\mu)|$',
            'mass':r'$m(\mu\mu)$ [GeV]','minpt':r'$\min(p_T(\mu_1),p_T(\mu_2))$ [GeV]',
            'maxeta':r'$\max(|\eta(\mu_1)|,|\eta(\mu_2)|)$'}
    for column,(observable,rows) in enumerate(report['spectra'].items()):
        finite=[row for row in rows if row['upper'] is not None]
        x=np.array([(r['lower']+r['upper'])/2 for r in finite]);dx=np.array([(r['upper']-r['lower'])/2 for r in finite])
        for component,label,color in [('born','Born','0.6'),('nlo','continuum order 1','C2'),('match',match_label,'C0'),('photos','PHOTOS','C1')]:
            values=np.array([r['densities_pb'][component] for r in finite])
            axes[0,column].errorbar(x,values[:,0],yerr=values[:,1],xerr=dx,fmt='.-',ms=3,lw=1,color=color,label=label)
        # Paired difference uncertainty comes from its own grouped-event
        # covariance, never subtraction of native histogram error bars.
        delta=np.array([r['densities_pb']['psdiff'] for r in finite])
        axes[1,column].errorbar(x,delta[:,0],yerr=delta[:,1],xerr=dx,fmt='o',ms=3,color='C3')
        axes[1,column].axhline(0,color='0.5',lw=.7);axes[1,column].set_xlabel(titles[observable])
        axes[0,column].set_ylabel('pb / unit');axes[1,column].set_ylabel('PHOTOS - match\npb / unit')
        for ax in axes[:,column]:ax.grid(alpha=.2)
        if len(finite)!=len(rows):
            overflow=rows[-1]['integrals_pb']['match']
            axes[0,column].text(.97,.96,f"Overflow > {rows[-1]['lower']:g}:\n{overflow[0]:.2g} ± {overflow[1]:.2g} pb",
                                ha='right',va='top',transform=axes[0,column].transAxes,fontsize=8)
    axes[0,0].legend(fontsize=8,loc='upper left' if hadronic else 'best')
    process={'wplus':r'$W^+$','wminus':r'$W^-$','z':r'$Z/\gamma^*$'}[report['process']]
    fig.suptitle(process+' at NLO QCD, 13 TeV; observed bare muons',y=.98)
    coupled_label='One-family' if report.get('flavour_matching',{}).get('profile','one-family')=='one-family' else ('HVP + LO QCD, pointlike boundary;' if report.get('flavour_matching',{}).get('profile')=='hadronic-pointlike' else 'Variable-flavour leptonic')
    footer=(('Pointlike boundary; QED evolution, strong FF evolution off;' if qcd_off else 'Pointlike hadronic FF boundary; QED kernels + LO QCD;')+' PHOTOS comparison diagnostic. Errors: integration statistics.' if hadronic else None)
    fig.text(.5,.015,footer or ('Fixed-alpha, no-pair comparison: analytic NLO + LL versus PHOTOS; errors are integration statistics.' if report.get('model')=='photonic' else coupled_label+' MSbar NLL versus fixed-alpha, no-pair PHOTOS: diagnostic comparison; errors are integration statistics.'),ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.04,1,.95));a.output.parent.mkdir(parents=True,exist_ok=True)
    for extension in ['png','svg']:fig.savefig(a.output.with_suffix('.'+extension),dpi=160)
if __name__=='__main__':main()
