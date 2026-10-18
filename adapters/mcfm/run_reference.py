#!/usr/bin/env python3
"""Run a native or FF-only NLO W/Z benchmark with MCFM 10.3.1."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import re
import math
import shutil
import configparser
from datetime import datetime,timezone

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as file:
        for chunk in iter(lambda:file.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def configuration(process,native,output,pdf,precision,alpha_dipole,seed=1045,calls_factor=1,call_scale=1,mu_r_factor=1,mu_f_factor=1,pdf_member=0):
    z=process=='z';proc={'wplus':1,'wminus':6,'z':31}[process]
    sections={
      'extra':{'nohistograms':'.false.'},
      'general':{'nproc':proc,'part':'nlo','runstring':output.name,'rundir':output,
        'sqrts':13000,'ih1':1,'ih2':1,'zerowidth':'.false.' if z else '.true.',
        'removebr':'.false.','ewcorr':'none'},
      'lhapdf':{'lhapdfset':pdf,'lhapdfmember':pdf_member,'dopdferrors':'.false.'},
      'scales':{'renscale':(91.1876 if z else 80.385)*mu_r_factor,'facscale':(91.1876 if z else 80.385)*mu_f_factor,
        'dynamicscale':'none','doscalevar':'.false.'},
      'masses':{'wmass':80.385,'zmass':91.1876},
      'masscuts':{'m34min':60 if z else 0,'m34max':120 if z and native else 13000},
      'cuts':{'makecuts':'.true.' if native else '.false.','ptleptmin':25,'etaleptmax':2.4,
        'ptlept2min':25 if z else 0,'etalept2max':2.4,'ptminmiss':0,'m34transmin':0,
        'Rjlmin':0,'Rllmin':0},
      'histogram':{'writetxt':'.true.','writetop':'.false.','newstyle':'.false.'},
      'dipoles':{'alpha':alpha_dipole},
      'integration':{'seed':seed,'usesobol':'.false.','precisiongoal':precision,
        'initcallslord':round(10000*calls_factor*call_scale),'initcallsnloreal':round(20000*calls_factor*call_scale),'initcallsnlovirt':round(10000*calls_factor*call_scale),
        'iterbatch1':10,'iterbatchwarmup':5,
        'readin':'.false.','writeintermediate':'.true.'}}
    # MCFM's parser REQUIRES indentation to associate keys with a section.
    return 'mcfm_version = 10.3\nwriterefs = .true.\n\n'+''.join(
      '['+name+']\n'+''.join('    '+k+' = '+str(v)+'\n' for k,v in values.items())+'\n'
      for name,values in sections.items())

def read_histogram(path):
    rows=[]
    for line in path.read_text().splitlines():
        if line.strip() and not line.startswith('#'):
            rows.append([float(x) for x in line.split()])
    return rows

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('host',type=Path)
    parser.add_argument('process',choices=['wplus','wminus','z'])
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--adapter',type=Path,default=Path('build-mcfm-adapter/mcfm_leptonff'))
    parser.add_argument('--native',action='store_true')
    parser.add_argument('--export-production',action='store_true',help='Write lossless grouped pre-FSR production records')
    parser.add_argument('--production-only',action='store_true',help='Export production without evaluating fragmentation or PHOTOS')
    parser.add_argument('--bins',type=int,default=4096)
    parser.add_argument('--threads',type=int,default=2)
    parser.add_argument('--alpha',type=float,default=1/137.035999084)
    parser.add_argument('--alpha-scheme',choices=['msbar','on-shell'])
    parser.add_argument('--flavours','--flavors',choices=['one-family','leptons','custom','hadronic-pointlike'],default='one-family')
    parser.add_argument('--flavour-file',type=Path)
    parser.add_argument('--threshold-factor',type=float,default=1)
    parser.add_argument('--dipole-alpha',type=float,default=1)
    photonic_modes=parser.add_mutually_exclusive_group()
    photonic_modes.add_argument('--photonic-ll',action='store_true')
    photonic_modes.add_argument('--photonic-nll',action='store_true',help='Fixed-alpha unitary-conditioned photonic NLL; requires hard dressing')
    parser.add_argument('--hard-samples',type=int,default=0)
    parser.add_argument('--dress-hard',action='store_true')
    parser.add_argument('--coupled-hard',action='store_true')
    parser.add_argument('--fraction-cells',action='store_true')
    parser.add_argument('--spectra',action='store_true',help='Observed W pT/mT/eta and charge-symmetric Z spectra')
    parser.add_argument('--mu-frag',type=float)
    parser.add_argument('--hard-soft',type=float,default=1e-5)
    parser.add_argument('--hard-collinear',type=float,default=1e-5)
    parser.add_argument('--photos-order',choices=['1','exponentiated'],default='exponentiated')
    parser.add_argument('--photos-cutoff',type=float)
    parser.add_argument('--photos-linear-scale',type=float,default=1)
    parser.add_argument('--precision',type=float,default=.02)
    parser.add_argument('--seed',type=int,default=1045)
    parser.add_argument('--calls-factor',type=int,default=1)
    parser.add_argument('--call-scale',type=float,default=1,help='Explicit benchmark call scaling; >=0.01')
    parser.add_argument('--pdf',default='NNPDF31_nnlo_as_0118')
    parser.add_argument('--pdf-member',type=int,default=0)
    parser.add_argument('--mu-r-factor',type=float,default=1)
    parser.add_argument('--mu-f-factor',type=float,default=1)
    parser.add_argument('--hadronic-input',type=Path)
    parser.add_argument('--hadronic-scale',type=float,default=2)
    parser.add_argument('--hadronic-shift',type=float,default=0)
    parser.add_argument('--alpha-s',type=float)
    parser.add_argument('--qcd-off',action='store_true')
    parser.add_argument('--log-span',type=float)
    args=parser.parse_args();root=Path(__file__).resolve().parents[2]
    photonic=args.photonic_ll or args.photonic_nll
    if args.production_only and (not args.export_production or args.native or args.spectra or args.hard_samples or photonic or args.coupled_hard):
        parser.error('--production-only requires --export-production and no decay/spectrum flags')
    if args.export_production and (args.native or args.threads!=1):
        parser.error('--export-production requires the adapter and --threads 1')
    if args.alpha_scheme is None:args.alpha_scheme='on-shell' if args.flavours in ['leptons','hadronic-pointlike'] else 'msbar'
    if args.flavours=='hadronic-pointlike':
        if args.alpha_scheme!='on-shell':parser.error('Hadronic pointlike profile requires on-shell input')
        if not math.isfinite(args.hadronic_scale) or not 1<=args.hadronic_scale<=4:parser.error('Hadronic scale must be 1..4 GeV')
        if not math.isfinite(args.hadronic_shift) or abs(args.hadronic_shift)>5:parser.error('Hadronic shift must be within five input standard deviations')
        if args.hadronic_input is None:args.hadronic_input=root/'configs/hadronic_alpha_pdg2025.dat'
    elif args.hadronic_input or args.qcd_off or args.alpha_s is not None or args.hadronic_shift!=0 or args.hadronic_scale!=2:
        parser.error('Hadronic and strong options require --flavours hadronic-pointlike')
    if args.alpha_s is not None and (not math.isfinite(args.alpha_s) or not 0<=args.alpha_s<=.3):parser.error('Invalid strong coupling input')
    if args.log_span is None:args.log_span=64 if args.flavours=='hadronic-pointlike' else 40
    maximum_span=128 if photonic else 80
    if not math.isfinite(args.log_span) or not 20<=args.log_span<=maximum_span:
        parser.error('Lattice log span must lie between 20 and '+str(maximum_span))
    if args.flavours!='one-family' and (not args.coupled_hard or args.native):
        parser.error('Variable flavours require --coupled-hard with the adapter')
    if args.flavours=='custom' and not args.flavour_file:parser.error('--flavours custom needs --flavour-file')
    if args.flavours!='custom' and args.flavour_file:parser.error('--flavour-file requires --flavours custom')
    if not math.isfinite(args.threshold_factor) or not .25<=args.threshold_factor<=4:
        parser.error('--threshold-factor must lie between .25 and 4')
    if args.flavours!='leptons' and args.threshold_factor!=1:parser.error('Use --threshold-factor only with --flavours leptons')
    if photonic and args.alpha_scheme!='msbar':parser.error('Photonic mode takes fixed alpha without input conversion')
    if args.photonic_nll and (args.native or not args.dress_hard or args.hard_samples<1):
        parser.error('--photonic-nll requires the adapter, --dress-hard and --hard-samples >=1')
    if args.pdf_member<0:parser.error('--pdf-member must be nonnegative')
    if any(not math.isfinite(v) or v<=0 for v in [args.mu_r_factor,args.mu_f_factor]):
        parser.error('QCD scale factors must be finite and positive')
    if not math.isfinite(args.call_scale) or args.call_scale<.01:parser.error('--call-scale must be finite and >=0.01')
    if args.spectra and args.threads!=1:parser.error('--spectra requires --threads 1 for event covariance')
    if args.calls_factor<1:parser.error('--calls-factor must be positive')
    if args.dress_hard and (not photonic or args.hard_samples<1):
        parser.error('--dress-hard requires --photonic-ll/--photonic-nll and --hard-samples >=1')
    if args.coupled_hard and (photonic or args.dress_hard or args.hard_samples<1):
        parser.error('--coupled-hard requires coupled evolution and --hard-samples >=1')
    if args.fraction_cells and not args.coupled_hard:
        parser.error('--fraction-cells requires --coupled-hard')
    if args.spectra and (not (args.coupled_hard or args.dress_hard) or args.native):
        parser.error('--spectra requires --coupled-hard or --photonic-ll --dress-hard')
    host=args.host.resolve();out=args.output.resolve()
    if out.exists():raise RuntimeError('Refusing to overwrite an existing run directory')
    exe=host/'build-gcc/mcfm' if args.native else args.adapter.resolve()
    build_manifest=exe.parent/'manifest.json'
    build_info=json.loads(build_manifest.read_text()) if not args.native and build_manifest.exists() else None
    if build_info and build_info['adapter_executable_sha256']!=digest(exe):
        raise RuntimeError('Adapter executable does not match its build manifest')
    if args.spectra and not (build_info or {}).get('iteration_audit'):
        raise RuntimeError('Spectra require a build with --iteration-audit for grouped-event covariance')
    if args.photonic_nll and 'photonic-nll-unitary-v1' not in (build_info or {}).get('capabilities',[]):
        raise RuntimeError('Photonic NLL requires a rebuilt adapter advertising photonic-nll-unitary-v1')
    if args.export_production and (not (build_info or {}).get('iteration_audit') or
        'production-archive-v1' not in (build_info or {}).get('capabilities',[])):
        raise RuntimeError('Production records require an audited adapter with production-archive-v1 capability')
    out.mkdir(parents=True)
    flavour_file=None
    if args.flavour_file:
        flavour_file=out/'flavours.dat';shutil.copyfile(args.flavour_file.resolve(),flavour_file)
    hadronic_file=None
    if args.hadronic_input:
        hadronic_file=out/'hadronic.dat';shutil.copyfile(args.hadronic_input.resolve(),hadronic_file)
    for name in ('process.DAT','br.sm1','br.sm2','PDFs'):
        (out/name).symlink_to(host/'build-gcc'/name)
    (out/'input.ini').write_text(configuration(args.process,args.native,out,args.pdf,args.precision,args.dipole_alpha,
                                            args.seed,args.calls_factor,args.call_scale,args.mu_r_factor,args.mu_f_factor,args.pdf_member))
    env=os.environ.copy();env.update(OMP_NUM_THREADS=str(args.threads),
      LHAPDF_DATA_PATH=str(host/'Bin/PDFs'),LEPTONFF_GRID_BINS=str(args.bins),LEPTONFF_ALPHA=str(args.alpha),LEPTONFF_LOG_SPAN=str(args.log_span),
      LEPTONFF_HADRONIC_SCALE=str(args.hadronic_scale),LEPTONFF_HADRONIC_SHIFT=str(args.hadronic_shift),LEPTONFF_QCD_OFF='1' if args.qcd_off else '0',
      LEPTONFF_ALPHA_SCHEME=args.alpha_scheme,LEPTONFF_FLAVOURS=args.flavours,
      LEPTONFF_THRESHOLD_FACTOR=str(args.threshold_factor),
      LEPTONFF_PT_MIN='25',LEPTONFF_ETA_MAX='2.4',LEPTONFF_NU_PT_MIN='0',LEPTONFF_MT_MIN='0',
      LEPTONFF_MASS_MIN='60',LEPTONFF_MASS_MAX='120',
      LEPTONFF_PHOTONIC_LL='1' if photonic else '0',
      LEPTONFF_PHOTONIC_NLL='1' if args.photonic_nll else '0',
      LEPTONFF_PRODUCTION_ONLY='1' if args.production_only else '0',
      LEPTONFF_PRODUCTION_ARCHIVE='',
      LEPTONFF_HARD_SAMPLES=str(args.hard_samples),LEPTONFF_HARD_SOFT=str(args.hard_soft),
      LEPTONFF_FRACTION_CELLS='1' if args.fraction_cells else '0',
      LEPTONFF_SPECTRA='1' if args.spectra else '0',
      LEPTONFF_COUPLED_HARD='1' if args.coupled_hard else '0',
      LEPTONFF_DRESS_HARD='1' if args.dress_hard else '0',
      LEPTONFF_HARD_COLLINEAR=str(args.hard_collinear),LEPTONFF_PHOTOS_ORDER=args.photos_order,
      LEPTONFF_MU_FRAG=str(args.mu_frag if args.mu_frag is not None else (91.1876 if args.process=='z' else 80.385)))
    if flavour_file:env['LEPTONFF_FLAVOUR_FILE']=str(flavour_file)
    else:env.pop('LEPTONFF_FLAVOUR_FILE',None)
    env['LEPTONFF_PHOTOS_CUTOFF']=str(args.photos_cutoff if args.photos_cutoff is not None else 1e-5)
    if hadronic_file:env['LEPTONFF_HADRONIC_FILE']=str(hadronic_file)
    else:env.pop('LEPTONFF_HADRONIC_FILE',None)
    if args.alpha_s is not None:env['LEPTONFF_ALPHA_S_REFERENCE']=str(args.alpha_s)
    else:env.pop('LEPTONFF_ALPHA_S_REFERENCE',None)
    env['LEPTONFF_PHOTOS_LINEAR_SCALE']=str(args.photos_linear_scale)
    pdf_hashes={p.name:digest(p) for p in sorted((host/'Bin/PDFs'/args.pdf).glob('*')) if p.is_file()}
    if args.export_production:
        host_config=configparser.ConfigParser(interpolation=None)
        host_config.read_string('[version]\n'+(out/'input.ini').read_text())
        physics={section:dict(host_config[section]) for section in host_config.sections()
                 if section not in ['integration','dipoles','histogram','extra']}
        for key in ['runstring','rundir']:physics['general'].pop(key,None)
        physics.update(pdf_files=pdf_hashes,host_executable=build_info['host_executable_sha256'],
                       branching_inputs={name:digest(out/name) for name in ['br.sm1','br.sm2']})
        physics_id=hashlib.sha256(json.dumps(physics,sort_keys=True).encode()).hexdigest()
        env.update(LEPTONFF_PRODUCTION_ARCHIVE=str(out/'production.lff'),
                   LEPTONFF_ARCHIVE_PROCESS=str({'wplus':1,'wminus':6,'z':31}[args.process]),
                   LEPTONFF_ARCHIVE_PHYSICS_ID=physics_id,LEPTONFF_ARCHIVE_RUN_ID=digest(out/'input.ini'),
                   LEPTONFF_ARCHIVE_RANDOM_SEED=str(args.seed),LEPTONFF_ARCHIVE_FULL_ANGLES='1')
    sys.path.insert(0,str(root/'python'))
    from check_reference import source_digest
    adapter_files=sorted((root/'adapters/mcfm').glob('*'))
    meta={'status':'running','started_utc':datetime.now(timezone.utc).isoformat(),
      'process':args.process,'native':args.native,'command':[str(exe),'input.ini'],
      'cwd':str(out),'source_sha256':build_info['source_sha256'] if build_info else source_digest(root),
      'working_tree_source_sha256':source_digest(root),'executable_sha256':digest(exe),
      'config_sha256':digest(out/'input.ini'),
      'adapter_build':build_info,
      'pdf_sha256':pdf_hashes,
      'adapter_source_sha256':{str(p.relative_to(root)):digest(p) for p in adapter_files if p.is_file()},
      'environment':{k:v for k,v in env.items() if k.startswith('LEPTONFF_') or k in ['OMP_NUM_THREADS','LHAPDF_DATA_PATH']},
      'theory':'NLO QCD host; W pole / off-shell Z-gamma; one-family MSbar FF-only convolution',
      'pdf_note':'Available NNLO PDF set used for host-interface validation; matrix elements are NLO QCD',
      'requested_inputs':{'pdf':args.pdf,'pdf_member':args.pdf_member,'mu_r_factor':args.mu_r_factor,
                          'mu_f_factor':args.mu_f_factor,'call_scale':args.call_scale,'calls_factor':args.calls_factor},
      'flavour_matching':{'profile':args.flavours,'threshold_factor':args.threshold_factor,
                          'alpha_input_scheme':args.alpha_scheme,'alpha_input':args.alpha,
                          'custom_file_sha256':digest(flavour_file) if flavour_file else None,
                          'scope':'Pure-QED charged-fermion evolution; no hadronic vacuum-polarisation model or QCD evolution of quark FFs'},
      'coupling_contract':{'ff':('fixed-alpha photonic, no pairs/running' if photonic else 'one-family MSbar input at mu0=m_mu; no on-shell conversion'),
                           'photos':'fixed alpha, no pair production',
                           'hard_decay':'unitary inclusive decay; host width/branching normalisation'},
      'units':'MCFM histogram integrals in fb; multiply by 1e-3 for pb'}
    if photonic:
        meta['theory']='NLO QCD host; fixed-alpha photonic LL, continuum MSbar D1 restored; no pairs/running'
    if args.hard_samples:
        meta['theory']+='; additive spin-correlated unitary O(alpha) hard decay, leading power'
    if args.dress_hard:
        meta['theory']+='; LL dressing of H1 and one boundary globally; P1/pairs absent'
    if args.photonic_nll:
        meta['theory']='NLO QCD host; fixed-alpha unitary-conditioned photonic NLL; one time-like P1 insertion, boundary or finite spin-correlated H1 globally; continuum NLO overlap; no pairs/running'
        meta['photonic_nll']={'convention':'unitary-conditioned','kernel':'G-G(N=1)*delta',
            'subleading_insertions':'one globally','raw_source_preserved':True}
    if args.coupled_hard:
        meta['theory']+='; coupled one-family NLL hard convolution, all parents, common MSbar alpha(M), continuum overlap'
        meta['photos_comparison_note']='Diagnostic only: PHOTOS has fixed alpha(0) and no pair channels; no OS/MSbar conversion or physical flavour thresholds'
    if args.flavours!='one-family':
        meta['theory']=meta['theory'].replace('one-family','variable-flavour pure-QED')
        meta['coupling_contract']['ff']='NLL MSbar with finite time-like flavour thresholds; '+args.alpha_scheme+' input'
        meta['coupling_contract']['hard_decay']='unitary inclusive decay; same threshold-matched MSbar alpha at original hard mass; host width/branching normalisation'
        meta['photos_comparison_note']='Diagnostic: PHOTOS fixed alpha and no pairs; analytic multi-flavour running/pair content differs'
    if args.flavours=='hadronic-pointlike':
        meta['theory']+='; pointlike light-quark boundary; LO QCD quark/gluon evolution and physical HVP anchor; mixed kernels absent'
        meta['flavour_matching']['scope']='QED kernels with LO QCD evolution, sourced HVP anchor and explicit pointlike light-quark boundary; no mixed kernels or fitted hadronic FF'
        meta['coupling_contract']['ff']='QED NLL kernels plus LO QCD quark/gluon evolution; sourced OS hadronic input and MSbar conversion through alpha*alpha_s'
        meta['hadronic_input']={'file_sha256':digest(hadronic_file),'matching_scale':args.hadronic_scale,'input_standard_deviation_shift':args.hadronic_shift,'alpha_s_override':args.alpha_s,'qcd_evolution_off':args.qcd_off,'boundary':'explicit pointlike approximation, not a fitted hadronic FF'}
        if args.qcd_off:
            meta['theory']=meta['theory'].replace('LO QCD quark/gluon evolution and physical HVP anchor; mixed kernels absent',
                'pure-QED evolution with sourced fermionic HVP anchor; strong FF evolution disabled')
            meta['flavour_matching']['scope']='Pure-QED kernels, sourced HVP anchor and explicit pointlike light-quark boundary; strong FF evolution disabled'
            meta['coupling_contract']['ff']='QED NLL kernels; sourced OS hadronic input and fermionic MSbar conversion through alpha*alpha_s; no strong FF evolution'
    if args.fraction_cells:
        meta['theory']+='; conservative resolved-cell measurement, exact endpoint and number/energy moments'
    if args.production_only:
        meta['theory']='NLO QCD pre-FSR production archive only; no fragmentation or PHOTOS evaluated'
    if args.export_production:
        meta['production_archive']={'file':'production.lff','schema':1,'physics_id':physics_id,
            'run_id':digest(out/'input.ini'),'random_seed':args.seed,'units':'fb',
            'weights':['central'],'full_born_angles':True,'contains_qed':False,
            'scope':'Grouped integration contributions; central production weight only; arbitrary PDF/scale reweighting is not available'}
    path=out/'manifest.json';path.write_text(json.dumps(meta,indent=2)+'\n')
    with (out/'run.log').open('w') as log:
        code=subprocess.run(meta['command'],cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT).returncode
    meta.update(exit_code=code,status='completed' if code==0 else 'failed',
      finished_utc=datetime.now(timezone.utc).isoformat())
    hist={}
    for name in ['total_cross','ff_born','ff_expanded','ff_resummed','ff_matched','ff_shift','ff_beyond_order1','ff_order2',
                 'photos_born','photos_post','photos_shift','photos_minus_ff','photos_radiated',
                 'hard_projected','hard_geometry','hard_shift','decay_order1','decay_matched','photos_minus_decay',
                 'photos_minus_order1','boundary_dressing','projected_dressing','geometry_dressing',
                 'dressing_shift','dressed_matched','photos_minus_dressed','dressed_minus_order1',
                 'coupled_projected','coupled_geometry','coupled_shift','coupled_matched','coupled_order2',
                 'coupled_minus_order1','photos_minus_coupled']:
        files=list(out.glob('*_'+name+'.txt'))
        if len(files)==1:hist[name]=read_histogram(files[0])
    meta['histograms_fb']=hist
    meta['spectra_histograms_fb']={}
    if args.spectra:
        for file in out.glob('*_sp_*.txt'):
            label='sp_'+file.stem.split('_sp_',1)[1]
            meta['spectra_histograms_fb'][label]=read_histogram(file)
        if len(meta['spectra_histograms_fb'])!=136:
            meta['status']='failed';meta['error']='Incomplete physical spectrum histograms'
    logfile=(out/'run.log').read_text()
    resolved_grid=re.findall(r'grid=(\d+), log_span=([^,]+), mu_frag=',logfile)
    meta['resolved_grid']=None
    if resolved_grid:
        meta['resolved_grid']={'bins':int(resolved_grid[-1][0]),'log_span':float(resolved_grid[-1][1])}
        if meta['resolved_grid']!={'bins':args.bins,'log_span':args.log_span}:
            meta['status']='failed';meta['error']='Adapter did not apply requested lattice controls'
    if args.flavours!='one-family':
        resolved=re.findall(r'alpha_msbar_mu0=([^;]+); alpha_msbar_mu_frag=([^;]+); beta0_integral_mu_frag=([^;]+);',logfile)
        if not resolved:meta['status']='failed';meta['error']='Variable-flavour callback absent'
        else:meta['resolved_flavour_coupling']=dict(zip(['alpha_msbar_mu0','alpha_msbar_mu_frag','beta0_integral_mu_frag'],map(float,resolved[-1])))
    chi=re.findall(r'Maximum chisq/it over all contributions is\s+([\d.Ee+\-]+)',logfile)
    meta['maximum_chisq_per_iteration']=float(chi[-1]) if chi else None
    meta['integration_warnings']=[]
    if chi and float(chi[-1])>3:
        meta['integration_warnings'].append('Large final integration chi-squared; repeat/refine before interpreting errors')
    meta['validation_scope']='interface/FF-only pilot; not a finite-hard-matched physical NLL prediction'
    if args.hard_samples:
        meta['validation_scope']='additive NLO decay pilot; finite-geometry cutoffs/statistics require refinement; not hard-dressed NLL'
    if args.dress_hard:
        meta['validation_scope']='one globally dressed finite insertion on LL; exact NLO overlap; no P1/pairs, not full NLL'
    if args.photonic_nll:
        meta['validation_scope']='Unitary-conditioned photonic NLL, fixed alpha and no pairs; numerical refinement and statistical precision remain to establish'
        resolved=re.findall(r'leptonff photonic NLL: convention=unitary-conditioned; source_N1=([^;]+); raw_inclusive_factor=([^;]+); kernel_number_residual=([^;]+); evolved_number_residual=([^;]+); wrap_bound=([^;]+); tail_bound=([^\n]+)',logfile)
        if not resolved:
            meta['status']='failed';meta['error']='Requested photonic NLL callback absent'
        else:
            meta['photonic_nll']['diagnostics']=dict(zip(['source_number_moment','raw_inclusive_factor','kernel_number_residual','evolved_number_residual','wrap_bound','tail_bound'],map(float,resolved[-1])))
    if args.coupled_hard:
        meta['validation_scope']='One-family MSbar NLL hard/FF pilot; no physical flavours/thresholds or precision error certification'
    if args.flavours!='one-family':
        meta['validation_scope']='Variable-flavour pure-QED MSbar NLL with finite thresholds; no hadronic input, coloured-parent QCD evolution or precision certification'
    if args.flavours=='hadronic-pointlike':
        meta['validation_scope']=('Pure-QED evolution with sourced HVP anchor and explicit pointlike coloured boundary; infrared prescription and numerical precision require assessment' if args.qcd_off else
            'Optional LO QCD plus QED evolution with sourced HVP anchor and explicit pointlike coloured boundary; no complete combined-accuracy claim; mixed kernels are not prerequisites to the separate pure-QED benchmark')
        if 'leptonff hadronic pointlike:' not in logfile:meta['status']='failed';meta['error']='Hadronic pointlike callback absent'
    if args.coupled_hard and (not hist.get('coupled_matched') or 'leptonff coupled hard:' not in logfile):
        meta['status']='failed';meta.setdefault('error','Requested coupled hard callback or histograms absent')
    if args.spectra and 'leptonff spectra:' not in logfile:
        meta['status']='failed';meta.setdefault('error','Requested physical spectrum callback absent')
    if args.fraction_cells and 'leptonff fraction cells:' not in logfile:
        meta['status']='failed';meta.setdefault('error','Conservative-cell callback absent')
    if not args.native and not args.production_only and (not all('ff_'+s in hist for s in ['born','expanded','matched','shift'])
        or 'leptonff MCFM:' not in logfile or 'NaN in' in logfile):
        meta['status']='failed';meta.setdefault('error','FF callback/histograms absent; host may have bypassed the analysis hook')
    if args.hard_samples and (not hist.get('decay_matched') or 'leptonff finite decay:' not in logfile):
        meta['status']='failed';meta.setdefault('error','Requested finite-hard callback or histograms absent')
    if args.dress_hard and (not hist.get('dressed_matched') or 'leptonff photonic dressing:' not in logfile):
        meta['status']='failed';meta.setdefault('error','Requested hard dressing callback or histograms absent')
    if not args.production_only and build_info and build_info.get('photos_sha256') and ('leptonff PHOTOS 3.56 comparison:' not in logfile
        or not hist.get('photos_post') or not hist.get('photos_radiated') or
        (args.alpha>0 and hist['photos_radiated'][0][2]<=0)):
        meta['status']='failed';meta.setdefault('error','PHOTOS callback or radiation absent')
    if args.production_only:
        meta['validation_scope']='Production archive export; no QED prediction. Replay and covariance closure required before phenomenology.'
    if args.export_production:
        archive=out/'production.lff'
        complete=False
        if archive.exists() and archive.stat().st_size>=4:
            with archive.open('rb') as file:
                file.seek(-4,os.SEEK_END);complete=file.read()==b'END\n'
        if not complete or 'leptonff production archive:' not in logfile:
            meta['status']='failed';meta['error']='Production archive missing/incomplete'
        else:
            meta['production_archive'].update(sha256=digest(archive),bytes=archive.stat().st_size)
    path.write_text(json.dumps(meta,indent=2)+'\n')
    print(json.dumps({'status':meta['status'],'manifest':str(path),'histograms_fb':hist},indent=2))
    if meta['status']!='completed':raise SystemExit(1)

if __name__=='__main__':main()
