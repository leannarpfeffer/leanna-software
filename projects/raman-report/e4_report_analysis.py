"""Recompute the report's exact E4 boundary and regional spectral recipe.

This is intentionally NOT a generic boundary detector. Its row range, seven
transects, fit bounds and region distances were selected for one E4-small field.
``run_e4_report`` hashes the supplied WDF before reading it or creating outputs
and refuses every source except the byte-identical report acquisition. Renaming
that acquisition is fine; applying its coordinates to a different map is not.
For another acquisition, select and validate that map's own boundary/plateaus
and implement a separate explicitly documented recipe. Do not remove the hash
guard to force a new map through this one.

Every boundary seed is fitted again from the raw WDF. No retained row-fit CSV,
regional spectrum, or plot cache is used as an analysis input. This module
writes numerical CSV/JSON outputs only; it does not generate figures, fit PP
peak shapes, run significance tests, or infer absolute crystallinity.

The original E4 numerical definitions are retained, including conservative
spike detection, exact-bound band areas, native-pixel regional means, robust
error-function fits and the empirical two-side mixture with held-out markers.
The two side means contain both polymers; their fitted coefficients are not
concentrations. Apparent Raman width is not intrinsic interphase thickness.
Nonlinear-fit roundoff can vary slightly with NumPy/SciPy and BLAS settings;
retain full-precision numerical outputs rather than rounding to a stored result.
"""
from pathlib import Path
import csv
import hashlib
import json

import numpy as np
import pandas as pd
from scipy.interpolate import UnivariateSpline, RegularGridInterpolator
from scipy.ndimage import median_filter
from scipy.optimize import least_squares, lsq_linear
from scipy.signal import find_peaks
from scipy.special import erfc, expit
from scipy.stats import norm
from renishawWiRE import WDFReader

__all__ = ["run_e4_report"]
EXPECTED_E4_SHA256 = "bc1a23788ccd2fe3b47d05a67e8e8213ca9c551ea8de9d3f269a56d166f78b83"
RECIPE_ID = "E4-small_report_boundary_v1"

BANDS = {'PE1295': {'integration': [1283, 1308], 'left': [1265, 1278], 'right': [1310, 1318]}, 'PP809': {'integration': [800, 820], 'left': [775, 789], 'right': [822, 832]}, 'PE1130_check': {'integration': [1120, 1138], 'left': [1103, 1115], 'right': [1138, 1143]}, 'PP_broad_check': {'integration': [793, 858], 'left': [775, 789], 'right': [862, 874]}}

def band_weights(x, spec):
    """Exact-bound trapezoid area, minus a line between sideband means."""
    lo, hi = spec['integration']
    grid = np.r_[lo, x[(x > lo) & (x < hi)], hi]
    t = np.zeros(len(grid))
    d = np.diff(grid)
    t[:-1] += d / 2
    t[1:] += d / 2
    w = np.zeros(len(x))
    for xv, tv in zip(grid, t):
        k = np.clip(np.searchsorted(x, xv), 1, len(x) - 1)
        f = (xv - x[k - 1]) / (x[k] - x[k - 1])
        w[k - 1] += tv * (1 - f)
        w[k] += tv * f
    left = (x >= spec['left'][0]) & (x <= spec['left'][1])
    right = (x >= spec['right'][0]) & (x <= spec['right'][1])
    xl, xr = (x[left].mean(), x[right].mean())
    cr = ((hi ** 2 - lo ** 2) / 2 - xl * (hi - lo)) / (xr - xl)
    cl = hi - lo - cr
    w[left] -= cl / left.sum()
    w[right] -= cr / right.sum()
    assert abs(w.sum()) < 1e-08 and abs(w @ x) < 1e-06
    return (w, left, right, xl, xr)


def clean_spikes(s, x):
    """Conservatively remove only very large, narrow positive artifacts."""
    quiet = (x > 1510) & (x < 1740)
    diff = np.diff(s[:, quiet], axis=1)
    noise = 1.4826 * np.median(abs(diff - np.median(diff, axis=1)[:, None]), axis=1) / np.sqrt(2)
    residual = s - median_filter(s, size=(1, 11), mode='nearest')
    out = s.copy()
    corrected = np.zeros(s.shape, dtype=bool)
    events = []
    for i in range(len(s)):
        peaks, _ = find_peaks(residual[i], height=max(100, 20 * noise[i]))
        for k in peaks:
            a = b = int(k)
            while a > 1 and residual[i, a - 1] > 5 * noise[i]:
                a -= 1
            while b < len(x) - 2 and residual[i, b + 1] > 5 * noise[i]:
                b += 1
            if b - a + 1 > 5:
                continue
            out[i, a:b + 1] = np.interp(x[a:b + 1], x[[a - 1, b + 1]], s[i, [a - 1, b + 1]])
            corrected[i, a:b + 1] = True
            events.append({'spectrum_index': i, 'center_cm_1': float(x[k]), 'first_channel': a, 'last_channel': b, 'original_peak_counts': float(s[i, k])})
    return (out, corrected, events, noise)

WFACTOR=2*norm.ppf(.9)

def model(p,t,kind='erf'):
    low,high,center,scale=p
    transition=.5*erfc((t-center)/(np.sqrt(2)*scale)) if kind=='erf' else expit(-(t-center)/scale)
    return low+(high-low)*transition

def fit(t,v,kind='erf'):
    good=np.isfinite(v);t=t[good];v=v[good]
    low=np.median(v[-10:]);high=np.median(v[:10]);center=t[np.argmin(abs(v-(low+high)/2))]
    p0=[np.clip(low,.001,.349),np.clip(high,.351,.999),center,.5]
    result=least_squares(lambda p:model(p,t,kind)-v,p0,
        bounds=([0,.35,t.min(),.05],[.35,1,t.max(),2]),loss='soft_l1',f_scale=.035)
    residual=model(result.x,t,kind)-v
    return result.x,float(np.sqrt(np.mean(residual**2))),float(1-np.sum(residual**2)/np.sum((v-v.mean())**2))

def _run_interface(source, out):
    r=WDFReader(str(source));o=np.argsort(r.xdata);wn=r.xdata[o].astype(float)
    raw=r.spectra.reshape(-1,len(wn))[:,o].astype(float)
    s,corrected,events,noise=clean_spikes(raw,wn)
    weights={k:band_weights(wn,v)[0] for k,v in BANDS.items()}
    areas={k:s@v for k,v in weights.items()}
    p=np.maximum(areas['PE1295'],0);q=np.maximum(areas['PP809'],0);a1130=np.maximum(areas['PE1130_check'],0)
    support=(weights['PE1295']!=0)|(weights['PP809']!=0)
    masked=corrected[:,support].any(1)
    ux,uy=np.unique(r.xpos),np.unique(r.ypos);x=ux-ux[0];y=uy-uy[0]
    col,row=np.searchsorted(ux,r.xpos),np.searchsorted(uy,r.ypos)
    assert len(x)==91 and len(y)==81 and np.allclose(np.diff(x),.1) and np.allclose(np.diff(y),.1)
    # Rounded analytical coordinate grid avoids platform-dependent floating-point offsets.
    x=np.arange(len(ux))*.1;y=np.arange(len(uy))*.1
    vals={'main':p/(p+q),'1130':a1130/(a1130+q),'broad_PP':p/(p+np.maximum(areas['PP_broad_check'],0))}
    grids={}
    for name,v in vals.items():
        v=v.copy();v[masked]=np.nan;g=np.full((len(y),len(x)),np.nan);g[row,col]=v;grids[name]=g
    g=grids['main']
    rowfit=[]
    for j in range(10,71):
        par,rmse,r2=fit(x,g[j]);rowfit.append([y[j],*par,rmse,r2])
    rowfit=np.array(rowfit)
    curve=UnivariateSpline(rowfit[:,0],rowfit[:,3],s=len(rowfit)*.05**2)
    interpolators={k:RegularGridInterpolator((y,x),v,bounds_error=False,fill_value=np.nan) for k,v in grids.items()}
    distances=np.arange(-2.5,2.51,.1);offsets=np.arange(-.2,.201,.1)
    centers=np.arange(1.75,6.26,.75);results=[];profiles=[];fits=[];paths=[]
    for i,y0 in enumerate(centers):
        x0=float(curve(y0));slope=float(curve.derivative()(y0))
        nv=np.array([1,-slope])/np.sqrt(1+slope*slope);tv=np.array([slope,1])/np.sqrt(1+slope*slope)
        scans={}
        for name,interp in interpolators.items():
            strips=[]
            for off in offsets:
                pts=np.array([x0,y0])+distances[:,None]*nv+off*tv
                strips.append(interp(pts[:,::-1]))
            scans[name]=np.array(strips)
        assert np.isfinite(scans['main']).all()
        data=scans['main'].mean(0);par,rmse,r2=fit(distances,data)
        single,_,_=fit(distances,scans['main'][2]);alt1130,_,_=fit(distances,scans['1130'].mean(0))
        altbroad,_,_=fit(distances,scans['broad_PP'].mean(0));logistic,_,_=fit(distances,data,'logistic')
        low,high,shift,sigma=par;w=WFACTOR*sigma
        centerxy=np.array([x0,y0])+shift*nv;p90=centerxy-.5*w*nv;p10=centerxy+.5*w*nv
        result={'transect':i+1,'seed_y_um':y0,'seed_x_um':x0,'boundary_slope_dx_dy':slope,
            'normal_x':float(nv[0]),'normal_y':float(nv[1]),'center_x_um':float(centerxy[0]),'center_y_um':float(centerxy[1]),
            'lower_plateau_index':low,'upper_plateau_index':high,'normal_center_shift_um':shift,'sigma_um':sigma,
            'width_10_90_um':w,'max_index_gradient_per_um':(high-low)/(np.sqrt(2*np.pi)*sigma),
            'width_single_line_um':single[3]*WFACTOR,'width_1130_check_um':alt1130[3]*WFACTOR,
            'width_broad_PP_check_um':altbroad[3]*WFACTOR,'width_logistic_um':logistic[3]*2*np.log(9),
            'rmse_index':rmse,'r_squared':r2,'x_90_um':float(p90[0]),'y_90_um':float(p90[1]),
            'x_10_um':float(p10[0]),'y_10_um':float(p10[1])}
        results.append(result);profiles.append((data,scans));fits.append(par)
        paths.append(np.array([x0,y0])+distances[:,None]*nv)
    widths=np.array([r['width_10_90_um'] for r in results]);med=float(np.median(widths))
    summary={'source':'E4-small.wdf','sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
      'step_um':.1,'laser_nm':float(r.laser_length),'objective_magnification':100,
      'quantity':'Apparent lateral Raman-marker transition width; not intrinsic interphase thickness',
      'width_definition':'Distance over which fitted local Raman-index contrast falls from 90% to 10% of the difference between fitted plateaus',
      'transects':7,'transect_seed_y_um':centers.tolist(),'same_field_not_independent_replicates':True,
      'normal_sample_spacing_um':.1,'tangential_offsets_um':offsets.tolist(),'profile_normal_range_um':[-2.5,2.5],
      'boundary_seed_method':'Per-row erf fits at y=1.0..7.0 um; cubic smoothing spline of fitted midpoints with s=61*(0.05 um)^2',
      'fit_method':'Four-parameter complementary error function; soft_l1 robust least squares, f_scale=0.035 index units',
      'main_median_width_um':med,'main_IQR_um':np.quantile(widths,[.25,.75]).tolist(),
      'main_min_max_um':[float(widths.min()),float(widths.max())],
      'sensitivity_median_widths_um':{k:float(np.median([z[k] for z in results])) for k in
         ['width_single_line_um','width_1130_check_um','width_broad_PP_check_um','width_logistic_um']},
      'median_max_index_gradient_per_um':float(np.median([z['max_index_gradient_per_um'] for z in results])),
      'masked_spectrum_indices':np.where(masked)[0].tolist(),'PSF_correction':'Not performed; no measured spatial response supplied',
      'interpretation_limits':['Uncalibrated Raman marker ratio, not PE wt%.','Optical blur and 3D sample geometry can broaden the transition.',
          'Sampling interval is not optical resolution.','Observed width does not establish interdiffusion, molecular miscibility, adhesion or chain entanglement.'],
      'instrument_notes':{'WXIS_Motor_Positions_Pinhole':'Out','WXDM_Pinhole_In':0,
        'WXIS_Slits_Opening':'50 um','WXDM_InstrumentState_slitOpening':10.0,
        'interpretation':'Do not infer nonconfocal acquisition from Pinhole Out alone. Renishaw EasyConfocal can operate without a fixed pinhole; effective spatial response must be measured.'}}
    (out/'E4-small_interface_summary.json').write_text(json.dumps(summary,indent=2))
    with (out/'E4-small_interface_transects.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(results[0]));wr.writeheader();wr.writerows(results)
    with (out/'E4-small_interface_profiles.csv').open('w',newline='') as f:
        wr=csv.writer(f);wr.writerow(['transect','distance_from_seed_um','distance_from_fitted_midpoint_um','PE_index',
            'normalized_local_contrast','fitted_index','fitted_normalized_contrast','single_line_index','index_1130_check','index_broad_PP_check'])
        for j,((data,scans),par) in enumerate(zip(profiles,fits)):
            low,high,c,sigma=par;pred=model(par,distances)
            for k,t in enumerate(distances):wr.writerow([j+1,t,t-c,data[k],(data[k]-low)/(high-low),pred[k],
                (pred[k]-low)/(high-low),scans['main'][2,k],scans['1130'][:,k].mean(),scans['broad_PP'][:,k].mean()])
    np.savetxt(out/'E4-small_boundary_seed_fits.csv',rowfit,delimiter=',',header='y_um,low,high,x50_um,sigma_um,rmse,R2',comments='')

    return summary, rowfit, s, corrected, wn, r

PP_BASE={'integration':[790,860],'left':[765,785],'right':[865,885]}
SPECS={
    'PP809':dict(PP_BASE,integration=[800,819]),
    'PP841':dict(PP_BASE,integration=[833,849]),
    'PP_envelope':PP_BASE,
    'PE1295':BANDS['PE1295'],
    'PE1416':{'integration':[1408,1424],'left':[1398,1406],'right':[1426,1431]}
}
REGIONS={'PE_side':(-2.5,-1.5),'interface':(-.5,.5),'PP_side':(1.5,2.5)}

def _run_structure(source, out, rowfits, spectra, corrected, wn, r):
    xy=np.c_[r.xpos-min(r.xpos),r.ypos-min(r.ypos)]
    # Boundary seeds below were recomputed from this same raw WDF in memory.
    curve=UnivariateSpline(rowfits[:,0],rowfits[:,3],s=len(rowfits)*.05**2)
    cy=np.linspace(.5,7.5,1401);points=np.c_[curve(cy),cy]
    nearest=np.argmin(np.sum((xy[:,None,:]-points[None,:,:])**2,axis=2),axis=1)
    foot=cy[nearest]
    distance=np.sqrt(np.sum((xy-points[nearest])**2,axis=1))*np.sign(xy[:,0]-curve(xy[:,1]))
    masked=corrected[:,(wn>=760)&(wn<=1500)].any(1)
    valid=(foot>=1.375)&(foot<6.625)&~masked
    masks=[valid&(distance>=lo)&(distance<hi) for lo,hi in REGIONS.values()]
    region_spectra=np.array([spectra[m].mean(0) for m in masks])
    weights={name:band_weights(wn,spec)[0] for name,spec in SPECS.items()}
    areas={name:region_spectra@v for name,v in weights.items()}
    pp_ratio=areas['PP809']/areas['PP841']
    pp_fraction=areas['PP809']/areas['PP_envelope']
    def record(sp,**meta):
        a={name:float(sp@v) for name,v in weights.items()}
        return dict(meta,**a,PP809_over_PP841=a['PP809']/a['PP841'],
                    PP809_over_envelope=a['PP809']/a['PP_envelope'])
    region_records=[record(sp,region=name,n_pixels=int(mask.sum()),distance_low_um=REGIONS[name][0],
        distance_high_um=REGIONS[name][1]) for name,sp,mask in zip(REGIONS,region_spectra,masks)]
    pd.DataFrame(region_records).to_csv(out/'E4-small_structure_regions.csv',index=False)
    # Test an empirical two-side spectral mixture. Target bands are held out.
    # The two side means are mixed local spectra, not pure polymer endmembers.
    M=np.c_[region_spectra[0],region_spectra[2],np.ones(len(wn)),(wn-1100)/500]
    fitmask=(wn>=900)&(wn<=1250)
    fit=lsq_linear(M[fitmask],region_spectra[1,fitmask],bounds=([0,0,-np.inf,-np.inf],[np.inf]*4))
    predicted=M@fit.x
    mixture_metrics={}
    for name,v in weights.items():
        observed=float(region_spectra[1]@v);prediction=float(predicted@v)
        mixture_metrics[name]={'observed_area':observed,'predicted_area':prediction,
            'difference_percent_of_prediction':100*(observed/prediction-1)}
    observed_pe_ratio=areas['PE1416'][1]/areas['PE1295'][1]
    predicted_pe_ratio=(predicted@weights['PE1416'])/(predicted@weights['PE1295'])
    mixture_pe_difference=100*(observed_pe_ratio/predicted_pe_ratio-1)
    # PP-template subtraction is a sensitivity check, not a validated separation.
    subtraction=[];subtracted_spectra={}
    for name in ['PP809','PP841','PP_envelope']:
        v=weights[name];scales=region_spectra@v/(region_spectra[2]@v)
        sub=region_spectra[:2]-scales[:2,None]*region_spectra[2]
        ratios=sub@weights['PE1416']/(sub@weights['PE1295'])
        subtraction.append({'PP_scaling_band':name,'PE_side_index':ratios[0],
            'interface_index':ratios[1],'relative_change_percent':100*(ratios[1]/ratios[0]-1)})
        subtracted_spectra[name]=sub
    spectra_export={'wavenumber_cm_1':wn,**{name:sp for name,sp in zip(REGIONS,region_spectra)},
        'interface_mixture_prediction':predicted,'interface_mixture_residual':region_spectra[1]-predicted}
    for name,sub in subtracted_spectra.items():
        spectra_export[f'PE_side_minus_PP_scaled_by_{name}']=sub[0]
        spectra_export[f'interface_minus_PP_scaled_by_{name}']=sub[1]
    pd.DataFrame(spectra_export).to_csv(out/'E4-small_structure_mean_spectra.csv',index=False)
    pp_change=float(100*(pp_ratio[1]/pp_ratio[2]-1))
    summary={
        'recipe_id':RECIPE_ID,
        'source':source.name,
        'sha256':EXPECTED_E4_SHA256,
        'sampling_um':.1,
        'spatial_selection':'Nearest boundary point at y=1.375..6.625 um',
        'boundary_definition':'Relative PE-index midpoint spline recomputed from raw WDF row fits in this run',
        'boundary_sampling_for_nearest_point_um':.005,
        'regions_um':REGIONS,
        'band_definitions_cm_1':SPECS,
        'region_counts':{name:int(m.sum()) for name,m in zip(REGIONS,masks)},
        'PP809_over_841':dict(zip(REGIONS,pp_ratio.tolist())),
        'PP809_over_envelope':dict(zip(REGIONS,pp_fraction.tolist())),
        'PP_interface_ratio_change_percent':pp_change,
        'empirical_mixture_fit_window_cm_1':[900,1250],
        'empirical_mixture_coefficients_not_concentrations':fit.x.tolist(),
        'mixture_marker_checks':mixture_metrics,
        'PE_ratio_vs_mixture_prediction_percent':float(mixture_pe_difference),
        'PP_subtraction_sensitivity':subtraction,
        'spike_masked_spectrum_indices':np.where(masked)[0].tolist(),
        'absolute_crystallinity':'Not determined',
        'interpretation':'Small PP band-ratio increase; no clear field-averaged PE crystallinity-specific change beyond spectral mixing. Orientation, overlapping bands, local heterogeneity and optical averaging remain unresolved.',
    }
    (out/'E4-small_structure_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    return summary


def run_e4_report(source: Path, data_dir: Path) -> dict:
    """Write exact-report E4 numerical outputs, deriving all seeds from raw data.

    Parameters
    ----------
    source : Path
        Byte-identical report E4-small.wdf (the filename may differ).
    data_dir : Path
        Parent output directory. Creates ``E4-small_interface_data`` and
        ``E4-small_structure_data`` with the report-compatible filenames.

    Returns
    -------
    dict
        Recipe identity, source SHA-256, output filenames, and the two numerical
        summaries. Existing matching output files are replaced only after the
        source hash has passed. No retained analysis input is read.

    Raises
    ------
    ValueError
        If the raw file differs from the acquisition for which this spatial
        recipe was validated. A new boundary requires its own selection.
    """
    source=Path(source)
    data_dir=Path(data_dir)
    digest=hashlib.sha256()
    with source.open('rb') as handle:
        for chunk in iter(lambda:handle.read(1024*1024),b''):
            digest.update(chunk)
    actual=digest.hexdigest()
    if actual != EXPECTED_E4_SHA256:
        raise ValueError(
            'This E4 boundary recipe is valid only for the original report '
            f'E4-small.wdf (SHA-256 {EXPECTED_E4_SHA256}). '
            f'The supplied file has SHA-256 {actual}. '
            "Select and validate this new map's own boundary and plateau "
            'regions, then create a separate recipe; do not reuse these '
            'hardcoded spatial selections or remove this guard.'
        )
    interface_dir=data_dir/'E4-small_interface_data'
    structure_dir=data_dir/'E4-small_structure_data'
    interface_dir.mkdir(parents=True,exist_ok=True)
    structure_dir.mkdir(parents=True,exist_ok=True)
    boundary,rowfits,spectra,corrected,wn,reader=_run_interface(source,interface_dir)
    regional=_run_structure(source,structure_dir,rowfits,spectra,corrected,wn,reader)
    files=[
        interface_dir/'E4-small_interface_transects.csv',
        interface_dir/'E4-small_interface_profiles.csv',
        interface_dir/'E4-small_boundary_seed_fits.csv',
        interface_dir/'E4-small_interface_summary.json',
        structure_dir/'E4-small_structure_mean_spectra.csv',
        structure_dir/'E4-small_structure_regions.csv',
        structure_dir/'E4-small_structure_summary.json',
    ]
    return {'recipe_id':RECIPE_ID,'source_sha256':actual,
            'files':[str(p) for p in files],
            'interface':boundary,'structure':regional}
