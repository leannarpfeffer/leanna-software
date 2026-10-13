"""Shared plotting data, style, and helpers. Edit colors and font settings here.
Native arrays are retained. RAMAN_DATA_DIR and RAMAN_OUTPUT_DIR override the
package's data/ and output/ directories; make_all_figures.py also accepts CLI flags.
"""
from pathlib import Path
import sys, io, os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import patheffects
from matplotlib.patches import Rectangle, FancyBboxPatch, FancyArrowPatch
from matplotlib.ticker import MaxNLocator, ScalarFormatter
from scipy.special import erfc

ROOT=Path(__file__).resolve().parent
DATA=Path(os.environ.get('RAMAN_DATA_DIR', str(ROOT/'data'))).expanduser().resolve()
OUT=Path(os.environ.get('RAMAN_OUTPUT_DIR', str(ROOT/'output'))).expanduser().resolve()
IP=DATA/'E4-small_interface_data'
SP=DATA/'E4-small_structure_data'

def input_error(message):
 raise SystemExit(f'Cannot plot Raman figures: {message}\nSelected data directory: {DATA}\nNo fallback to the bundled data was used.')

def required_path(path):
 path=Path(path)
 if not path.is_file(): input_error(f'missing required input {path}')
 return path

def required_csv(path, columns):
 path=required_path(path)
 try: frame=pd.read_csv(path)
 except (OSError, ValueError, pd.errors.ParserError) as error:
  input_error(f'could not read {path}: {error}')
 missing=set(columns)-set(frame.columns)
 if missing: input_error(f'{path.name} lacks columns: {", ".join(sorted(missing))}')
 if frame.empty: input_error(f'{path.name} contains no rows')
 return frame

try: D=json.loads(required_path(DATA/'map_metadata.json').read_text(encoding='utf-8'))
except (OSError, ValueError) as error: input_error(f'could not read map_metadata.json: {error}')
if not isinstance(D, dict): input_error('map_metadata.json must contain a map-name dictionary')
required_maps=['E0-map2','E1-map','E2-map1','E3-map','E4-map','E5-map','PP-maps','PE-map']
required_maps += [f'E{i}-{suffix}' for i in range(6) for suffix in ('2-1.0','small')]
missing_maps=set(required_maps)-set(D)
if missing_maps: input_error('map_metadata.json lacks maps: '+', '.join(sorted(missing_maps)))
for name in D:
 if not isinstance(D[name], dict): input_error(f'metadata for {name} must be a dictionary')
 missing={'step','width','height','xmin','ymin'}-set(D[name])
 if missing: input_error(f'metadata for {name} lacks: {", ".join(sorted(missing))}')
 path=required_path(DATA/'maps'/f'{name}.npz')
 try:
  with np.load(path,allow_pickle=False) as z: D[name].update({k:z[k] for k in z.files})
 except (OSError, ValueError) as error: input_error(f'could not read {path}: {error}')
 if 'map_index' not in D[name] or np.ndim(D[name]['map_index'])!=2:
  input_error(f'{path.name} must contain a two-dimensional map_index array')
 if name in ('PP-maps','PE-map') and not {'x','corrected'} <= set(D[name]):
  input_error(f'{path.name} lacks the x/corrected reference spectrum arrays')
S=required_csv(DATA/'summary.csv', ['name','bounded_pool']).set_index('name')
MATCH=required_csv(DATA/'matched_footprints.csv', ['sample','r1_coarse_full','r1_coarse_roi','r1_fine']).set_index('sample')
OFF=required_csv(DATA/'sampling_offsets.csv', ['sample','step_um','r1'])
REG=required_csv(DATA/'regional_structure.csv', ['sample','region','PP809_841_shared_baseline'])
for sample in [f'E{i}' for i in range(6)]:
 if sample not in MATCH.index: input_error(f'matched_footprints.csv lacks {sample}')
 for step in (1,2):
  if OFF[(OFF['sample']==sample)&(OFF['step_um']==step)].empty:
   input_error(f'sampling_offsets.csv lacks {sample} at {step} µm')
 if not {'lower_PE_marker','higher_PE_marker'} <= set(REG.loc[REG['sample']==sample,'region']):
  input_error(f'regional_structure.csv lacks both marker quartiles for {sample}')
if not {'PP-maps','PE-map'} <= set(S.index): input_error('summary.csv lacks PP/PE controls')
E4_TRANSECTS=required_csv(IP/'E4-small_interface_transects.csv',
 ['transect','seed_x_um','seed_y_um','normal_x','normal_y','sigma_um','width_10_90_um','width_single_line_um'])
E4_PROFILES=required_csv(IP/'E4-small_interface_profiles.csv',
 ['transect','distance_from_fitted_midpoint_um','normalized_local_contrast'])
E4_SPECTRA=required_csv(SP/'E4-small_structure_mean_spectra.csv',
 ['wavenumber_cm_1','PP_side','interface','PE_side','interface_mixture_prediction'])
for column in ('width_10_90_um','width_single_line_um','sigma_um'):
 values=E4_TRANSECTS[column].to_numpy(dtype=float)
 if not np.isfinite(values).all() or (values<=0).any():
  input_error(f'E4-small_interface_transects.csv contains invalid {column} values')
if set(E4_TRANSECTS.transect) != set(E4_PROFILES.transect):
 input_error('E4 profile and transect identifiers do not match')
# All report inputs are checked before a plotting script writes an output.
try: OUT.mkdir(parents=True,exist_ok=True)
except OSError as error: raise SystemExit(f'Cannot create output directory {OUT}: {error}')
PP='#0072B2'; PE='#D55E00'; INTERFACE='#8E5AA7'; GRAY='#6B747C'
MAP='cividis'
plt.rcParams.update({'font.family':'serif','font.serif':['Nimbus Roman','Times New Roman','DejaVu Serif'],
 'font.size':9,'axes.titlesize':9.5,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,
 'text.color':'black','axes.labelcolor':'black','axes.edgecolor':'black',
 'axes.spines.top':False,'axes.spines.right':False,'mathtext.fontset':'stix',
 'savefig.facecolor':'white','axes.titleweight':'normal','figure.facecolor':'white'})
FILES=[]
def atomic_write(path,blob):
 temp=path.with_name(path.name+'.tmp')
 with open(temp,'wb') as f: f.write(blob); f.flush(); os.fsync(f.fileno())
 os.replace(temp,path)
 assert path.read_bytes()==blob
def save(fig,name):
 buf=io.BytesIO();fig.savefig(buf,format='png',dpi=400,bbox_inches='tight',pad_inches=.035)
 atomic_write(OUT/(name+'.png'),buf.getvalue())
 for extension in ('pdf','svg'):
  buffer=io.BytesIO();fig.savefig(buffer,format=extension,bbox_inches='tight',pad_inches=.035)
  atomic_write(OUT/(name+'.'+extension),buffer.getvalue())
 FILES.append(name);plt.close(fig)
 print('Saved',name,'as PNG, PDF and SVG in',OUT)
def clean(ax):
 ax.grid(axis='x',color='.90',linewidth=.5);ax.set_axisbelow(True)
def mapplot(ax,d,small=False):
 s=d['step'];w=d['width'];h=d['height']
 im=ax.imshow(d['map_index'],origin='lower',extent=[-s/2,w+s/2,-s/2,h+s/2],cmap=MAP,vmin=0,vmax=1,interpolation='nearest',rasterized=True)
 bad=~np.isfinite(d['map_index']);y,x=np.where(bad)
 ax.scatter(x*s,y*s,s=3 if s<.2 else 6,c='#909090',marker='x',linewidths=.35)
 ax.set_xlim(-s/2,w+s/2);ax.set_ylim(-s/2,h+s/2)
 ax.set_aspect('equal');ax.xaxis.set_major_locator(MaxNLocator(3));ax.yaxis.set_major_locator(MaxNLocator(3))
 ax.set_xlabel('x (µm)',labelpad=1.8);ax.set_ylabel('y (µm)',labelpad=1.8)
 ax.tick_params(pad=1.5,length=2.5,width=.6)
 for z in ax.spines.values():z.set_visible(True);z.set_linewidth(.6)
 return im


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

def local_corrected(s, x, spec):
    _, left, right, xl, xr = band_weights(x, spec)
    yl, yr = (s[:, left].mean(1), s[:, right].mean(1))
    base = yl[:, None] + (yr - yl)[:, None] * (x - xl) / (xr - xl)
    return s - base
