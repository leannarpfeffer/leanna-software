"""Plot arbitrary compatible outputs of wdf_analysis.analyze_wdf.

Requires only retained map/mean arrays, never the full spectral cube or report
roles. Q is an uncalibrated Raman marker response, not PE mass fraction.
"""
from pathlib import Path
import io
import os
import re
import textwrap
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

STYLE = {'font.family': 'serif', 'font.serif': ['Nimbus Roman', 'Times New Roman', 'DejaVu Serif'],
         'font.size': 9, 'axes.titlesize': 9.5, 'axes.labelsize': 9,
         'xtick.labelsize': 8, 'ytick.labelsize': 8, 'text.color': 'black',
         'axes.labelcolor': 'black', 'axes.spines.top': False, 'axes.spines.right': False,
         'mathtext.fontset': 'stix', 'figure.facecolor': 'white', 'savefig.facecolor': 'white'}


def _save(fig, output_dir, stem):
    """Write all formats atomically and return their paths."""
    paths = []
    for extension in ('png', 'pdf', 'svg'):
        buffer = io.BytesIO()
        fig.savefig(buffer, format=extension, dpi=400, bbox_inches='tight', pad_inches=.08)
        path = output_dir / f'{stem}.{extension}'
        temporary = path.with_name(path.name + '.tmp')
        with temporary.open('wb') as stream:
            stream.write(buffer.getvalue()); stream.flush(); os.fsync(stream.fileno())
        os.replace(temporary, path)
        paths.append(str(path))
    plt.close(fig)
    return paths


def plot_raw_maps(data: dict, output_dir: Path) -> list[str]:
    """Save per-file native-Q maps and mean spectra, plus a multi-file overview.

    Raw means include all acquired spectra; cleaned means use cosmic-screened
    eligible spectra. The overview shows median and 10th–90th percentiles of
    valid native-pixel Q values: spatial distributions, not confidence intervals.
    Invalid maps are labeled without assigning a response value.
    """
    output_dir = Path(output_dir).expanduser().resolve()
    if not isinstance(data, dict) or not data:
        raise ValueError('plot_raw_maps requires at least one analyzed-map dictionary.')
    prepared = []
    for index, (key, d) in enumerate(data.items(), 1):
        if not isinstance(d, dict):
            raise ValueError(f'{key}: expected an analyzed-map dictionary.')
        name = str(d.get('name', key))
        try:
            q = np.asarray(d['map_index'], dtype=float).copy()
            x = np.asarray(d['x'], dtype=float)
            raw = np.asarray(d['raw_mean'], dtype=float)
            mean = np.asarray(d['mean'], dtype=float)
            step, width, height = (float(d[k]) for k in ('step', 'width', 'height'))
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError(f'{name}: missing or invalid map geometry/mean-spectrum arrays: {error}') from error
        if q.ndim != 2 or not q.size or x.ndim != 1 or x.size < 2 or raw.shape != x.shape or mean.shape != x.shape:
            raise ValueError(f'{name}: expected a 2-D map and matching 1-D spectral arrays.')
        if not np.isfinite([step, width, height]).all() or step <= 0 or min(width, height) < 0:
            raise ValueError(f'{name}: invalid micrometre geometry.')
        if not np.isfinite(x).all() or not (np.diff(x) > 0).all():
            raise ValueError(f'{name}: Raman-shift axis must be finite and increasing.')
        pixels = d.get('pixels')
        eligible = int(pixels['eligible'].sum()) if pixels is not None and 'eligible' in pixels else None
        if eligible == 0:
            q[:] = np.nan
            mean = np.full_like(mean, np.nan)
        q[~np.isfinite(q) | (q < 0) | (q > 1)] = np.nan
        safe = re.sub(r'[^A-Za-z0-9._-]+', '_', name).strip('._')[:72] or 'map'
        prepared.append((f'{index:03d}_{safe}', name, q, x, raw, mean, step, width, height, eligible))
    output_dir.mkdir(parents=True, exist_ok=True)
    paths, overview = [], []
    with plt.rc_context(STYLE):
        cmap = plt.get_cmap('cividis').copy(); cmap.set_bad('#E5E5E5')
        for stem, name, q, x, raw, mean, step, width, height, eligible in prepared:
            fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.8, 3.65), gridspec_kw={'width_ratios': [1, 1.2]})
            fig.subplots_adjust(left=.09, right=.99, bottom=.22, top=.78, wspace=.40)
            fig.suptitle(textwrap.fill(name, 78), fontsize=10, y=.98)
            im = ax.imshow(q, origin='lower', extent=[-step/2, width+step/2, -step/2, height+step/2],
                           cmap=cmap, vmin=0, vmax=1, interpolation='nearest', aspect='equal')
            ax.set(xlabel='x (µm)', ylabel='y (µm)', title='a  Native Raman marker map')
            ax.xaxis.set_major_locator(MaxNLocator(4)); ax.yaxis.set_major_locator(MaxNLocator(4))
            valid = np.isfinite(q)
            if not valid.any():
                ax.text(.5, .5, 'No eligible spectra' if eligible == 0 else 'No valid marker pixels',
                        transform=ax.transAxes, ha='center', va='center', fontsize=9)
            cb = fig.colorbar(im, ax=ax, orientation='horizontal', pad=.30, fraction=.07, ticks=[0, .5, 1])
            cb.set_label(r'$Q=A^{+}_{1295}/(A^{+}_{1295}+A^{+}_{809})$', labelpad=2)
            if np.isfinite(raw).any():
                bx.plot(x, raw, color='#737A80', lw=.9, label='Raw mean: all spectra')
            if np.isfinite(mean).any():
                bx.plot(x, mean, color='#0072B2', lw=1, label='Cleaned mean: eligible spectra')
            else:
                bx.text(.5, .85, 'No usable cleaned mean', transform=bx.transAxes, ha='center', fontsize=8.5)
            bx.set(xlabel=r'Raman shift (cm$^{-1}$)', ylabel='Mean intensity (counts)', title='b  Mean spectra')
            bx.xaxis.set_major_locator(MaxNLocator(4))
            if bx.lines:
                low, high = bx.get_ylim(); bx.set_ylim(low, high + .25*(high-low))
                bx.legend(frameon=False, fontsize=7.4, loc='upper right')
            fig.text(.50, .02, f'Uncalibrated marker response; native step {step:g} µm. '
                     f'Valid map pixels: {valid.sum():,}/{q.size:,}.\n'
                     'Gray map pixels are excluded. Step size does not establish optical resolution.',
                     ha='center', fontsize=8)
            paths.extend(_save(fig, output_dir, f'{stem}_map_and_spectra'))
            overview.append((name, np.quantile(q[valid], [.1, .5, .9]) if valid.any() else None))
        if len(overview) > 1:
            fig, ax = plt.subplots(figsize=(6.8, max(2.7, .38*len(overview)+1.5)))
            fig.subplots_adjust(left=.36, right=.97, bottom=.31, top=.88)
            for i, (name, values) in enumerate(overview):
                if values is None:
                    ax.text(.5, i, 'No valid marker pixels', ha='center', va='center', color='#737A80', fontsize=8)
                else:
                    ax.plot(values[[0, 2]], [i, i], color='#0072B2', lw=1.6)
                    ax.scatter(values[1], i, s=24, color='#D55E00', zorder=3)
            ax.set_yticks(range(len(overview)), [textwrap.fill(name, 28) for name, _ in overview])
            ax.set_xlim(0, 1); ax.set_ylim(len(overview)-.5, -.5)
            ax.set_xlabel('Valid native-pixel Q (uncalibrated)')
            ax.set_title('Raman marker response across measured fields', pad=10)
            ax.grid(axis='x', color='.90', lw=.5); ax.set_axisbelow(True)
            fig.text(.5, .025, 'Points: median Q. Lines: 10th–90th percentiles over valid pixels.\n'
                     'These describe spatial distributions, not confidence intervals or PE mass fractions.',
                     ha='center', fontsize=8)
            paths.extend(_save(fig, output_dir, 'All_maps_native_Q_distribution'))
    return paths
