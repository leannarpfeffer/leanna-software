"""Raw-WDF Raman analysis used by the report; indices are not composition.

Only complete, regular, square-spaced two-dimensional maps with count-valued
intensities, a Raman-shift axis in cm^-1 and spatial coordinates in micrometres
are supported. Native
pixels are retained without spatial smoothing or binning. Mean-spectrum ALS
correction is for display; quantitative areas use the local baselines below.
"""
from pathlib import Path
import datetime
import hashlib
import re
from zoneinfo import ZoneInfo, ZoneInfoNotFoundError
import numpy as np
import pandas as pd
from scipy.ndimage import median_filter
from scipy.linalg import solveh_banded
from scipy.stats import spearmanr
from renishawWiRE import WDFReader
from renishawWiRE.types import DataType, UnitType

BANDS = {
    '1295': ((1285, 1302), [(1268, 1280), (1310, 1318)]),
    '1330': ((1320, 1340), [(1310, 1318), (1344, 1352)]),
    '1128': ((1123, 1133), [(1116, 1121), (1135, 1139)]),
    '974': ((963, 984), [(947, 956), (987, 992)]),
    '1062': ((1059, 1069), [(1051, 1057), (1072, 1078)]),
    '809': ((798, 820), [(779, 791), (824, 831)]),
    '841': ((833, 853), [(824, 831), (860, 871)]),
    '1416': ((1409, 1423), [(1396, 1404), (1425, 1431)]),
    'twist': ((1280, 1314), [(1263, 1272), (1316, 1320)]),
}
BOUNDED_BANDS = {
    'PE1295': {'integration': [1283, 1308], 'left': [1265, 1278], 'right': [1310, 1318]},
    'PP809': {'integration': [800, 820], 'left': [775, 789], 'right': [822, 832]},
    'PE1130_check': {'integration': [1120, 1138], 'left': [1103, 1115], 'right': [1138, 1143]},
    'PP_broad_check': {'integration': [793, 858], 'left': [775, 789], 'right': [862, 874]},
}
RATIOS = {'r1': ('1295', '1330'), 'r2': ('1128', '974'),
          'r3': ('1062', '974'), 'pp_order': ('809', '841'),
          'pe_order': ('1416', 'twist')}
SUMMARY_KEYS = ['name', 'file', 'sample', 'nx', 'ny', 'n', 'step', 'width',
                'height', 'xmin', 'ymin', 'xmax', 'ymax', 'acquired_local',
                'duration_s', 'laser', 'complete', 'uuid', 'sha256']
REGIONAL_COLUMNS = ['sample', 'region', 'n_pixels', 'selection_score_low_quartile',
    'selection_score_high_quartile', 'r1_pooled', 'r2_pooled', 'PP809_841_local',
    'PP809_841_shared_baseline', 'PE1416_1295', 'PE1416_area', 'PE1295_area',
    'PP809_area', 'PP841_area']


def _axis(x):
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or x.size < 3 or not np.isfinite(x).all() or not (np.diff(x) > 0).all():
        raise ValueError('The spectral axis must be finite, unique and strictly increasing.')
    return x


def _window(x, bounds, minimum=2):
    lo, hi = bounds
    if lo >= hi or x[0] > lo or x[-1] < hi:
        raise ValueError(f'Spectral coverage does not contain required window {lo}–{hi} cm^-1.')
    mask = (x >= lo) & (x <= hi)
    if mask.sum() < minimum:
        raise ValueError(f'Too few spectral channels in required window {lo}–{hi} cm^-1.')
    return mask


def band_weights(x, band, anchors):
    """Native-channel trapezoid area with a least-squares sideband baseline."""
    x = _axis(x)
    b = _window(x, band)
    a = np.zeros(len(x), dtype=bool)
    for anchor in anchors:
        a |= _window(x, anchor)
    xb = x[b]
    tw = np.empty(len(xb))
    tw[0] = (xb[1] - xb[0]) / 2
    tw[-1] = (xb[-1] - xb[-2]) / 2
    tw[1:-1] = (xb[2:] - xb[:-2]) / 2
    w = np.zeros(len(x))
    w[b] = tw
    design = np.column_stack((np.ones(len(x)), (x - 1300) / 100))
    w[a] -= (tw @ design[b]) @ np.linalg.pinv(design[a])
    return w


def exact_band_weights(x, spec):
    """Exact-bound trapezoid minus a line between two sideband means.

    Returns (weights, left_mask, right_mask, left_mean_x, right_mean_x),
    matching the retained local E4 method.
    """
    x = _axis(x)
    lo, hi = spec['integration']
    _window(x, (lo, hi))
    left = _window(x, spec['left'])
    right = _window(x, spec['right'])
    grid = np.r_[lo, x[(x > lo) & (x < hi)], hi]
    t = np.zeros(len(grid))
    delta = np.diff(grid)
    t[:-1] += delta / 2
    t[1:] += delta / 2
    w = np.zeros(len(x))
    for xv, tv in zip(grid, t):
        k = np.clip(np.searchsorted(x, xv), 1, len(x) - 1)
        f = (xv - x[k - 1]) / (x[k] - x[k - 1])
        w[k - 1] += tv * (1 - f)
        w[k] += tv * f
    xl, xr = x[left].mean(), x[right].mean()
    if not xr > xl:
        raise ValueError('Baseline sidebands must have increasing, distinct centers.')
    cr = ((hi ** 2 - lo ** 2) / 2 - xl * (hi - lo)) / (xr - xl)
    cl = hi - lo - cr
    w[left] -= cl / left.sum()
    w[right] -= cr / right.sum()
    return w, left, right, xl, xr


def despike(s, x):
    """Replace large positive one/two-channel events; return clean, noise, flags."""
    x = _axis(x)
    s = np.asarray(s, dtype=float)
    if s.ndim != 2 or s.shape[1] != len(x) or not np.isfinite(s).all():
        raise ValueError('Spectra must be a finite two-dimensional array matching the spectral axis.')
    quiet = _window(x, (1500, 1700), minimum=8)
    def estimate(a):
        second = np.diff(a[:, quiet], n=2, axis=1)
        return np.maximum(1., 1.4826 * np.median(np.abs(second), axis=1) / np.sqrt(6))
    smooth = median_filter(s, size=(1, 5), mode='nearest')
    noise = estimate(s)
    threshold = np.maximum(np.maximum(12 * noise[:, None], 5 * np.sqrt(np.maximum(s, 0))), 40.)
    candidate = (s - smooth) > threshold
    fixed = np.zeros_like(candidate)
    clean = s.copy()
    for i in np.flatnonzero(candidate.any(axis=1)):
        ids = np.flatnonzero(candidate[i])
        for group in np.split(ids, np.where(np.diff(ids) > 1)[0] + 1):
            if 0 < len(group) <= 2:
                fixed[i, group] = True
                clean[i, group] = smooth[i, group]
    return clean, estimate(clean), fixed


def als(y, lam=1e5, p=.001, niter=12):
    """Asymmetric least-squares display baseline, preserving the original method."""
    y = np.asarray(y, dtype=float)
    if y.ndim != 1 or len(y) < 4 or not np.isfinite(y).all():
        raise ValueError('ALS requires at least four finite spectral values.')
    if lam <= 0 or not 0 < p < 1 or niter < 1:
        raise ValueError('ALS requires positive lambda, 0 < p < 1, and at least one iteration.')
    n = len(y)
    diag = np.full(n, 6.)
    diag[:2], diag[-2:] = [1, 5], [5, 1]
    off = np.full(n - 1, -4.)
    off[[0, -1]] = -2
    weights = np.ones(n)
    for _ in range(niter):
        ab = np.zeros((3, n))
        ab[2] = weights + lam * diag
        ab[1, 1:] = lam * off
        ab[0, 2:] = lam
        baseline = solveh_banded(ab, weights * y, lower=False, check_finite=False)
        weights = np.where(y > baseline, p, 1 - p)
    return baseline


def _divide(a, b, positive=False):
    """Safe scalar division; signed nonzero denominators remain valid for pooled indices."""
    if not np.isfinite(a) or not np.isfinite(b) or (b <= 0 if positive else b == 0):
        return float('nan')
    return float(a / b)


def _add_ratios(df):
    for key, (a, b) in RATIOS.items():
        numerator, denominator = df['A' + a].to_numpy(), df['A' + b].to_numpy()
        df[key] = np.divide(numerator, denominator, out=np.full(len(df), np.nan), where=denominator > 0)
        df['valid_' + key] = (denominator / np.sqrt(df['V' + b]) >= 5) & np.isfinite(df[key]) & df.eligible
        df['detected_' + key] = (numerator / np.sqrt(df['V' + a]) >= 3) & df['valid_' + key]
    return df


def pool(df, key):
    """Ratio of summed areas over eligible pixels; never the mean of pixel ratios."""
    if 'eligible' in df:
        df = df[df.eligible]
    if not len(df):
        return float('nan')
    a, b = RATIOS[key]
    return _divide(df['A' + a].sum(), df['A' + b].sum())


def _rho(a, b):
    a, b = np.asarray(a), np.asarray(b)
    good = np.isfinite(a) & np.isfinite(b)
    if good.sum() < 5 or np.ptp(a[good]) == 0 or np.ptp(b[good]) == 0:
        return float('nan')
    return float(spearmanr(a[good], b[good]).statistic)


def _geometry(reader):
    if not getattr(reader, 'is_completed', False):
        raise ValueError('Incomplete WDF acquisition: partial maps are unsupported.')
    count = int(reader.count)
    shape = getattr(reader, 'map_shape', None)
    if shape is None or len(shape) != 2 or min(shape) < 2 or count < 4:
        raise ValueError('A complete 2D map is required; single spectra and line scans are unsupported.')
    nx, ny = map(int, shape)
    if count != nx * ny:
        raise ValueError('Map spectrum count does not match its declared dimensions.')
    if getattr(reader, 'xpos_unit', None) != UnitType.Micron or getattr(reader, 'ypos_unit', None) != UnitType.Micron:
        raise ValueError('Map X and Y coordinates must be in micrometres.')
    xp, yp = np.asarray(reader.xpos, dtype=float), np.asarray(reader.ypos, dtype=float)
    if xp.shape != (count,) or yp.shape != (count,) or not (np.isfinite(xp).all() and np.isfinite(yp).all()):
        raise ValueError('Missing, non-finite, or mismatched map coordinates.')
    info = getattr(reader, 'map_info', {})
    sx, sy = abs(float(info.get('x_pad', 0))), abs(float(info.get('y_pad', 0)))
    if not sx > 0 or not sy > 0 or not np.isfinite([sx, sy]).all() or not np.isclose(sx, sy, rtol=1e-5, atol=1e-7):
        raise ValueError('Map must have positive, equal X/Y sampling intervals.')
    ix = np.rint((xp - xp.min()) / sx).astype(int)
    iy = np.rint((yp - yp.min()) / sy).astype(int)
    # Tolerance covers stored float32 step precision, not visibly irregular grids.
    tolerance = max(1e-5, sx * 1e-4)
    if not np.allclose(xp, xp.min() + ix * sx, rtol=0, atol=tolerance) or not np.allclose(yp, yp.min() + iy * sy, rtol=0, atol=tolerance):
        raise ValueError('Coordinates do not lie on the declared regular map grid.')
    if ix.min() != 0 or iy.min() != 0 or ix.max() != nx - 1 or iy.max() != ny - 1 or len(set(zip(ix, iy))) != count:
        raise ValueError('Map has duplicate, missing, or out-of-grid coordinate positions.')
    return nx, ny, sx, xp, yp, ix, iy


def _metadata(reader, path):
    meta = dict(acquired_local=None, acquired_utc=None, duration_s=None, uuid=None,
                accumulations=getattr(reader, 'accumulation_count', None),
                laser=float(reader.laser_length) if reader.laser_length else None,
                scan_type=str(reader.scan_type), measurement_type=str(reader.measurement_type),
                spectral_axis_type=str(reader.xlist_type), spectral_axis_unit=str(reader.xlist_unit),
                intensity_unit=str(reader.spectral_unit), application_name=reader.application_name,
                application_version=list(reader.application_version), acquisition_title=reader.title)
    try:
        from wdf import Wdf
        with Wdf(str(path)) as w:
            start = datetime.datetime(1601, 1, 1, tzinfo=datetime.timezone.utc) + datetime.timedelta(microseconds=w.hdr.time_start_ / 10)
            try:
                local_time = start.astimezone(ZoneInfo('America/Chicago')).isoformat()
            except ZoneInfoNotFoundError:
                local_time = None
                meta['optional_metadata_note'] = 'Chicago time-zone database unavailable; UTC acquisition time retained.'
            meta.update(acquired_utc=start.isoformat(), acquired_local=local_time,
                        duration_s=(w.hdr.time_end_ - w.hdr.time_start_) / 1e7,
                        uuid=w.hdr.uuid_.hex(), accumulations=w.hdr.naccum)
    except (ImportError, AttributeError, OSError, ValueError, OverflowError) as exc:
        meta['optional_metadata_note'] = f'Extended WDF metadata unavailable: {type(exc).__name__}'
    return meta


def analyze_wdf(path: Path, name=None, sample=None):
    """Analyze one supported map. No image block or naming convention is required.

    Metadata lives at the dictionary's top level for compatibility with the
    figure loaders. Arrays include x, raw_mean, mean, corrected, normalized,
    map_index and clean. pixels is a native-pixel pandas DataFrame.
    """
    path = Path(path)
    if not path.is_file():
        raise ValueError(f'WDF file not found: {path}')
    try:
        reader = WDFReader(str(path))
    except Exception as exc:
        raise ValueError(f'Unable to read WDF {path.name}: {exc}') from exc
    try:
        if reader.xlist_unit != UnitType.RamanShift or reader.xlist_type != DataType.Frequency:
            raise ValueError('The WDF spectral axis must be Raman shift in cm^-1, not wavelength, energy, or absolute wavenumber.')
        if reader.spectral_unit != UnitType.Counts:
            raise ValueError('WDF intensities must be stored as counts because the cosmic-screening and noise criteria use count-scale thresholds. Use an original count-valued WDF, not normalized or arbitrary-intensity spectra.')
        nx, ny, step, xp, yp, ix, iy = _geometry(reader)
        order = np.argsort(reader.xdata)
        x = _axis(np.asarray(reader.xdata, dtype=float)[order])
        weights = {key: band_weights(x, *definition) for key, definition in BANDS.items()}
        exact = {key: exact_band_weights(x, definition)[0] for key, definition in BOUNDED_BANDS.items()}
        _window(x, (1500, 1700), minimum=8)
        if np.isclose(step, .1, atol=1e-6):
            _window(x, (765, 785))
            _window(x, (865, 885))
        source = np.asarray(reader.spectra, dtype=float)
        if source.size != reader.count * len(x):
            raise ValueError('Spectral data size does not agree with map count and axis length.')
        raw = source.reshape(reader.count, len(x))[:, order]
        clean, noise, spikes = despike(raw, x)
        eligible = ~spikes[:, (x >= 700) & (x <= 1500)].any(axis=1)
        if not eligible.any():
            raise ValueError('No eligible spectra remain after conservative cosmic-event screening.')
        canonical = str(name) if name is not None else re.sub(r'\(\d+\)$', '', path.stem)
        if sample is None:
            match = re.match(r'^(E\d+|PP|PE)(?=$|[-_\s(])', canonical, re.I)
            sample = match.group(1).upper() if match else canonical
        d = dict(name=canonical, file=path.name, path=str(path.resolve()),
                 sha256=hashlib.sha256(path.read_bytes()).hexdigest(), sample=str(sample),
                 nx=nx, ny=ny, n=int(reader.count), step=round(step, 6),
                 width=(nx - 1) * step, height=(ny - 1) * step,
                 xmin=float(xp.min()), ymin=float(yp.min()), xmax=float(xp.max()), ymax=float(yp.max()),
                 complete=True, x=x, raw_mean=raw.mean(axis=0), mean=clean[eligible].mean(axis=0), clean=clean)
        d.update(_metadata(reader, path))
        df = pd.DataFrame(dict(ix=ix, iy=iy, stage_x_um=xp, stage_y_um=yp,
                              noise_counts=noise, despiked_channels=spikes.sum(axis=1), eligible=eligible))
        for key, w in weights.items():
            df['A' + key], df['V' + key], df['raw_A' + key] = clean @ w, noise ** 2 * (w @ w), raw @ w
        df = _add_ratios(df)
        a, b = clean @ exact['PE1295'], clean @ exact['PP809']
        ap, bp = np.maximum(a, 0), np.maximum(b, 0)
        den = ap + bp
        sig = noise * np.linalg.norm(exact['PE1295'] + exact['PP809'])
        df['PE1295_common'], df['PP809_common'], df['combined_snr'] = a, b, den / sig
        df['pe_bounded'] = np.divide(ap, den, out=np.full(len(df), np.nan), where=den > 0)
        df['valid_bounded'] = eligible & (den > 0) & (den / sig >= 5)
        if not df.valid_bounded.any():
            raise ValueError('No pixels have usable combined PE/PP marker signal (SNR >= 5); a chemical map cannot be reported.')
        # No PP-only threshold is applied to the displayed bounded chemical map.
        image = np.full((ny, nx), np.nan)
        image[iy, ix] = np.where(df.valid_bounded, df.pe_bounded, np.nan)
        d['pixels'], d['map_index'] = df, image
        d['corrected'] = d['mean'] - als(d['mean'])
        peak = d['corrected'][(x >= 1325) & (x <= 1335)].max()
        d['normalized'] = d['corrected'] / peak if peak > 0 else np.full(len(x), np.nan)
        d['regional_spectra'] = {}
        return d
    finally:
        handle = getattr(reader, 'file_obj', None)
        if handle is not None:
            handle.close()


def summary(d):
    """Native-pixel map summary; intentionally omits legacy binned metrics."""
    df = d['pixels']
    row = {key: d.get(key) for key in SUMMARY_KEYS}
    row.update(despiked_spectra=int((df.despiked_channels > 0).sum()),
               despiked_channels=int(df.despiked_channels.sum()), n_eligible=int(df.eligible.sum()),
               cosmic_excluded_spectra=int((~df.eligible).sum()), noise_median=float(df.noise_counts.median()),
               pp_snr_median=float(np.median(df.A1330 / np.sqrt(df.V1330))))
    good = df.valid_r1 & df.valid_r2
    row['cross_marker_rho_native'] = _rho(df.loc[good, 'r1'], df.loc[good, 'r2'])
    row['cross_marker_native_n'] = int(good.sum())
    for key, (a, b) in RATIOS.items():
        values = df.loc[df['valid_' + key], key]
        row.update({key + '_pooled': pool(df, key), key + '_n_valid': len(values),
                    key + '_median': float(values.median()), key + '_q10': float(values.quantile(.1)),
                    key + '_q90': float(values.quantile(.9)),
                    key + '_raw_pooled': _divide(df['raw_A' + a].sum(), df['raw_A' + b].sum()),
                    key + '_detected_pct': 100 * float(df['detected_' + key].mean())})
    eligible = df.eligible.to_numpy()
    a, b = np.maximum(df.PE1295_common.to_numpy(), 0), np.maximum(df.PP809_common.to_numpy(), 0)
    row.update(valid_bounded_n=int(df.valid_bounded.sum()), valid_bounded_pct=100 * float(df.valid_bounded.mean()),
               bounded_pool=_divide(a[eligible].sum(), (a[eligible] + b[eligible]).sum()),
               old_r1_discard_pct=100 * float((~df.valid_r1).mean()))
    return row


def regional_rows(d):
    """Composition-marker quartiles and PP/PE regional ratios for native fine maps.

    Applies to any 0.1 µm map, including generically named files. Quartiles alone
    do not validate chemically distinct regions. Returns an empty list for
    coarser maps or insufficient eligible composition-marker signal.
    """
    if not np.isclose(d['step'], .1, atol=1e-6):
        return []
    df, x, clean = d['pixels'], d['x'], d['clean']
    ca, cb = np.maximum(df.A1295.to_numpy(), 0), np.maximum(df.A1330.to_numpy(), 0)
    cw = band_weights(x, *BANDS['1295']) + band_weights(x, *BANDS['1330'])
    cs = np.divide(ca, ca + cb, out=np.full(len(ca), np.nan), where=(ca + cb) > 0)
    good = df.eligible.to_numpy() & ((ca + cb) / (df.noise_counts.to_numpy() * np.linalg.norm(cw)) >= 5)
    if not good.any():
        return []
    q = np.quantile(cs[good], [.25, .75])
    # Tied scores can overlap; retain the exact method and expose cutoffs/counts.
    masks = {'lower_PE_marker': good & (cs <= q[0]), 'higher_PE_marker': good & (cs >= q[1])}
    rows = []
    for label, mask in masks.items():
        spectrum = clean[mask].mean(axis=0)
        d['regional_spectra'][label] = spectrum
        areas = {key: float(spectrum @ band_weights(x, *definition)) for key, definition in BANDS.items()}
        base = {'left': [765, 785], 'right': [865, 885]}
        p809 = spectrum @ exact_band_weights(x, dict(base, integration=[800, 819]))[0]
        p841 = spectrum @ exact_band_weights(x, dict(base, integration=[833, 849]))[0]
        rows.append(dict(sample=d['sample'], region=label, n_pixels=int(mask.sum()),
                         selection_score_low_quartile=float(q[0]), selection_score_high_quartile=float(q[1]),
                         r1_pooled=_divide(areas['1295'], areas['1330']), r2_pooled=_divide(areas['1128'], areas['974']),
                         PP809_841_local=_divide(areas['809'], areas['841']), PP809_841_shared_baseline=_divide(p809, p841),
                         PE1416_1295=_divide(areas['1416'], areas['1295']), PE1416_area=areas['1416'],
                         PE1295_area=areas['1295'], PP809_area=areas['809'], PP841_area=areas['841']))
    return rows


def _fine_report_maps(data):
    for sample in (f'E{i}' for i in range(6)):
        d = data.get(sample + '-small')
        if d is not None and d['sample'] == sample and np.isclose(d['step'], .1, atol=1e-6):
            yield sample, d


def matched(data):
    """Recorded-stage-footprint comparisons for available E0–E5 report roles."""
    columns = ['sample', 'coarse', 'fine', 'n_coarse_roi', 'n_fine', 'fine_span_um2']
    for key in RATIOS:
        columns += [key + suffix for suffix in ['_coarse_full', '_coarse_roi', '_fine', '_roi_vs_full_pct', '_fine_vs_roi_pct']]
    rows = []
    for sample, fine in _fine_report_maps(data):
        for d in data.values():
            if d['sample'] != sample or not np.isclose(d['step'], 1, atol=1e-6):
                continue
            p = d['pixels']
            good = ((p.stage_x_um >= fine['xmin'] - .05) & (p.stage_x_um <= fine['xmax'] + .05)
                    & (p.stage_y_um >= fine['ymin'] - .05) & (p.stage_y_um <= fine['ymax'] + .05))
            if not good.any():
                continue
            row = dict(sample=sample, coarse=d['name'], fine=fine['name'], n_coarse_roi=int(good.sum()),
                       n_fine=fine['n'], fine_span_um2=fine['width'] * fine['height'])
            for key in RATIOS:
                full, roi, f = pool(p, key), pool(p[good], key), pool(fine['pixels'], key)
                row.update({key + '_coarse_full': full, key + '_coarse_roi': roi, key + '_fine': f,
                            key + '_roi_vs_full_pct': 100 * (_divide(roi, full) - 1),
                            key + '_fine_vs_roi_pct': 100 * (_divide(f, roi) - 1)})
            rows.append(row)
    return pd.DataFrame(rows, columns=columns)


def subsample(data):
    """Point-thinning over every 1/2 µm grid offset in available report fine maps."""
    rows = []
    for sample, d in _fine_report_maps(data):
        p = d['pixels']
        for factor in (10, 20):
            for oy in range(factor):
                for ox in range(factor):
                    group = p[(p.ix % factor == ox) & (p.iy % factor == oy)]
                    rows.append(dict(sample=sample, step_um=.1 * factor, x_offset_pixels=ox,
                                     y_offset_pixels=oy, n=len(group), r1=pool(group, 'r1'), r2=pool(group, 'r2')))
    return pd.DataFrame(rows, columns=['sample', 'step_um', 'x_offset_pixels', 'y_offset_pixels', 'n', 'r1', 'r2'])
