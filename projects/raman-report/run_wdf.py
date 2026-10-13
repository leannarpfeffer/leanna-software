"""Read raw Renishaw WDF maps, analyze native pixels, and save figures/tables.

No arguments opens a file chooser. --report reproduces the report from its
checksum-identified inputs; ordinary mode accepts new compatible PP/PE maps.
"""
import argparse
import datetime
import hashlib
import importlib.metadata
import io
import json
import os
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def write(path, content):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    blob = content.encode('utf-8') if isinstance(content, str) else content
    temp = path.with_name(path.name + '.tmp')
    with temp.open('wb') as stream:
        stream.write(blob)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temp, path)
    if path.read_bytes() != blob:
        raise OSError(f'Write verification failed: {path}')


def json_value(value):
    # Strict JSON uses null for non-finite values, including absent ratios.
    if isinstance(value, dict):
        return {str(k): json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_value(v) for v in value]
    if hasattr(value, 'item'):
        return json_value(value.item())
    if isinstance(value, float) and not (-float('inf') < value < float('inf')):
        return None
    return value


def write_json(path, value):
    write(path, json.dumps(json_value(value), indent=2, allow_nan=False) + '\n')


def sha256(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def choose_files():
    try:
        import tkinter as tk
        from tkinter import filedialog
        window = tk.Tk()
        window.withdraw()
        try:
            selected = filedialog.askopenfilenames(
                title='Select raw Raman WDF maps',
                filetypes=[('Renishaw WDF', '*.wdf *.WDF'), ('All files', '*.*')])
        finally:
            window.destroy()
    except (ImportError, RuntimeError, OSError) as error:
        raise ValueError('File chooser unavailable. Use: python run_wdf.py "path/to/map.wdf" --output "new_results"') from error
    except Exception as error:
        # Includes TclError on headless hosts; CLI does not need Tk.
        raise ValueError('File chooser unavailable. Supply WDF paths or --input-dir on the command line.') from error
    if not selected:
        raise ValueError('No files selected; no analysis was started.')
    return list(map(Path, selected))


def discover(paths):
    found = set()
    for source in paths:
        source = source.expanduser().resolve()
        if source.is_dir():
            found.update(p.resolve() for p in source.rglob('*')
                         if p.is_file() and p.suffix.lower() == '.wdf')
        elif source.is_file() and source.suffix.lower() == '.wdf':
            found.add(source)
        else:
            raise ValueError(f'Not a WDF file or input directory: {source}')
    if not found:
        raise ValueError('No .wdf files found. Extract ZIP archives before running.')
    return sorted(found, key=str)


def select_inputs(paths, report):
    by_hash = {}
    for path in paths:
        by_hash.setdefault(sha256(path), path)
    if report:
        config = json.loads((ROOT / 'report_inputs.json').read_text(encoding='utf-8'))
        missing = [r['source_file'] for r in config['inputs']
                   if r['required'] and r['sha256'] not in by_hash]
        if missing:
            raise ValueError('Report reproduction requires the original files identified by SHA-256. Missing or changed: '
                             + ', '.join(missing) + '. Select their enclosing folder; use ordinary mode for new datasets.')
        return [(r['name'], r['sample'], by_hash[r['sha256']], r['sha256'])
                for r in config['inputs'] if r['sha256'] in by_hash]
    selected, used = [], set()
    for digest, path in by_hash.items():
        name = re.sub(r'[^A-Za-z0-9_.-]+', '_', path.stem).strip('._') or 'map'
        if name.lower() in used:
            name += '_' + digest[:10]
        used.add(name.lower())
        # Generic labels come from the selected filename, with no sample grouping.
        selected.append((name, name, path, digest))
    return selected


def export_map(d, data_dir):
    import gzip
    import numpy as np
    name = d['name']
    buffer = io.BytesIO()
    arrays = {k: d[k] for k in ('x', 'raw_mean', 'mean', 'corrected', 'map_index')}
    arrays.update({'region_' + k: v for k, v in d['regional_spectra'].items()})
    np.savez_compressed(buffer, **arrays)
    write(data_dir / 'maps' / (name + '.npz'), buffer.getvalue())
    write(data_dir / 'pixels' / (name + '_pixels.csv.gz'),
          gzip.compress(d['pixels'].to_csv(index=False).encode('utf-8'), mtime=0))
    excluded = {'pixels', 'clean', 'regional_spectra'} | set(arrays) | {'normalized'}
    return {k: v for k, v in d.items() if k not in excluded}


def process(selected, output, report):
    import pandas as pd
    import wdf_analysis as analysis
    from raw_figures import plot_raw_maps
    data_dir = output / 'data'
    data, metadata, regional = {}, {}, []
    manifest = dict(status='running', started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                    mode='report' if report else 'generic', inputs=[], software={}, outputs=[])
    for package in ('numpy', 'pandas', 'scipy', 'matplotlib', 'renishawWiRE', 'renishaw-wdf'):
        try:
            manifest['software'][package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            manifest['software'][package] = 'not installed (optional metadata reader)'
    manifest['software']['python'] = sys.version
    manifest['interpretation'] = 'Marker responses are uncalibrated. Sampling step is not optical resolution. E4 widths are apparent Raman-marker transition widths.'
    manifest['method'] = dict(native_pixels=True, spatial_smoothing=False,
                              bounded_combined_snr_minimum=5,
                              cosmic_affected_fingerprint_spectra_excluded=True,
                              regional_quartile_score='positive A1295/(positive A1295+positive A1330)')
    try:
        write_json(output / 'run_manifest.json', manifest)
        for index, (name, sample, source, digest) in enumerate(selected, 1):
            print(f'[{index}/{len(selected)}] Reading {source.name}', flush=True)
            d = analysis.analyze_wdf(source, name=name, sample=sample)
            if d['sha256'] != digest:
                raise ValueError(f'Input changed during analysis: {source}')
            regional.extend(analysis.regional_rows(d))
            d.pop('clean', None)  # Release full spectra after regional means.
            metadata[name] = export_map(d, data_dir)
            data[name] = d
            manifest['inputs'].append(dict(name=name, source=str(source), sha256=digest, spectra=d['n']))
        write_json(data_dir / 'map_metadata.json', metadata)
        write_json(data_dir / 'source_manifest.json', manifest['inputs'])
        write_json(data_dir / 'band_definitions.json', dict(local_least_squares=analysis.BANDS,
                                                          bounded_exact=analysis.BOUNDED_BANDS))
        write(data_dir / 'summary.csv', pd.DataFrame([analysis.summary(d) for d in data.values()]).to_csv(index=False))
        write(data_dir / 'regional_structure.csv', pd.DataFrame(regional, columns=analysis.REGIONAL_COLUMNS).to_csv(index=False))
        if report:
            write(data_dir / 'matched_footprints.csv', analysis.matched(data).to_csv(index=False))
            write(data_dir / 'sampling_offsets.csv', analysis.subsample(data).to_csv(index=False))
            from e4_report_analysis import run_e4_report
            manifest['e4_analysis'] = run_e4_report(Path(data['E4-small']['path']), data_dir)
        plot_raw_maps(data, output / 'raw_figures')
        if report:
            subprocess.run([sys.executable, str(ROOT / 'make_all_figures.py'),
                            '--data-dir', str(data_dir), '--output-dir', str(output / 'figures')], check=True)
        manifest['status'] = 'complete'
        manifest['spectra_total'] = sum(d['n'] for d in data.values())
        manifest['outputs'] = [str(p.relative_to(output)) for p in sorted(output.rglob('*')) if p.is_file()]
    except Exception as error:
        manifest.update(status='failed', error=f'{type(error).__name__}: {error}')
        raise
    finally:
        manifest['finished_utc'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        write_json(output / 'run_manifest.json', manifest)
    return manifest


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('inputs', nargs='*', type=Path, help='WDF files or folders (searched recursively)')
    parser.add_argument('--input-dir', type=Path, help='Additional folder containing raw WDF files')
    parser.add_argument('--output', type=Path, help='New or empty result folder; default: runs/<UTC timestamp>')
    parser.add_argument('--report', action='store_true', help='Reproduce all seven report/supplement figures from the original checksum-identified inputs')
    args = parser.parse_args(argv)
    output = None
    started = False
    try:
        paths = args.inputs + ([args.input_dir] if args.input_dir else [])
        selected = select_inputs(discover(paths or choose_files()), args.report)
        stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S_%fZ')
        output = (args.output or ROOT / 'runs' / stamp).expanduser().resolve()
        if output.exists() and (not output.is_dir() or any(output.iterdir())):
            raise ValueError(f'Output must be a new or empty directory: {output}')
        output.mkdir(parents=True, exist_ok=True)
        started = True
        manifest = process(selected, output, args.report)
    except Exception as error:
        print(f'Cannot complete Raman analysis: {error}', file=sys.stderr)
        if isinstance(error, ImportError):
            print('Install dependencies: python -m pip install -r requirements.txt', file=sys.stderr)
        if started and (output / 'run_manifest.json').is_file():
            print(f'See {output / "run_manifest.json"}; this run is incomplete.', file=sys.stderr)
        return 1
    print(f'Completed {len(selected)} maps / {manifest["spectra_total"]:,} spectra. Results: {output}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
