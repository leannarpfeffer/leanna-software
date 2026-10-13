# Raman report analysis and plotting

Thirteen Python scripts for analyzing raw PP/PE Raman maps and reproducing the E-series Raman figures. This public directory contains **code and reproducibility metadata only**. Raw WDF acquisitions, derived NPZ/CSV datasets, and generated figures must be supplied or generated separately.

## Install and analyze raw maps

Use Python 3.11 or newer. From the repository root:

```bash
cd projects/raman-report
python -m pip install -r requirements.txt
python run_wdf.py "path/to/map.wdf" --output "new_results"
```

Other input options:

```bash
python run_wdf.py --input-dir "raw_data_folder" --output "new_results"
python run_wdf.py
```

The last command opens a file chooser requiring Tkinter and a desktop session. CLI input works without a desktop. On Windows, `py` can replace `python`.

Folders are searched recursively; multiple file/folder arguments and uppercase `.WDF` extensions work. Extract ZIP archives first. Byte-identical duplicates are analyzed once. Generic mode labels maps by filename without inferring E0–E5 pairings. Raw inputs are read only. Output must be a new or empty folder; the default is `runs/<UTC timestamp>`.

**Supported maps:** complete regular 2D grids with equal X/Y pitch, micrometre coordinates, count-valued intensity, and a Raman-shift axis in cm⁻¹. Required bands, baseline windows, and noise coverage span approximately 765–1700 cm⁻¹. Rectangular fields are supported; partial maps, line scans, single spectra, irregular grids, and unsupported units are rejected. Microscope images are unnecessary. This is a PP/PE marker method, not a general chemical-identification tool.

Outputs include native-pixel CSV tables, NPZ maps/mean spectra, summary and regional tables, source/acquisition metadata, PNG/PDF/SVG figures, and `run_manifest.json`. Check that the manifest says `complete` before using a run. Generated provenance records contain local input paths and remain separate from this public source directory.

## Reproduce the report figures

With the original raw acquisitions:

```bash
python run_wdf.py --report --input-dir "raw_wdf_folder" --output "report_results"
```

[report_inputs.json](report_inputs.json) identifies 20 required maps by SHA-256: E0–E5 at 2, 1, and 0.1 µm spacing, plus PP/PE controls. Two extra E4 maps are optional. Renamed byte-identical files work; missing or changed required files cause an error. Use generic mode for new acquisitions.

Report mode freshly computes the map metrics, matched footprints, grid-offset sensitivity, and original E4 boundary/regional-spectra recipe before plotting. That E4 recipe is hash-guarded and specific to one acquisition, not an automatic boundary detector for new specimens. The workflow does not create a Word document, reproduce original Figure 1 or DSC/XRD Figures 4A–4E, or rerun the separate significance tests.

To replot a completed report run:

```bash
python make_all_figures.py --data-dir "report_results/data" --output-dir "replotted_figures"
```

Alternatively, point `--data-dir` to the complete `data` directory from the separately supplied original analysis package. No such datasets are published here. Running without an external data path fails unless that directory has been supplied locally. Individual figure scripts accept `RAMAN_DATA_DIR` and `RAMAN_OUTPUT_DIR`; their shared loader requires the complete report data layout. Missing inputs never trigger a fallback to another dataset.

The scripts are unchanged from the standalone package. References to “included data” or `README.txt` in original docstrings refer to that package; use these instructions for this public checkout.

## Scripts

| Script | Role |
|---|---|
| [run_wdf.py](run_wdf.py) | Input selection, raw workflow, exports, and run manifest |
| [wdf_analysis.py](wdf_analysis.py) | WDF validation, preprocessing, band metrics, regions, and sampling comparisons |
| [e4_report_analysis.py](e4_report_analysis.py) | Original E4 boundary and empirical spectral-mixture recipe |
| [raw_figures.py](raw_figures.py) | Generic native maps, mean spectra, and map distributions |
| [plot_common.py](plot_common.py) | Shared report-data loading, validation, styling, and exports |
| [make_all_figures.py](make_all_figures.py) | Generate all five report and two supplementary figures |
| [figure_2_regional_spectra.py](figure_2_regional_spectra.py) | Figure 2: references, regional PP ratios, and E4 spectra |
| [figure_3A_maps_E0_E2.py](figure_3A_maps_E0_E2.py) | Figure 3A: native E0–E2 maps |
| [figure_3B_maps_E3_E5.py](figure_3B_maps_E3_E5.py) | Figure 3B: native E3–E5 maps |
| [figure_3C_sampling_effects.py](figure_3C_sampling_effects.py) | Figure 3C: field selection and grid-offset sensitivity |
| [figure_3D_E4_interface_width.py](figure_3D_E4_interface_width.py) | Figure 3D: E4 profiles and apparent 10–90% widths |
| [supplement_PP_PE_control_maps.py](supplement_PP_PE_control_maps.py) | Neat PP/PE reference maps |
| [supplement_Raman_workflow.py](supplement_Raman_workflow.py) | Raman-focused study diagram, distinct from original Figure 1 |

## Interpretation and reproducibility

Native maps are not spatially smoothed or rebinned. Cosmic-event screening and combined marker signal-to-noise determine valid pixels; weak PP signal alone does not exclude a PE-rich region. The displayed marker is:

```text
Q = max(A1295, 0) / [max(A1295, 0) + max(A809, 0)]
```

**Q and the other spectral ratios are not calibrated PE mass fractions or crystallinity.** Local baseline definitions are explicit in the source. Pooled ratios use summed eligible band areas; ALS mean-spectrum correction is for display. Raw and cleaned means differ through both replacement and screening. Fine-map quartiles use a separate 1295/1330-based marker; tied quartiles can overlap.

Sampling pitch is not optical resolution. Stage-coordinate matching is not image registration. Grid-offset percentiles are not confidence intervals. E4 widths are apparent Raman-marker transitions, not deconvolved intrinsic interphase thicknesses; seven profiles along one boundary are not independent specimens.

[VALIDATION.txt](VALIDATION.txt) records prior raw-data reproduction and input checks. [tested_versions.txt](tested_versions.txt) records the tested environment; it is not a pip requirements file. Fonts, dependency versions, and optimizer/BLAS roundoff may cause small visual or numerical differences. The validation datasets are not distributed here.

This directory follows the MIT license at the repository root.
