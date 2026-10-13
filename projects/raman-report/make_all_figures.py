"""Create report and supplementary figures from bundled or separately analyzed data.

Examples:
    python make_all_figures.py
    python make_all_figures.py --data-dir analysis/data --output-dir analysis/figures
Individual figure scripts honor RAMAN_DATA_DIR and RAMAN_OUTPUT_DIR as well.
"""
import argparse
import importlib
import os
from pathlib import Path
SCRIPTS = ['figure_2_regional_spectra', 'figure_3A_maps_E0_E2', 'figure_3B_maps_E3_E5', 'figure_3C_sampling_effects', 'figure_3D_E4_interface_width', 'supplement_PP_PE_control_maps', 'supplement_Raman_workflow']
def main(argv=None):
    parser=argparse.ArgumentParser(description='Plot all Raman report figures. Custom data must contain the complete report data layout; bundled results are never substituted for missing inputs.')
    parser.add_argument('--data-dir',type=Path,help='Directory containing map_metadata.json, maps/, report CSVs and E4 analysis subdirectories (default: RAMAN_DATA_DIR or bundled data/)')
    parser.add_argument('--output-dir',type=Path,help='Destination for PNG, PDF and SVG files (default: RAMAN_OUTPUT_DIR or output/)')
    args=parser.parse_args(argv)
    if args.data_dir is not None: os.environ['RAMAN_DATA_DIR']=str(args.data_dir.expanduser().resolve())
    if args.output_dir is not None: os.environ['RAMAN_OUTPUT_DIR']=str(args.output_dir.expanduser().resolve())
    # Set paths before importing figure modules, which load their shared data.
    for name in SCRIPTS:
        importlib.import_module(name).main()
if __name__ == '__main__':
    main()
