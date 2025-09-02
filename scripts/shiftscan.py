#!/usr/bin/env python3

import argparse, sys, os, gc
from pathlib import Path
import pandas as pd
import numpy as np
from DSF_functions import (read_default_metadata, roche_import, biorad_import,roche_import_platewise, biorad_import_platewise,concatenate_files,default_analysis)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="DSF Data Processing Script")
    parser.add_argument("-i", "--input_dir", required=True, help="Full file path to directory with input files")
    parser.add_argument("-m", "--metadata", required=True, help="Full file path to metadata file")
    parser.add_argument("-o", "--output_dir", required=True, help="Full path to desired output directory")
    parser.add_argument("-p", "--processors", type=int, default=4, help="No. of processors you want to use")
    parser.add_argument("-f", "--file_origin", choices=['RocheLightCycler', 'BioRadOpticonMonitor'], default='RocheLightCycler', help="Instrument used to generate raw output.")
    parser.add_argument("-d", "--delimiter", default="\t", help="Separator character in raw data input file")
    parser.add_argument("-c", "--control_cols", default="1,2", help="The column numbers of your controls (comma-separated) or path to a well mapping file.")
    parser.add_argument("-s", "--smoothing_factor", type=float, default=0.0005, help="Desired smoothing factor")
    parser.add_argument("-n", "--normalization", choices=['y', 'n'], default="y", help="Should data be normalized (y/n)")
    parser.add_argument("-x", "--failed_control_wells", type=int, default=8, help="The number of controls allowed to fail before a plate is failed")
    parser.add_argument("-t", "--ctrl_tm_cutoff", type=float, default=1.5, help="Maximum z-score of ctrl curve melting temps")
    parser.add_argument("-a", "--ctrl_amp_cutoff", type=float, default=2.0, help="Maximum z-score of ctrl curve amplitude")
    parser.add_argument("-u", "--max_amp_cutoff", type=float, default=6.0, help="Maximum relative amplitude of curves allowed")
    parser.add_argument("-l", "--min_amp_cutoff", type=float, default=0.2, help="Minimum relative amplitude of curves allowed")
    parser.add_argument("--disk", action='store_true', help="Use disk-based (plate-wise) processing for lower RAM usage. Default is RAM-based.")
    parser.add_argument("--only_tm", action='store_true', help="Flag to enable only Tm calling mode, skipping QC and reports.")
    parser.add_argument("--dose_response", action='store_true', help="Flag to enable dose response mode")

    args = parser.parse_args()

    # Get output folder prepped
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output will be saved to: {output_dir.resolve()}")

    #Read in metadata
    dose_response = ("yes" if args.dose_response else "no")
    map_raw_df = read_default_metadata(args.metadata, args.delimiter, dose_response)

    # Define a dictionary of parameters to pass to the analysis function
    analysis_params = {
        'control_cols': args.control_cols,
        'delimiter': args.delimiter,
        'normalization': args.normalization,
        'smoothing_factor': args.smoothing_factor,
        'processors': args.processors,
        'only_tm': args.only_tm,
        'dose_response':args.dose_response,
        'failed_control_wells': args.failed_control_wells,
        'ctrl_tm_cutoff': args.ctrl_tm_cutoff,
        'ctrl_amp_cutoff': args.ctrl_amp_cutoff,
        'min_amp_cutoff': args.min_amp_cutoff,
        'max_amp_cutoff': args.max_amp_cutoff
    }

    if args.disk:
        print("Starting DISK mode (low RAM, plate-by-plate)")
        plates = [file for file in Path(args.input_dir).glob('*') if file.suffix in ['.txt', '.csv']]
        if not plates:
            print(f"No .txt or .csv files found in specified input folder: {args.input_dir}")
            sys.exit()
        
        for i, plate_file in enumerate(plates):
            print(f"Processing Plate {i+1}/{len(plates)}: {plate_file.name}")
            
            if args.file_origin == 'RocheLightCycler':
                semifinal_df = roche_import_platewise(plate_file, args.delimiter)
            else:
                semifinal_df = biorad_import_platewise(plate_file, args.delimiter)


            tmp_map = map_raw_df[['Assay_Plate', 'Source_Plate']].drop_duplicates()
            semifinal_df2 = pd.merge(semifinal_df, tmp_map, on='Assay_Plate', how='inner')
            if args.dose_response:
                final_df = pd.merge(semifinal_df2, map_raw_df[['Unique_key', 'Compound', 'Fraction','Concentration','Replicate']].drop_duplicates(), on='Unique_key', how='left')
            else: 
                final_df = pd.merge(semifinal_df2, map_raw_df[['Unique_key', 'Compound', 'Fraction']].drop_duplicates(), on='Unique_key', how='left')
                
            # RUN THE ANALYSIS
            results = default_analysis(final_df, **analysis_params)
            
            # Write intermediate files
            if args.only_tm:
                results['Tm_df'].to_csv(output_dir / f"Only_Tm_values_plate_{i+1}.txt", sep="\t", index=False)
                results['curves_df'].to_csv(output_dir / f"Only_Tm_curves_plate_{i+1}.txt", sep="\t", index=False)
            elif args.dose_response:
                results['Tm_df'].to_csv(output_dir / f"DR_final_results_{i+1}.txt", sep="\t", index=False)
                results['curves_df'].to_csv(output_dir / f"DR_final_curves_{i+1}.txt", sep="\t", index=False)
            else:
                results['Tm_df'].to_csv(output_dir / f"Final_results_plate_{i+1}.txt", sep="\t", index=False)
                results['curves_df'].to_csv(output_dir / f"Final_curves_plate_{i+1}.txt", sep="\t", index=False)
                results['plate_report_df'].to_csv(output_dir / f"Plate_report_plate_{i+1}.txt", sep="\t", index=False)
                results['problems_df'].to_csv(output_dir / f"Potential_problems_plate_{i+1}.txt", sep="\t", index=False)
        
        if args.only_tm:
            print("Concatenating all plate data")
            tm_files = [file for file in os.listdir(output_dir) if (file.endswith(".txt"))&file.startswith("Only_Tm_values")]
            curves_files = [file for file in os.listdir(output_dir) if (file.endswith(".txt"))&file.startswith("Only_Tm_curves")]
            concatenate_files(tm_files, output_dir, 'Only_Tm_values.txt')
            concatenate_files(curves_files, output_dir, 'Only_Tm_curves.txt')
            for file_path in tm_files:
                os.remove(os.path.join(output_dir,file_path))
            for file_path in curves_files:
                os.remove(os.path.join(output_dir,file_path))
        elif args.dose_response:
            print("Concatenating all plate data")
            tm_files = [file for file in os.listdir(output_dir) if (file.endswith(".txt"))&file.startswith("DR_final_results_")]
            curves_files = [file for file in os.listdir(output_dir) if (file.endswith(".txt"))&file.startswith("DR_final_curves_")]
            concatenate_files(tm_files, output_dir, 'DR_final_results.txt')
            concatenate_files(curves_files, output_dir, 'DR_final_curves.txt')
            for file_path in tm_files:
                os.remove(os.path.join(output_dir,file_path))
            for file_path in curves_files:
                os.remove(os.path.join(output_dir,file_path))
        else:
            print("Concatenating all plate data")
            Tm_files = [file for file in os.listdir(output_dir) if (file.endswith(".txt"))&file.startswith("Final_results_plate_")]
            curve_files = [file for file in os.listdir(output_dir) if (file.endswith(".txt"))&file.startswith("Final_curves_plate_")]
            plate_report_files = [file for file in os.listdir(output_dir) if (file.endswith(".txt"))&file.startswith("Plate_report_plate_")]
            well_error_count_files = [file for file in os.listdir(output_dir) if (file.endswith(".txt"))&file.startswith("Potential_problems_plate_")]
            concatenate_files(Tm_files, output_dir, 'Final_results.txt')
            concatenate_files(curve_files, output_dir, 'Final_curves.txt')
            concatenate_files(plate_report_files, output_dir, 'Plate_report.txt')
            concatenate_files(well_error_count_files, output_dir, 'Potential_problems.txt')
            for file_path in Tm_files:
                os.remove(os.path.join(output_dir,file_path))
            for file_path in curve_files:
                os.remove(os.path.join(output_dir,file_path))
            for file_path in plate_report_files:
                os.remove(os.path.join(output_dir,file_path))
            for file_path in well_error_count_files:
                os.remove(os.path.join(output_dir,file_path))
    else:
        print("Starting RAM mode (fast, all-in-memory)")
        files = [file for file in Path(args.input_dir).glob('*') if file.suffix in ['.txt', '.csv']]
        if not files:
            print(f"No .txt or .csv files found in specified input folder: {args.input_dir}")
            sys.exit()
            
        print(f"Importing data from {len(files)} files...")
        if args.file_origin == 'RocheLightCycler':
            semifinal_df = roche_import(files, args.delimiter)
        else:
            semifinal_df = biorad_import(files, args.delimiter)

        tmp_map = map_raw_df[['Assay_Plate', 'Source_Plate']].drop_duplicates()
        semifinal_df2 = pd.merge(semifinal_df, tmp_map, on='Assay_Plate', how='inner')

        if args.dose_response:
            final_df = pd.merge(semifinal_df2, map_raw_df[['Unique_key', 'Compound', 'Fraction','Concentration','Replicate']].drop_duplicates(), on='Unique_key', how='left')
        else: 
            final_df = pd.merge(semifinal_df2, map_raw_df[['Unique_key', 'Compound', 'Fraction']].drop_duplicates(), on='Unique_key', how='left')
        
        del semifinal_df, semifinal_df2, tmp_map, map_raw_df
        gc.collect()
        
        # RUN THE ANALYSIS
        results = default_analysis(final_df, **analysis_params)
        
        print("Writing final output files...")
        if args.only_tm:
            results['Tm_df'].to_csv(output_dir / "Only_Tm_values.txt", sep="\t", index=False)
            results['curves_df'].to_csv(output_dir / "Only_Tm_curves.txt", sep="\t", index=False)
        elif args.dose_response:
            results['Tm_df'].to_csv(output_dir / "DR_final_results.txt", sep="\t", index=False)
            results['curves_df'].to_csv(output_dir / "DR_final_curves.txt", sep="\t", index=False)
        else:
            results['Tm_df'].to_csv(output_dir / "Final_results.txt", sep="\t", index=False)
            results['curves_df'].to_csv(output_dir / "Final_curves.txt", sep="\t", index=False)
            if results['plate_report_df'] is not None:
                results['plate_report_df'].to_csv(output_dir / "Plate_report.txt", sep="\t", index=False)
            if results['problems_df'] is not None:
                results['problems_df'].to_csv(output_dir / "Potential_problems.txt", sep="\t", index=False)

    print("\nAnalysis complete!")