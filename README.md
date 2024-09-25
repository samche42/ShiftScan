# ShiftScan
Analysis and visualization of high-throughput DSF data

#### Installation and usage video
A very brief video of ShiftScan installation (conda already installed) is available on YouTube: https://youtu.be/shN0vUHEGOk

#### Step 0: Installation

Download/Clone all files from https://github.com/samche42/ShiftScan/tree/main/scripts with ```git clone https://github.com/samche42/ShiftScan.git```

Move into the scripts subfolder in the Shiftscan folder (```cd ShiftScan/scripts/ ```)

Create DSF conda anv using ```conda env create --file=shiftscan.yml```

Activate the environment using ```conda activate shiftscan```

#### Step 1: Raw data format

Currently, the pipeline is designed to take in files with the following raw data format:

| X     | A1: Sample 1 | X     | A2: Sample 2 | X     | A3: Sample 3 | X     | A4: Sample 4 |
|-------|--------------|-------|--------------|-------|--------------|-------|--------------|
| 20.05 | 3.052574772  | 20.05 | 3.229382978  | 20.05 | 3.659788069  | 20.05 | 5.17510602   |
| 20.29 | 3.024562299  | 20.29 | 3.230742839  | 20.29 | 3.63512805   | 20.29 | 5.138402879  |
| 20.68 | 3.024019192  | 20.68 | 3.21544673   | 20.68 | 3.645488649  | 20.68 | 5.126487043  |
| ...   | ...          | ...   | ...          | ...   | ...          | ...   | ...          |
| 84.25 | 4.55696771   | 84.25 | 4.988596048  | 84.25 | 5.40864947   | 84.25 | 5.445997799  |
| 84.49 | 4.552895136  | 84.49 | 4.995465156  | 84.49 | 5.397218846  | 84.49 | 5.445511766  |
| 84.76 | 4.560226904  | 84.76 | 4.991807996  | 84.76 | 5.404517462  | 84.76 | 5.471008942  |

Please ensure your data follows the same format for a seamless analysis.

#### Step 2: Metadata format

Metadata must be stored in the following format (Do NOT include control wells in the metadata):

| ASSAY_PLATE    | SOURCE_PLATE | WELL | COMPOUND    | FRACTION |
|----------------|--------------|------|-------------|----------|
| RBD_040323_01  | Collection1  | P11  | compound_x  | 3        |
| RBD_040323_01  | Collection1  | G18  | compound_y  | 1        |
| RBD_040323_02  | Collection2  | G19  | compound_z  | 0        |
| RBD_040323_02  | Collection2  | G20  | compound_q  | 1        |

 - The 'ASSAY_PLATE' values **MUST** match the names of your raw data files. E.g. One of the input data files was called "RBD_040323_01.txt"
 - The 'SOURCE_PLATE' should be the name of the collection your samples originated from.
 - Do not include control wells in the metadata.
 - If your compounds are pure (I.e. not fractions), leave the 'FRACTION' column as zeros all the way down. 

#### Step 3: Running the analysis

The analysis pipeline can be run with :

```python3 multiprocessor_main.py -i example_input/ -m example_metadata/metadata.txt -o example_output/```

where the parameters are:

- ```-i``` The path to the folder where all raw input files are located
- ```-m``` A tab-delimited file with metadata for all experimental wells. If a compound value is left blank, it is assumed that the well is empty and is not assessed.
- ```-o``` The path to the desired output folder. If it does not exist, a folder will be created with the name specified in this path

Additionally, there are some additional parameters you can provide if you would like to tweak how the data is processed

- ```-c``` A comma-delimited list of which columns contain your controls (Default:  "1,2")
- ```-p``` The number of processors you wish to use (Default: 4)
- ```-x``` The maximum number of control wells allowed to fail per plate. E.g. At default, if 9 control wells fail, the plate is failed. (Default: 8)
- ```-t``` The maximum z-score of control melting temperatures tolerated. I.e. At default, if the z-score of a control well melting temp is 1.6 the well is failed. (Default: 1.5)
- ```-a``` The maximum z-score of control melting curve amplitudes tolerated. I.e. At default, if the z-score of a control well amplitude is 3.1 the well is failed. (Default: 3)
- ```-u``` The maximum relative amplitude of experimental wells allowed (relative to control average) (Default: 6)
- ```-l``` The minimum relative amplitude of experimental wells allowed (relative to control average) (Default: 0.25)
- ```-s``` The desired smoothing factor for smoothing of raw data (default = 0.0005)
- ```-n``` Whether the input data should be normalized per plate. Input options are "y" or "n". (Default: y)

**Note: You can use the ```-h``` flag to list these options in the command line.

The processing of data can take some time depending on the number of plates included in the analysis. If you're using 4 CPUs, processing 100 plates (384-well) takes around 7 minutes. Please see the associated manuscript for additional details on performance. 

Once complete, 4 files would have been generated in the specified output directory:
 - "Final_curves.txt" includes all coordinates for original and cleaned/sliced curves
 - "Final_results.txt" has results from all calculations, including final melting temperatures, amplitudes, failures, reasons for failures etc
 - "Plate_report.txt" is a small table listing which plates passed or failed. For any plate in which 8 or more control wells failed, the entire plate is labelled a failure
 - "Potential_problems.txt" lists any wells that have failed consistently across 3 or more plates. This has no bearing on any of the results and is meant to serve in an informative capacity. i.e. if the same well is failing in several plates, there may be a pipetting issue.

#### Step 4: Visualization

This step is optional. The visualization script should be run using the ```i``` parameter to point to the directory with the four files generated from the previous step. The script is run with:

```python3 visualization.py -i path/to/output/from/previous/step```

e.g.

```python3 visualization.py -i example_output/ ```

This will start up a local Dash server, with a message like so:

>Dash is running on http://127.0.0.1:8050/
>
> * Serving Flask app 'visualization'
> * Debug mode: on

Navigate to the link provided in your web browser
