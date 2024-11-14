# ShiftScan
Analysis and visualization of high-throughput DSF data

> [!IMPORTANT]
> If you do not have Python or miniconda installed, please [follow these instructions](Installing_Python_and_conda.md) to get those set up prior to installing ShiftScan

### Step 1: Installing ShiftScan

1. Navigate to your Downloads folder, or wherever you would like to store the ShiftScan tool

2. Download/Clone all files from the ShiftScan repo with: ```git clone https://github.com/samche42/ShiftScan.git```

**Note:** You can install git using ```conda install -y git```

3. Move into the scripts subfolder in the Shiftscan folder (```cd ShiftScan/scripts/ ```)

4. Create ShiftScan conda environemnet using ```conda env create --file=shiftscan.yml```

5. Activate the environment using ```conda activate shiftscan```

6. You can test that everything works by running: 

**Mac users:** ```python3 multiprocessor_main.py -i example_input/ -m example_metadata/metadata.txt -o example_output/```.

**Windows users:** ```python multiprocessor_main.py -i example_input/ -m example_metadata/metadata.txt -o example_output/```. 

At first, it will seem like nothing is happening, but after about 10 seconds, several messages detailing the pipeline steps will show up.

7. Once the analysis is complete, you can test the visualization tool with:
 
**Mac users:** ```python3 visualization.py -i example_output/```.

**Windows users:** ```python visualization.py -i example_output/```.

Navigate to the address that the Dash is running on (e.g. http://127.0.0.1:8050/) and you can play around and get familiar with the visualizations.
<br/><br/>

#### Installation and usage video
A very brief video of ShiftScan installation (conda already installed) and usage is available on YouTube: https://youtu.be/shN0vUHEGOk

___

### Step 2: Raw data format

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

> [!IMPORTANT]
> Please ensure your data follows the same format for a seamless analysis.

___
### Step 3: Metadata format

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

___
### Step 4a: Running the analysis (default RAM-intensive mode)

The analysis pipeline can be run with :

**Mac users:** ```python3 multiprocessor_main.py -i example_input/ -m example_metadata/metadata.txt -o example_output/```.

**Windows users:** ```python multiprocessor_main.py -i example_input/ -m example_metadata/metadata.txt -o example_output/```.

where the parameters are:

- ```-i``` The path to the folder where all raw input files are located
- ```-m``` A tab-delimited file with metadata for all experimental wells. If a compound value is left blank, it is assumed that the well is empty and is not assessed.
- ```-o``` The path to the desired output folder. If it does not exist, a folder will be created with the name specified in this path

> [!NOTE]
> If you only want the Tm values estimated and no further comparison or analysis, you can add the ```--only_tm``` flag, which will stop the pipeline early and generate a single output file with The estimated Tm value per sigmoidal region detected. This output will **not** work with the default companion tool (See later for details for visualization). Example code to use:
> 
> **Mac users:** ```python3 multiprocessor_main.py -i example_input/ -m example_metadata/metadata.txt -o example_output/ --only_tm```
> **Windows users:** ```python multiprocessor_main.py -i example_input/ -m example_metadata/metadata.txt -o example_output/ --only_tm```

Additionally, there are some additional parameters you can provide if you would like to tweak how the data is processed

- ```-c``` A comma-delimited list of which columns contain your controls (Default:  "1,2")
- ```-p``` The number of processors you wish to use (Default: 4)
- ```-x``` The maximum number of control wells allowed to fail per plate. E.g. At default, if 9 control wells fail, the plate is failed. (Default: 8)
- ```-t``` The maximum z-score of control melting temperatures tolerated. I.e. At default, if the z-score of a control well melting temp is 1.6 the well is failed. (Default: 1.5)
- ```-a``` The maximum z-score of control melting curve amplitudes tolerated. I.e. At default, if the z-score of a control well amplitude is 2.1 the well is failed. (Default: 2)
- ```-u``` The maximum relative amplitude of experimental wells allowed (relative to control average) (Default: 6)
- ```-l``` The minimum relative amplitude of experimental wells allowed (relative to control average) (Default: 0.2)
- ```-s``` The desired smoothing factor for smoothing of raw data (default = 0.0005)
- ```-n``` Whether the input data should be normalized per plate. Input options are "y" or "n". (Default: y)

**Note: You can use the ```-h``` flag to list these options in the command line.

Data processing can take some time depending on the number of plates in the analysis. If you're using 4 CPUs, processing 100 plates (384-well) takes around 7 minutes. Please see the associated manuscript for additional details on performance. 

> [!IMPORTANT]
>Several warnings will potentially be raised during running. These are all handled within the script. Unless they crash the program, these have been handled and are **NOT** a concern. As soon as I figure out out to stop those from being printed out (but still raised for the script to deal with the problem)I'll fix that

Once complete, 4 files would have been generated in the specified output directory:
 - "Final_curves.txt" includes all coordinates for original and cleaned/sliced curves
 - "Final_results.txt" has the results from all calculations, including final melting temperatures, amplitudes, failures, reasons for failures, etc
 - "Plate_report.txt" is a small table listing which plates passed or failed. For any plate in which 8 or more control wells failed, the entire plate is labelled a failure
 - "Potential_problems.txt" lists any wells that have consistently failed across 3 or more plates. This has no bearing on any of the results and is meant to serve in an informative capacity. For example, if the same well is failing in several plates, there may be a pipetting issue.

### Step 4b: Running the analysis (Disk-intensive mode)
If you are limited by available RAM but still need to run hundreds of plates, you may opt to use the disk-intensive mode of ShiftScan. 

The analysis pipeline is largely identical to the RAM-intensive mode in terms of usage, all you have to do is change the name of the script:

```python3 multiprocessor_main_platewise.py -i example_input/ -m example_metadata/metadata.txt -o example_output/```

All available parameters and the option to use the ```--only_tm flag``` are the same as in the RAM-intensive mode. Note that the disk-intensive mode will be slower than the RAM-intensive mode. 

___
### Step 5: Visualization

This step is optional. The visualization script should be run using the ```i``` parameter to point to the directory with the four files generated from the previous step. The script is run with:

```python3 visualization.py -i path/to/output/from/previous/step```

e.g.

**Mac users:** ```python3 visualization.py -i example_output/```.

**Windows users:** ```python visualization.py -i example_output/```.

> [!IMPORTANT]
> If you have opted to use the ```--only_tm``` flag, you cannot use the default visualization tool. Instead, you must run:
>
> **Mac users:** ```python3 visualization_only_tm.py -i example_output/ ```
> **Windows users:** ```python visualization_only_tm.py -i example_output/ ```

In either case, a local Dash server will be fired up, with a message like so:

```
Dash is running on http://0.0.0.0:8050/

 * Serving Flask app 'visualization'
 * Debug mode: off
WARNING: This is a development server. Do not use it in a production deployment. Use a production WSGI server instead.
 * Running on all addresses (0.0.0.0)
 * Running on http://127.0.0.1:8050
 * Running on http://10.162.203.135:8050
Press CTRL+C to quit
```

Copy and paste the address (e.g. ```http://127.0.0.1:8050/```) into your web browser. 
