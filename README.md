# DSF_pipeline
All scripts for analysing and visualizing DSF data

### Step 1: Raw data format
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

### Step 2: Metadata format
Metadata must be stored in the following format (Do not include control wells in the metadata, these are assumed blank):

| ASSAY_PLATE    | SOURCE_PLATE | WELL | COMPOUND    | FRACTION |
|----------------|--------------|------|-------------|----------|
| RBD_040323_01  | Collection1  | P11  | compound_x  | 3        |
| RBD_040323_01  | Collection1  | G18  | compound_y  | 1        |
| RBD_040323_01  | Collection2  | G19  | compound_z  | 0        |
| RBD_040323_01  | Collection2  | G20  | compound_q  | 1        |

Here, the 'ASSAY_PLATE' values **MUST** match the names of your raw data files. 'SOURCE_PLATE' should be the name of the collection your samples originated from. Do not include control wells in the metadata. If your compounds are not the result of fractionation or are pure, leave this column as zeros all the way down. 

### Step 3: Raw data concatenation and processing
To analyze many files at once, the first step is to concatenate the data with the file_concatenator.py script. You will need to navigate to the folder in which your raw data files are stored and then run:

```python3 file_concatenator.py -i input_directory -b Analysis1```

-i is the fuull path to where your raw data files are
-b is a name that you would like your input file to be called

The script will find all files with a "*.txt" extension, transform them to data frames with an additional column 'Origin of data' that will have the file name listed. It is therefore **very important** that your file names match the name listed under 'ASSAY_PLATE' in the metadata file, as this is how the information is linked between the two tables. 
A new file with a "_concatenated.txt" suffix will be created in your current directory, so in this example, our file would be called "Analysis1_concatenated.txt". This will serve as the data input for the next step.

### Step 4: Running the analysis
The analysis pipeline can be run with :

```python3 multiprocessor_main.py -c Analysis1_concatenated.txt -m metadata_file.txt -p num_of_processors```

This can take some time depending on the number of plates included in the analysis. At the time of testing, using 6 processors on a standard MacBook Pro, it took ~1 hour to analyze 100 384-well plates.

Once complete, 4 files would have been generated in the active directory:
 - 
