# Installing Python and/or miniconda on Mac

1. Navigate to your terminal (Applications -> Utilities -> Terminal)

2. Type in ```python3 --version``` and hit enter. If the message ```command not found: python3``` pops up, you'll need to install Python. If something like ```Python 3.13.0``` pops up, you can skip ahead to installing miniconda.
___

## Installing Python

1. Navigate to www.python.org/downloads/ in your web browser.

2. Click on the yellow download button to start the download.

3. Once the installer is downloaded (probably in your Downloads folder), double-click on the file to start the installer. Continue through the prompts, agree to the terms, and then hit the Install button. You can move the installer to the trash once the installation is complete

4. Double click on the Install Certificates.command icon to complete the installation (this should open a terminal window with a bunch of installation/set-up text). You can close the window once it says the process has finished successfully.

On to the next step: Installing miniconda!

___

## Installing miniconda
1. Open a new terminal window

2. Next, we're going to go to the [Miniconda website](https://docs.anaconda.com/miniconda/install/) to find the correct version of miniconda for your computer. Click on the drop-down for MacOS and pick the tab that is right for you (I.e. M1/2 chip, or Intel).

3. Copy the correct terminal command into the terminal and hit enter. This will download the installer.

4. Next, type in the following and then hit enter: ```bash miniconda3/miniconda.sh``` This will run the installer - wait for it to complete before continuing with the next step.
  
5. To fire up conda, type in the following and then hit enter: ```source ~/miniconda3/bin/activate```. You will see the word "(base)" appear on the left hand side of your active line.

8. Next, to install ```git```, type in ```conda install -y git``` and then press enter. 

9. You are ready to install ShiftScan!

___

## Installing ShiftScan

1. Navigate to your Downloads folder, or wherever you would like to store the ShiftScan tool

2. Download/Clone all files from the ShiftScan repo with: 

```git clone https://github.com/samche42/ShiftScan.git```

> [!NOTE]
> You can install git using ```conda install -y git```

3. Move into the scripts subfolder in the Shiftscan folder (```cd ShiftScan/scripts/ ```)

4. Create ShiftScan conda environment using ```conda env create --file=shiftscan.yml```

5. Activate the environment using ```conda activate shiftscan```

6. You can test that everything works by running: 

```python3 shiftscan.py -i example_input/ -m example_metadata/metadata.txt -o example_output/```.

At first, it will seem like nothing is happening, but after about 10 seconds, several messages detailing the pipeline steps will show up.

7. Once the analysis is complete, you can test the visualization tool with:
 
```python3 shiftscan_viewer.py -i example_output/```.

8. Navigate to the address that the Dash is running on (e.g. http://127.0.0.1:8050/) and you can play around and get familiar with the visualizations.

9. When you're done, close the window and then shut down the Dash server with Ctrl+C
<br/><br/>

___

## Using ShiftScan

1. Each time you run ShiftScan, you'll need to ensure that the ShiftScan conda environment is activated. To do this, you can run
```conda activate shiftscan```

2. Navigate to the scripts folder in the ShiftScan folder, e.g.:
```cd Downloads/ShiftScan/scripts```

3. Run ShiftScan using the locations of your input files as input, e.g.:

```python3 shiftscan.py -i complete/path/to/your/own/data -m complete/path/to/your/metdata/file/metdata.txt -o complete/path/to/wherever/you/want/to/store/output```

4. Similarly, when it's finished running, you run the visualization (from the ShiftScan/scripts folder) by running:
```python3 shiftscan_viewer.py -i /path/to/wherever/you/put/the/output```
