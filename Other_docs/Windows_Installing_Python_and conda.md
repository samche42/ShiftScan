# Installing Python and/or miniconda on Windows

1. Navigate to your terminal (Start menu -> Programs -> Windows Powershell) OR just type "Powershell" into the search bar

2. Open Powershell by clicking on the program 

3. In the newly opened Powershell terminal, type in ```python --version``` and hit enter.

4. If the message ```Python was not found; run without arguments ...``` pops up, you'll need to install Python. If something like ```Python 3.13.0``` pops up, you can skip ahead to installing miniconda.

___
## Installing Python

1. Navigate to www.python.org/downloads/ in your web browser.

2. Click on the yellow download button to start the download.

3. Once the installer is downloaded (probably in your Downloads folder), double-click on the file to start the installer.

   a. Check the box that says to add Python to your PATH

   b. Continue through the prompts

   c. Agree to the terms

   d. Click the Install button.
   
   e. You can move the installer to the trash once the installation is complete.

5. Open a new PowerShell window and type in ```python --version``` to confirm that Python has been installed successfully

On to the next step: Installing miniconda!

___

## Installing miniconda 

1. Open a new PowerShell terminal window, type in the following and then hit enter: ```curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe -o miniconda.exe```. This downloads the miniconda installer.

2. Next, we need to actually run the installer. Type in the following and then hit enter: ```Start-Process -FilePath ".\miniconda.exe" -ArgumentList "/S" -Wait```

3. When the installation is complete, you can remove the installer with: ```del miniconda.exe```

4. Close the PowerShell terminal

5. Next, open the **“Anaconda Powershell Prompt”** *instead of* the regular PowerShell terminal

6. Finally, to install ```git```, type in ```conda install -y git``` and then press enter.

___

## Installing ShiftScan

1. In the Anaconda Powershell Prompt terminal window navigate to your Downloads folder (or wherever you would like to store the ShiftScan tool). You can do this by typing in: ```cd .\Downloads\```

2. Next, let's use git to clone the ShiftScan repo: ```git clone https://github.com/samche42/ShiftScan.git```

3. Great! Next step, navigate into the ShiftScan folder and then into the scripts subfolder with: ```cd \ShiftScan\scripts```

4. Now, we'll need to set up the conda environment for ShiftScan: ```conda env create --file=shiftscan.yml```

5. Activate the environment using ```conda activate shiftscan```

6. Test that the script is working with ```python multiprocessor_main.py -i example_input -m example_metadata/metadata.txt -o example_output```

> [!IMPORTANT]
> The first time you run ShiftScan it may take a long time and appear to hang. It is running and everything is fine, promise! The first time it's run on a new system, it creates a ```__pycache__``` folder, where it stores all it's important imports and ensures rapid runs for every subsequent usage. 

___

## Running ShiftScan

Any time that you want to run ShiftScan, you'll have to:
1. Activate the conda environment with ```conda activate shiftscan```.

2. Navigate into the ShiftScan/scripts folder (unless you added this folder to your PATH)

3. Run the program pointing to your data: e.g. ```python multiprocessor_main.py -i /path/to/my/input/data -m path/to/my/metdata.txt -o path/to/wherever/you/want/to/store/output```

4. To run the visualization tool: ```python visualization.py -i path/to/wherever/you/stored/the/output```

> [!NOTE]
> If the provided address (http://0.0.0.0:8050) does not work. Please try entering ```http://localhost:8050``` instead.
