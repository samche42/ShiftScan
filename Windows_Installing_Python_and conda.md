# Installing Python and/or miniconda on Windows

1. Navigate to your terminal (Start menu -> Programs -> Windows Powershell) OR just type "Powershell" into the search bar
2. Type in ```python --version``` and hit enter. If the message ```Python was not found; run without arguments ...``` pops up, you'll need to install Python. If something like ```Python 3.13.0``` pops up, you can skip ahead to installing miniconda.

___
## Installing Python

1. Navigate to www.python.org/downloads/ in your web browser.

2. Click on the yellow download button to start the download.

3. Once the installer is downloaded (probably in your Downloads folder), double-click on the file to start the installer. Check the box that says to add Python to your PATH and then continue through the prompts, agree to the terms and then hit the Install button. You can move the installer to the trash once the installation is complete.

4. Open a new PowerShell window and type in ```py --version``` to confirm that Python has been installed successfully

On to the next step: Installing miniconda!

___

## Installing miniconda 

1. Open a new PowerShell terminal window, type in the following and then hit enter: ```curl https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe -o miniconda.exe```

2. Next, type in the following and then hit enter: ```Start-Process -FilePath ".\miniconda.exe" -ArgumentList "/S" -Wait```

3. Next, remove the installer with: ```del miniconda.exe```

4. After installing, you'll want to open the “Anaconda Powershell Prompt” instead of the regular PowerShell terminal

5. Finally, to install ```git```, type in ```conda install -y git``` and then press enter.
 
6. You are ready to [install ShiftScan](README.md)!
