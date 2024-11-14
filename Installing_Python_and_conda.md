# Installing Python and/or miniconda on Mac

1. Navigate to your terminal (Applications -> Utilities -> Terminal)

2. Type in ```python --version``` and hit enter. If the message ```command not found: python``` pops up, you'll need to install Python. If something like ```Python 3.13.0``` pops up, you can skip ahead to installing miniconda.
___
## Installing Python

1. Navigate to www.python.org/downloads/ in your web browser.

2. Click on the yellow download button to start the download.

3. Once the installer is downloaded (probably in your Downloads folder), double click on the file to start the installer. Continue through the prompts, agree to the terms and then hit the Install button. You can move the installer to the trash once the installation is complete

4. Double click on the Install Certificates.command icon to complete the installation (this should open a terminal windows with a bunch of installation/set-up text). You can close the window once it says that the process has finished successfully.

On to the next step: Installing miniconda!

___

## Installing miniconda
1. Open a new terminal window, type in the following and then hit enter: ```mkdir -p ~/miniconda3```

2. Next, type in the following and then hit enter: ```curl https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh -o ~/miniconda3/miniconda.sh``` This will download the installer.

3. Next, type in the following and then hit enter: ```bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3``` This will run the installer - wait for it to complete before continuing with the next step.

4. Finally, type in the following and then hit enter: ```rm ~/miniconda3/miniconda.sh``` this will remove the installer

5. To fire up conda, type in the following and then hit enter: ```source ~/miniconda3/bin/activate```. You will see the word "(base)" appear on the lefthand side of your active line.

6. If you would like conda to be active whenever you open a terminal window (I prefer it this way, but up to you!), you can run ```conda init --all```. Then close and re-open your terminal.

7. Next, to install ```git```, type in ```conda install -y git``` and then press enter. 

8. You are ready to [install ShiftScan](README.md)!

# Installing Python and/or miniconda on Windows

1. Navigate to your terminal (Start menu -> Programs -> Powershell)
2. Type in ```py --version``` and hit enter. If the message ```Python was not found; run without arguments ...``` pops up, you'll need to install Python. If something like ```Python 3.13.0``` pops up, you can skip ahead to installing miniconda.

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
