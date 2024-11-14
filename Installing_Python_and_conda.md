# Installing Python and/or miniconda on Mac

1. Navigate to your terminal (Applications -> Utilities -> Terminal)

2. Type in ```python --version``` and hit enter. If the message ```command not found: python``` pops up, you'll need to install Python. If something like ```Python 3.10.12``` pops up, you can skip ahead to installing miniconda.
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

7. Next, to install ```git```, you can should type in ```conda install -y git``` and then press enter. 

8. You are ready to [install ShiftScan](README.md)!
