# Installing Python and/or miniconda on Mac

1. Navigate to your terminal (Applications -> Utilities -> Terminal)

2. Type in ```python3 --version``` and hit enter. If the message ```command not found: python3``` pops up, you'll need to install Python. If something like ```Python 3.13.0``` pops up, you can skip ahead to installing miniconda.
___
## Installing Python

1. Navigate to www.python.org/downloads/ in your web browser.

2. Click on the yellow download button to start the download.

3. Once the installer is downloaded (probably in your Downloads folder), double click on the file to start the installer. Continue through the prompts, agree to the terms and then hit the Install button. You can move the installer to the trash once the installation is complete

4. Double click on the Install Certificates.command icon to complete the installation (this should open a terminal windows with a bunch of installation/set-up text). You can close the window once it says that the process has finished successfully.

On to the next step: Installing miniconda!

___

## Installing miniconda
1. Open a new terminal window

2. Next, we're going to go to the [Miniconda website](https://docs.anaconda.com/miniconda/install/) to find the correct version of miniconda for your computer. Click on the drop down for MacOS and pick the tab that is right for you (I.e. M1/2 chip, or Intel).

3. Copy the correct terminal command into the terminal and hit enter. This will download the installer.

4. Next, type in the following and then hit enter: ```bash miniconda3/miniconda.sh``` This will run the installer - wait for it to complete before continuing with the next step.
  
5. To fire up conda, type in the following and then hit enter: ```source ~/miniconda3/bin/activate```. You will see the word "(base)" appear on the lefthand side of your active line.

8. Next, to install ```git```, type in ```conda install -y git``` and then press enter. 

9. You are ready to [install ShiftScan](README.md)!
