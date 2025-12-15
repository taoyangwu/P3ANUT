# P3ANUT: Python Pipeline for Peptide Analysis through a Normative Unified Toolset

## Prerequisites
### 0. Prerequisites libaries or softwares
- Conda (Miniconda or Anaconda) for environment management

### 1. Creating Conda Environment - Recommended

#### 1. Install Conda (Miniconda or Anaconda)
   1.1 Visit the [Miniconda documentation](https://www.anaconda.com/docs/getting-started/miniconda/main) or [Anaconda.com/download](https://www.anaconda.com/download).
   
   1.2 Download and install the appropriate version for your operating system.
#### 2. Create a Conda Environment
  ##### 2.1 Open a Terminal
  Make sure Conda is available in your PATH. If you are on Windows, use the “Anaconda Prompt” or “Miniconda Prompt.”
    
  ##### 2.2 Create a New Environment
  ###### 2.2.1 Option A: From an environment file (recommended if you have P3ANUT_Env.yml)
  This command will create an environment named P3ANUT using all required dependencies as defined in the YAML file.

``` bash 
conda env create -n P3ANUT --file P3ANUT_Env.yml
```
  ###### 2.2.2 Option B: From scratch
  This installs the core libraries needed for P3ANUT. Adjust or add additional packages as required.

``` bash
conda create -n P3ANUT python=3.12
conda activate P3ANUT
conda install pyyaml numpy pandas matplotlib scikit-learn logomaker anytree conda-forge::python-levenshtein
```

#### 3. Verify Your Environment
  Make sure the P3ANUT environment is listed, and confirm it’s activated (indicated by (P3ANUT) at the beginning of your terminal prompt).
  ```bash
  conda env list
  ```

#### 4. Actavate the environment
```bash
conda activate P3ANUT
```

### (OPTIONAL Alterinative to Conda) Install Required Packages (Alternative: pip)
#### 1. Install python 3.12 or check the version is installed with one of the following commands. 
```bash
python3.12 --version
python3 --version
```

#### 2. Install the dependencies via pip
```bash
pip install -r requirements.txt
```
#### 3. Make sure you are inside the P3ANUT environment before running this command.

## Running the scripts
### Running the gui

1. Run the Unified GUI
  To launch the main GUI, run:
  
  ```bash
  python unifiedGUI.py
  ```

This graphical interface provides access to all features of P3ANUT.

### How to operate the Paired Assembler GUI
1. Click load a File in the Forward Read section and select a forward FASTQ file. These are denoted by a _001 in the file name.
2. Click load a File in the Reverse Read section and select a reverse FASTQ file. These are denoted by a _002 in the file name.
3. Select an output Files. This is done through the OS file explorer.
4. Click Run Paired Assembler. This will run the paired assembler and output the results to the selected output file.

### How to operate the Sequence Counter GUI
1. Select an Clustering or counting method from the dropdown menu in the Method section.
2. Select an encoding method from the dropdown menu in the Encoding section.
3. Click load a File in the FASTQ File section and select an json or csv file. This file should contain the output of the Paired Assembler.
4. Select an output folder for the results. Multiple files will be outputted to this folder.
5. Select whether to count or cluster the amino acids or the nucleotides sequences.
6. Click Run Sequence Counter. This will run the sequence counter and output the results to the selected output folder.

### How to operate the run unified GUI
1. Click add file and select a csv file to add to the list of files to run. To remove the file from the list, first click the name of the file then click remove file.
2. Add remaining files to the list.
3. Click output file to select the output file for the results.
4. Click run to unify the files. This will output the results to the selected output file.

### How to operate the Volcano Plot
1. Select a csv file to load into the Volcano Plot for File A and File B.
2. Click Create Plot to generate the plot.
3. Adjust the p-value and ratio sliders to adjust the plot.
4. Export the graph or quadrant data to a csv file. This is done through clicking the associated button.

### How to operate the Upset Plot
1. Import all comparision files using the add file button. To remove a file click on the file name and then click remove file.
2. Click run to generate the plot.
3. To export the plot click the export graph button.
4. To export the data click the output file button. In order to select which set associations to export click the associated file names. When they turn green then they are included within the export. The file size should update to reflect the number of files within the file set.

### How to operate the Ranking Plot
1. Import File 1 and File 2 using the select file Buttons
2. Click Create Plot to generate the plot.
3. Adjust the number of point to compare. **You must click create graph button to see updated graph**
4. Graph axis and which file to count using the associate buttons
5. Export the graph or quadrant data to a csv file. This is done through clicking the associated button.

### Loading advanced settings
1. Click the advanced settings button to load the advanced settings. This will bring up a new window.
2. Click Load YAML to load a yaml file. This will load the settings into the advanced settings window.
3. Select the desired script from the dropdown menu on the top of the page.
4. Scroll through the settings and adjust as needed.
5. Click Save YAML to save the settings to a yaml file.
6. Click Save AS to create a new yaml file.
7. Click Apply changes to apply the settings to the script.
8. Click Close to close the advanced settings window without returning the settings to the script.

### Software Development Team
This software application is developed by Ethan Koland, Liam Tucker, Jasmyn Gooding, and Taoyang Wu.

### Reference
Please cite the associated paper written by By:

Uncovering and Correcting Errors in Peptide Phage Display Library Sequencing

by

Liam Tucker, Ethan Koland, Jasmyn Gooding, Hassan Boudjelal, Derek T. Warren, David Baker, Maria Marin, Taoyang Wu, Chris J. Morris.









# Encountered Errors
- [x] Error: 'module _tkinter not found' 
  - MAC Solution: 'brew install python-tk'
