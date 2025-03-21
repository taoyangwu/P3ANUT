# P3ANUT: Python Pipeline for Phage Analysis through a Normative Unified Toolset

## How to run software
Prerequisites
1.Python 3.12.3 (recommended- The project was developed using Python 3.12.3. It is recommended to use this version to run the project; to ensure compatibility.)

2.Conda (Miniconda or Anaconda) for environment management

2. Install Conda (Miniconda or Anaconda)
Visit the Miniconda documentation or Anaconda.com/download.
(Optional) Submit your email to register with Anaconda.
Download and install the appropriate version for your operating system.

4. Create a Conda Environment
3.1 Open a Terminal
Make sure Conda is available in your PATH. If you are on Windows, use the “Anaconda Prompt” or “Miniconda Prompt.”

3.2 Create a New Environment
Option A: From an environment file (recommended if you have P3ANUT_Env.yml)

'''bash
conda env create -n P3ANUT --file P3ANUT_Env.yml--'''

This command will create an environment named P3ANUT using all required dependencies as defined in the YAML file.

Option B: From scratch

'''bash
conda create -n P3ANUT python=3.12
conda activate P3ANUT
conda install pyyaml numpy pandas matplotlib scikit-learn logomaker anytree conda-forge::python-levenshtein--
This installs the core libraries needed for P3ANUT. Adjust or add additional packages as required.

3.3 Verify Your Environment
'''bash
conda env list--
Make sure the P3ANUT environment is listed, and confirm it’s activated (indicated by (P3ANUT) at the beginning of your terminal prompt).

4. Install Required Packages (Alternative: pip)
If you prefer to install dependencies via pip instead of conda:

'''bash
pip install -r requirements.txt--
Make sure you are inside the P3ANUT environment before running this command.

5. Run the Software
5.1 Unified GUI
To launch the main GUI, run:

'''bash
python unifiedGUI.py--

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

Liam Tucker†1, Ethan Koland†2, Jasmyn Gooding†2, Hassan Boudjelal1, Derek T. Warren1, David Baker3, Maria Marin4, Taoyang Wu2, Chris J. Morris1,5.

† Corresponding authors

1 School of Pharmacy, University of East Anglia, Norwich, Norfolk, NR47TJ, United Kingdom

2 School of Computing Sciences, University of East Anglia, Norwich, Norfolk, NR47TJ, United Kingdom

3 Quadrum Institute of Bioscience, Norwich, Norfolk, NR47UQ, United Kingdom

4 School of Chemistry, University of East Anglia, Norwich, Norfolk, NR47TJ, United Kingdom 

5 School of Pharmacy, University Collage London, Greater London, WC1N 1AX, United Kingdom







# Encountered Errors
- [x] Error: 'module _tkinter not found' 
  - MAC Solution: 'brew install python-tk'
