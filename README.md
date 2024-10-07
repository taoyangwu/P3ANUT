# P3ANUT

## How to run overall project
1. The project was developed using Python 3.12.3. It is recommended to use this version of Python to run the project.
2. Install the required packages within the requirements.txt file
3. Run the unifiedGUI.py file. This is a GUI that contains all the features of the project.

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







# Encountered Errors
- [x] Error: 'module _tkinter not found' 
  - MAC Solution: 'brew install python-tk'