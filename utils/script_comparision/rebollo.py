
import matlab.engine
import os
import shutil
import numpy as np

def rebollo(f):
    eng = matlab.engine.start_matlab()
    eng.cd("RebolloScripts", nargout=0)
    dataPath =  file
    filedir, filename = os.path.split(dataPath)

    #Check if there is a existing file, delete it
    outputFolderReplacement = os.path.join(filedir, filename.replace(".fastq", "_BC"))
    shutil.rmtree(outputFolderReplacement, ignore_errors=True)


    result= eng.Step1('inname', filename, 'indir', "../" + filedir, nargout=0)

    #Change the directory to the output directory
    outputFolder = os.path.join(filedir, filename.replace(".fastq", "_BC"))


    files = [os.path.join(outputFolder, x) for x in os.listdir(outputFolder)]
    fileSizes = [os.path.getsize(x) for x in files]
    idexOfMax = files[np.argmax(fileSizes)]

    print(idexOfMax)
    filedir, filename = os.path.split(idexOfMax)
    result = eng.Step2('inname', filename, 'indir', "../" + filedir, 'start','TCTTGT','end','TTCGAT',
                        'uplimit',16,'DOWNlimit',8,'fixerr',1000,'badmax',2,'ACGTOnly',False,nargout=0)

    #Find the amount of equences in the file
    outputFolder = os.path.join(filedir, f"Translation_{filename.replace('.txt', '')}")
    outputFileName = f"fixerrTranslated_{filename.replace('.txt', '')}_GOOD.txt"

    print(outputFolder, outputFileName)
    
    shutil.rmtree(outputFolderReplacement, ignore_errors=True)