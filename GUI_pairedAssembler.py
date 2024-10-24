import tkinter as tk
from tkinter import ttk
from configPopupV2 import parameterMenu
import csv
from multiprocessedPairAssembler import parse
import json
from pathlib import Path
from runningPopUP import runningPopUp


#---------------------------------------------------------#
#Class: app
#Description: Main class for the GUI, Used for stand alone operation and testing
#---------------------------------------------------------#
class app(tk.Tk):
    def __init__(self, controller):
        super().__init__()
        
        self.geometry("800x600")
        
        frame = pairedAssemblerFrame(self)
        frame.pack(fill="both", expand=True)
        
        self.mainloop()

#---------------------------------------------------------#
#Class: pairedAssemblerFrame
#Description: Frame for the paired assembler GUI
#---------------------------------------------------------#
class pairedAssemblerFrame(tk.Frame):
    def __init__(self, controller):
        super().__init__(controller)
        
        #Variables
        self.forwardFile = ""
        self.reverseFile = ""
        self.outPutFile = ""
        
        self.configFile = ""
        self.parameters = {}
        

        #Grid Configuration
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=3)
        self.grid_columnconfigure(1, weight=1)
        
        #-------#
        self.graphFrame = tk.LabelFrame(self, text="Graph", background="#DDDEE5")
        self.graphFrame.grid(row=0, column=0, sticky="nsew")
        self.graphFrame.grid_columnconfigure(0, weight=1)
        
        #-------#
        
        self.graph = ScrollableFrame(self.graphFrame, csvPath = Path(__file__).parent / "log.csv")
        self.graph.grid(row=0, column=0, sticky="nsew")
        
        #-------#
        
        self.parametersFrame = tk.LabelFrame(self, text="Parameters", background="#C3C3C3")
        self.parametersFrame.grid(row=0, column=1, sticky="nsew")
        self.parametersFrame.grid_columnconfigure(0, weight=1)
        
        #-------#
        
        self.forwardFrame = tk.LabelFrame(self.parametersFrame, text="Forward File", background="#6C8193")
        self.forwardFrame.grid(row=0, column=0, sticky="nsew")
        self.forwardFrame.grid_columnconfigure(0, weight=1)
        
        #-------#
        
        self.forwardFileButton = tk.Button(self.forwardFrame, text="Load File", command=self.loadForwardFile)
        self.forwardFileButton.grid(row=0, column=0, sticky="nsew")
        
        #-------#
        
        self.reverseFrame = tk.LabelFrame(self.parametersFrame, text="Reverse File", background="#7C8594")
        self.reverseFrame.grid(row=1, column=0, sticky="nsew")
        self.reverseFrame.grid_columnconfigure(0, weight=1)
        
        #-------#
        
        self.reverseFileButton = tk.Button(self.reverseFrame, text="Load File", command=self.loadReverseFile)
        self.reverseFileButton.grid(row=0, column=0, sticky="nsew")
        
        #-------#
        
        self.outputFrame = tk.LabelFrame(self.parametersFrame, text="Output File", background="#6C8193")
        self.outputFrame.grid(row=2, column=0, sticky="nsew")
        self.outputFrame.grid_columnconfigure(0, weight=1)
        
        #-------#
        
        self.outputFileButton = tk.Button(self.outputFrame, text="Output File", command=self.output)
        self.outputFileButton.grid(row=0, column=0, sticky="nsew")
        
        #-------#
        
        self.parameterFrame = tk.LabelFrame(self.parametersFrame, text="Settings", background="#7C8594")
        self.parameterFrame.grid(row=3, column=0, sticky="nsew")
        self.parameterFrame.grid_columnconfigure(0, weight=1)
        
        #-------#
        
        self.parameterButton = tk.Button(self.parameterFrame, text="Advanced Settings", command=self.openParameterMenu)
        self.parameterButton.grid(row=0, column=0, sticky="nsew")
        
        #-------#
        
        self.runButton = tk.Button(self.parametersFrame, text="Run", command=self.run)
        self.runButton.grid(row=4, column=0, sticky="nsew")
        
        
        
    def openParameterMenu(self):
        parameterMenu(self, parameters=self.parameters)
        
    #---------------------------------------------------------#
    #Function: loadForwardFile
    #Description: Function to load the forward file, separate from the reverse file due to hardcoded file references
    #Inputs: None
    #Outputs: None
    #---------------------------------------------------------#
    def loadForwardFile(self):
        fileName = tk.filedialog.askopenfilename(initialdir = "/",
                                            title = "Select a File",
                                            filetypes = (("FASTQ Files", "*.fastq"),("TXT File","*.txt")))
        
        if(fileName == ""):
            return
        else:
            self.forwardFile = fileName
            self.forwardFileButton["text"] = self.forwardFile.split("/")[-1]
        
        
    #---------------------------------------------------------#
    #Function: loadReverseFile
    #Description: Function to load the reverse file
    #Inputs: None
    #Outputs: None
    #---------------------------------------------------------# 
    def loadReverseFile(self):
        fileName= tk.filedialog.askopenfilename(initialdir = "/",
                                            title = "Select a File",
                                            filetypes = (("FASTQ Files", "*.fastq"),("TXT File","*.txt")))
        
        if(fileName == ""):
            return
        else:
            self.reverseFile = fileName
            self.reverseFileButton["text"] = self.reverseFile.split("/")[-1]
    
    #---------------------------------------------------------#
    #Function: outputs
    #Description: Function to open a file dialog to select the output file
    #Inputs: None
    #Outputs: None
    #---------------------------------------------------------# 
    def output(self):
        filename = tk.filedialog.asksaveasfilename(filetypes=[("Json Files", "*.json")])
        
        if(filename == ""):
            return
        else:
            self.outPutFile = filename
            self.outputFileButton["text"] = self.outPutFile.split("/")[-1]
        
    #---------------------------------------------------------#
    #Function: run
    #Description: Function to run the paired assembler using the loaded files and parameters
    #Inputs: None
    #Outputs: None
    #---------------------------------------------------------#
    def run(self):
        
        #Check that all the parameters are loaded
        if(self.forwardFile == "" or self.outPutFile == ""):
            return
        
        print("Running")
        csvColNames = ["timeOfRun","forwardFile","ReverseFile","forwardTime","ReverseTime","forwardFileCount","ReverseFileCount","uniqueCount","mergeTime","finalCount","totalTime", "droppedCount"]
        
        outPutData = {}
        
        trimmedParametersPaired = {}
        trimmedParametersGlobal = {}
        
        for key,value in self.parameters.get("pairedAssembler", {}).items():
            trimmedParametersPaired[key] = value["Value"]
            
        for key,value in self.parameters.get("global", {}).items():
            trimmedParametersGlobal[key] = value["Value"]        
        
        runningPopup = runningPopUp(self)

        meta = parse(self.forwardFile, self.reverseFile, outPutData, trimmedParametersGlobal, **trimmedParametersPaired)
        runningPopup.finishedProgram()
        
        
        with open(self.outPutFile, "w") as file:
            json.dump(outPutData, file, indent=4) 
            
        t = []
        for i in csvColNames:
            t.append(meta.get(i, "N/A"))
        
        self.graph.addData(t)
        pass

class ScrollableFrame(ttk.Frame):
    def __init__(self, container, csvPath = ""):
        super().__init__(container)
        canvas = tk.Canvas(self)
        scrollbar = ttk.Scrollbar(self, orient="vertical", command=canvas.yview)
        scrollbar_h = ttk.Scrollbar(self, orient="horizontal", command=canvas.xview)
        self.scrollable_frame = tableEntry(canvas, csvPath)
        

        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )

        canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")

        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.configure(xscrollcommand=scrollbar_h.set)

        scrollbar_h.pack(side="bottom", fill="x")
        scrollbar.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)
        
        
        
    def changeData(self, data):
        self.scrollable_frame.changeData(data)
        
    def addData(self, data):
        self.scrollable_frame.addData(data)
        
class tableEntry(ttk.Frame):
    def __init__(self, container, csvPath, maxLength = 100):
        super().__init__(container)
        
        self.data = []
        self.maxLength = maxLength
        
        self.i = 1
        self.aviableColumns = ["timeOfRun","forwardFile","ReverseFile","forwardTime","ReverseTime","forwardFileCount","ReverseFileCount","uniqueCount","mergeTime","finalCount","totalTime", "droppedCount"]
        self.usedColumns = [0,9, 11, 7, 10,1,2]
        #We want datetime, forward, reverse, and finalcount and time
        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=1)
        self.grid_columnconfigure(3, weight=1)
        self.grid_columnconfigure(4, weight=1)
        
        
        
        self.data = []
        
        #If the csv file is not found -> create it
        if(not Path(csvPath).exists()):
            csvPath = Path(__file__).parent / "log.csv"
            with open(csvPath, "w") as file:
                writer = csv.writer(file)
                writer.writerow(self.aviableColumns)
        
        with open(csvPath, "r") as file:
            reader = csv.reader(file)
            for i in range(maxLength):
                try:
                    t = next(reader)
                    self.data.insert(1, [t[j] for j in self.usedColumns])
                except:
                    break
                
        self.draw()
            
        
    def changeMaxLength(self, maxLength):
        self.maxLength = maxLength
        if(maxLength < len(self.data)):
            self.data = self.data[:maxLength]
        
    def addData(self, newData):
        
        if(len(newData) > 1):
            newData = [newData]
        
        for d in newData:
            self.data.insert(1, [str(d[j]) for j in self.usedColumns])
            
        self.changeMaxLength(self.maxLength)
        self.draw()
                
                
    def draw(self):
        for i in range(len(self.data)):
            for j in range(len(self.data[i])):
                shortendData = self.data[i][j].split("/")[-1]
                ttk.Label(self, text=shortendData).grid(row=i, column=j, sticky="nsew")
            
        
                


def main():
    app(None)

if(__name__ == "__main__"):
    main()
