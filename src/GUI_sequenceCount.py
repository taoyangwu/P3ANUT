import tkinter as tk
from tkinter import ttk
from configPopupV2 import parameterMenu
import os
import csv
import json
from sequenceCounter import countCSVFile, countJsonFile
from runningPopUP import runningPopUp


class sequenceCountFrame(tk.Frame):
    def __init__(self, controller):
        super().__init__(controller)
        
        self.forwardFile = ""
        self.reverseFile = ""
        self.outPutFile = ""
        
        self.configFile = ""
        self.parameters = {}
        
        self.countingOptions = ["direct", "OPTICS", "DBSCAN"]
        self.encodingOption = ["BLSOUM", "ONEHOT"]
        self.runOptions = ["AMINO", "DNA", "BOTH"]
        
        self.currentCountingOption = tk.StringVar(value = self.countingOptions[0])
        self.currentEncodingOption = tk.StringVar(value = self.encodingOption[0])
        self.currentRunOption = tk.StringVar(value = self.runOptions[0])
        
        self.currentOutput = tk.StringVar(value = "sequenceCountOutputs")
        
       
        
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=3)
        self.grid_columnconfigure(1, weight=1)
        
        
        
        
        # self.graphFrame = tk.LabelFrame(self, text="Graph", background="blue")
        # self.graphFrame.grid(row=0, column=0, sticky="nsew")
        # self.graphFrame.grid_columnconfigure(0, weight=1)
        
        # self.graph = ScrollableFrame(self.graphFrame, csvPath = "/Users/ethankoland/Desktop/3rd Year Project/code/log.csv")
        # self.graph.grid(row=0, column=0, sticky="nsew")
        
        
        self.parametersFrame = tk.LabelFrame(self, text="Parameters", background="#C3C3C3")
        self.parametersFrame.grid(row=0, column=1, sticky="nsew")
        self.parametersFrame.grid_columnconfigure(0, weight=1)
        
        
        self.methodFrame = tk.LabelFrame(self.parametersFrame, text="Method Selection", background="#6C8193")
        self.methodFrame.grid(row=0, column=0, sticky="nsew")
        self.methodFrame.grid_columnconfigure(0, weight=1)
        
        self.methodOption = tk.OptionMenu(self.methodFrame, self.currentCountingOption, *self.countingOptions)
        self.methodOption.grid(row=0, column=0, sticky="nsew")
        
        self.encodingFrame = tk.LabelFrame(self.parametersFrame, text="Encoding", background="#7C8594")
        self.encodingFrame.grid(row=1, column=0, sticky="nsew")
        self.encodingFrame.grid_columnconfigure(0, weight=1)
        
        self.encodingOption = tk.OptionMenu(self.encodingFrame, self.currentEncodingOption, *self.encodingOption)
        self.encodingOption.grid(row=0, column=0, sticky="nsew")
        
        
        self.dataFrame = tk.LabelFrame(self.parametersFrame, text="Data File", background="#6C8193")
        self.dataFrame.grid(row=2, column=0, sticky="nsew")
        self.dataFrame.grid_columnconfigure(0, weight=1)
        
        self.dataFileButton = tk.Button(self.dataFrame, text="Load File", command=self.loadDataFile)
        self.dataFileButton.grid(row=0, column=0, sticky="nsew")
        
        
        self.outputFrame = tk.LabelFrame(self.parametersFrame, text="Output Folder", background="#7C8594")
        self.outputFrame.grid(row=3, column=0, sticky="nsew")
        self.outputFrame.grid_columnconfigure(0, weight=1)
        
        self.outputFileButton = tk.Button(self.outputFrame, text=self.currentOutput.get(), command=self.output)
        self.outputFileButton.grid(row=0, column=0, sticky="nsew")
        
        self.runFrame = tk.LabelFrame(self.parametersFrame, text="Run Options", background="#6C8193")
        self.runFrame.grid(row=4, column=0, sticky="nsew")
        self.runFrame.grid_columnconfigure(0, weight=1)
        
        self.runOption = tk.OptionMenu(self.runFrame, self.currentRunOption, *self.runOptions)
        self.runOption.grid(row=0, column=0, sticky="nsew")
        
        
        self.parameterFrame = tk.LabelFrame(self.parametersFrame, text="Settings", background="#7C8594")
        self.parameterFrame.grid(row=5, column=0, sticky="nsew")
        self.parameterFrame.grid_columnconfigure(0, weight=1)
        
        self.parameterButton = tk.Button(self.parameterFrame, text="Advanced Settings", command=self.openParameterMenu)
        self.parameterButton.grid(row=0, column=0, sticky="nsew")
        
        self.runButton = tk.Button(self.parametersFrame, text="Run", command=self.run)
        self.runButton.grid(row=6, column=0, sticky="nsew")
        
        self.dataViewer = dataViewer(self)
        self.dataViewer.grid(row=0, column=0, sticky="nsew")
        
        

        
    
        
        
        
        
        
        
    def openParameterMenu(self):
        parameterMenu(self, parameters=self.parameters)
        
    def loadDataFile(self):
        fileName = tk.filedialog.askopenfilename(initialdir = "/",
                                            title = "Select a File",
                                            filetypes = (("CSV Files", "*.csv"),("JSON File","*.json")))
        
        if(fileName == ""):
            return
        else:
            self.forwardFile = fileName
            self.dataFileButton["text"] = self.forwardFile.split("/")[-1]
        
    
        
    def output(self):
        filename = tk.filedialog.askdirectory()
        
        if(filename == ""):
            return
        else:
            self.outPutFile = filename
            self.outputFileButton["text"] = self.outPutFile.split("/")[-1]
            self.currentOutput.set(self.outPutFile)
            
            self.dataViewer.changeFolerLocation(self.outPutFile)
        
    def run(self):
        #load the parameters
        print("Running")
        csvColNames = ["timeOfRun", "forwardFile", "ReverseFile", "forwardTime", "ReverseTime", "forwardFileCount", "ReverseFileCount", "uniqueCount", "multiprocessTime", "finalCount", "totalTime"]
        
        trimmedParametersPaired = {}
        trimmedParametersGlobal = {}
        
        for key,value in self.parameters.get("pairedAssembler", {}).items():
            trimmedParametersPaired[key] = value["Value"]
            
        for key,value in self.parameters.get("global", {}).items():
            trimmedParametersGlobal[key] = value["Value"]
            
        trimmedParametersPaired["countingMethod"] = self.currentCountingOption.get()
        trimmedParametersPaired["encodingMethod"] = self.currentEncodingOption.get()
            
        if(self.forwardFile == "" or self.outPutFile == ""):
            return
        
        runPopup = runningPopUp(self)
        
        if(self.forwardFile.endswith(".csv")):
            countCSVFile(self.forwardFile, self.outPutFile, **trimmedParametersPaired)
        elif(self.forwardFile.endswith(".json")):
            
            runAmino = True if self.currentRunOption.get() in ["AMINO", "BOTH"] else False
            runDNA = True if self.currentRunOption.get() in ["DNA", "BOTH"] else False
            
            countJsonFile(self.forwardFile, runAmino, runDNA, self.outPutFile,  **trimmedParametersPaired)
            
        runPopup.finishedProgram()
        
        self.dataViewer.changeFolerLocation(self.outPutFile)
        
        
        pass
    
class dataViewer(tk.Frame):
    #Have a top frame that will show the logo plots and a bottom that will show the sequences
    def __init__(self, container, csvPath = ""):
        super().__init__(container)
        
        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)
        
        self.grid_columnconfigure(0, weight=1)
        
        self.graphFrame = tk.LabelFrame(self, text="Graph", background="#DDDEE5")
        self.graphFrame.grid(row=0, column=0, sticky="nsew")
        self.graphFrame.grid_columnconfigure(0, weight=1)
        self.graphFrame.grid_rowconfigure(0, weight=1)
        
        self.graph = ScrollableFrame(self.graphFrame, csvPath = csvPath)
        self.graph.grid(row=0, column=0, sticky="nsew")
        
        self.sequenceFrame = tk.LabelFrame(self, text="Sequences", background="#DDDEE5")
        self.sequenceFrame.grid(row=1, column=0, sticky="nsew")
        self.sequenceFrame.grid_columnconfigure(0, weight=1)
        
        self.imgViewer = logoViewer(self.sequenceFrame)
        self.imgViewer.grid(row=0, column=0, sticky="nsew")
        
        # self.sequence = ScrollableFrame(self.sequenceFrame, csvPath = csvPath)
        # self.sequence.grid(row=0, column=0, sticky="nsew")
        
    def changeFolerLocation(self, newFolderLocation):
        self.graph.changeOutputFolder(newFolderLocation)
        self.imgViewer.changeOutputFolder(newFolderLocation)

class ScrollableFrame(ttk.Frame):
    def __init__(self, container, csvPath = "", outputFolder =  ""):
        super().__init__(container)
        
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=2)
        
        self.outputFolder = outputFolder
        self.folderSelections = ["None"]
        self.currentSelection = tk.StringVar(value = csvPath)
        
        self.dropDownFrame = tk.LabelFrame(self, text="Output Selection", background="#7C8594")
        self.dropDownFrame.grid(row=0, column=0, sticky="nsew")
        
       
        self.canvasFrame = ttk.Frame(self)
        self.canvasFrame.grid(row=1, column=0, sticky="nsew")
        
        
        canvas = tk.Canvas(self.canvasFrame, background="#C3C3C3")
        
        scrollbar = ttk.Scrollbar(self.canvasFrame, orient="vertical", command=canvas.yview)
        scrollbar_h = ttk.Scrollbar(self.canvasFrame, orient="horizontal", command=canvas.xview)
        self.tableFrame = tableEntry(canvas, self.currentSelection.get())
        

        self.tableFrame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )
        
        self.dropDown = tk.OptionMenu(self.dropDownFrame, self.currentSelection, *self.folderSelections, command=self.changeData)
        self.dropDown.grid(row=0, column=0, sticky="nsew")
        
        if(outputFolder != ""):
            self.changeOutputFolder(outputFolder)
        

        canvas.create_window((0, 0), window=self.tableFrame, anchor="nw")

        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.configure(xscrollcommand=scrollbar_h.set)

        scrollbar_h.pack(side="bottom", fill="x")
        scrollbar.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)
        
        
        
    def changeData(self, data):
          
        print("Changing Data", data)
        self.tableFrame.changeData(os.path.join(self.outputFolder, data))
        
    def addData(self, data):
        self.tableFrame.addData(data)
        
    def changeOutputFolder(self, folder):
        self.outputFolder = folder
        self.folderSelections = []
        
        #Clear the options from the drop down
        self.dropDown["menu"].delete(0, "end")
        
        for file in os.listdir(self.outputFolder):
            if file.endswith(".csv"):
                self.folderSelections.append(file)
                #Add the file to the drop down
                
                self.dropDown["menu"].add_command(label=file, command=tk._setit(self.currentSelection, file, self.changeData))
        
        
class tableEntry(ttk.Frame):
    def __init__(self, container, csvPath, maxLength = 100):
        super().__init__(container)
        
        self.data = []
        self.maxLength = maxLength
        
        self.i = 1
        self.aviableColumns = ["Sequence", "m_index", "s_index"]
        self.usedColumns = [0,1,2]
        #We want datetime, forward, reverse, and finalcount and time
        self.grid_columnconfigure(0, weight=1)
        self.grid_columnconfigure(1, weight=1)
        self.grid_columnconfigure(2, weight=1)
        
        self.data = []
        
        #Check to see if the csv path exists
        if(csvPath == ""):
            self.data.append(["Sequence", "m_index", "s_index"])
        else:
        
            with open(csvPath, "r") as file:
                reader = csv.reader(file)
                for i in range(maxLength):
                    try:
                        t = next(reader)
                        # self.data.insert(0, [t[j] for j in self.usedColumns])
                        self.data.append([t[j] for j in self.usedColumns])
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
            self.data.append([d[j] for j in self.usedColumns])
            
        self.changeMaxLength(self.maxLength)
        self.draw()
                
                
    def draw(self):
        for i in range(len(self.data)):
            for j in range(len(self.data[i])):
                print(i, j, self.data[i][j])
                ttk.Label(self, text=self.data[i][j]).grid(row=i, column=j, sticky="nsew")
                
    def changeData(self, csvPath):
        print("Changing Data", csvPath)
        #Clear the data
        if(not os.path.exists(csvPath)):
            print("File does not exist")
            return
        
        self.data = []
        with open(csvPath, "r") as file:
            reader = csv.reader(file)
            for i in range(self.maxLength):
                try:
                    t = next(reader)
                    print(t)
                    self.data.append([t[j] for j in self.usedColumns])
                except:
                    break
                
        self.draw()
        
        
        
            
            
class logoViewer(tk.Frame):
    def __init__(self, container, logoDirectory = ""):
        super().__init__(container)
        
        self.grid_columnconfigure(0, weight=1)
        
        #Have a drop down menu to select the logo
        #How to pass into the directroy
        self.avaiableLogos = [None]
        self.currentLogo = tk.StringVar(value = "")
        self.logoDirectory = logoDirectory
        
        self.dropDownFrame = tk.LabelFrame(self, text="Logo Selection", background="#6C8193")
        self.dropDownFrame.grid(row=0, column=0, sticky="nsew")
        
        self.dropDown = tk.OptionMenu(self.dropDownFrame, self.currentLogo, *self.avaiableLogos)
        self.dropDown.grid(row=0, column=0, sticky="nsew")
        
        self.imageLable = tk.Label(self, background="#DDDEE5", width=self.winfo_width())
        self.imageLable.grid(row=1, column=0, sticky="ew")
        
        self.loadLogoSelection()
        # self.loadLogo("/Users/ethankoland/Desktop/3rd Year Project/code/sequenceCountOutputs/AminoLogo_10.png")
        
    def loadLogoSelection(self):
        
        #Clear the drop down
        self.dropDown["menu"].delete(0, "end")
        
        if(self.logoDirectory == ""):
            return
        
        #Loop through the files in the logo directory and add them to the avaiableLogos
        for file in os.listdir(self.logoDirectory):
            if file.endswith(".png"):
                self.avaiableLogos.append(file)
                self.dropDown["menu"].add_command(label=file, command=tk._setit(self.currentLogo, file, self.loadLogo))
                
                
                
    def loadLogo(self, img):
        print("Loading Logo", img)
        #Clear the canvas
        # self.imageCanvas.delete("all")
        #Load the image
        t = os.path.join(self.logoDirectory, img)
        self.imgObj = tk.PhotoImage(file=os.path.join(self.logoDirectory, img), width=self.winfo_width())
        self.imageLable.configure(image=self.imgObj)
        
    def changeOutputFolder(self, folder):
        self.logoDirectory = folder
        self.loadLogoSelection()
        
        
    
                


def main():
    root = tk.Tk()
    
    
    root.geometry("800x600")
    
    mainFrame = sequenceCountFrame(root)
    mainFrame.pack()
    
    root.mainloop()
    

if(__name__ == "__main__"):
    main()
