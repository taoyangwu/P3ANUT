import tkinter as tk
from tkinter import ttk
from configPopupV2 import parameterMenu
import os
import pandas as pd


import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle
import numpy as np

from runningPopUP import runningPopUp

#---------------------------------------------------------#
#Class: app
#Description: Sample class to create a tkinter application
#---------------------------------------------------------#
class app(tk.Tk):
    def __init__(self, controller):
        super().__init__()
        
        self.geometry("800x600")
        
        up = upsetPlotFrame(self)
        up.pack(fill="both", expand=True)
        
        self.mainloop()
    

#---------------------------------------------------------#
#Class: upsetPlotFrame
#Description: Class that creates and maintains the UI for the upset plot
#---------------------------------------------------------#
class upsetPlotFrame(tk.Frame):
    def __init__(self, controller):
        super().__init__(controller)
        
        self.fileNamesDict = {}
        self.key = 0
    
        
        self.outPutFile = ""
        
        #Array of tk Boolean vars
        self.fileSelectionVars = []
        
        
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=3)
        self.grid_columnconfigure(1, weight=1)
    
        
        self.parametersFrame = tk.LabelFrame(self, text="Parameters", background="#C3C3C3")
        self.parametersFrame.grid(row=0, column=1, sticky="nsew")
        self.parametersFrame.grid_columnconfigure(0, weight=1)
        
        self.filesFrame = tk.LabelFrame(self.parametersFrame, text="Files", background="#6C8193")
        self.filesFrame.grid(row=0, column=0, sticky="nsew")
        self.filesFrame.grid_columnconfigure(0, weight=1)
        self.filesFrame.grid_columnconfigure(1, weight=1)
        
        self.fileListBox = tk.Listbox(self.filesFrame)
        self.fileListBox.grid(row=1, column=0, columnspan=2, sticky="nsew")
        
        self.fileListBox.bind("<<ListboxSelect>>", self.listBoxSelect)
        
        self.addFileButton = tk.Button(self.filesFrame, text="Add File", command=self.addFile)
        self.addFileButton.grid(row=0, column=0, sticky="nsew")
        
        self.removeFileButton = tk.Button(self.filesFrame, text="Remove File", command=self.removeFile)
        self.removeFileButton.grid(row=0, column=1, sticky="nsew")
        
        self.outputFrame = tk.LabelFrame(self.parametersFrame, text="Output", background="#7C8193")
        self.outputFrame.grid(row=1, column=0, sticky="nsew")
        self.outputFrame.grid_columnconfigure(0, weight=1)
        
        self.outputCheckBoxes = tk.Listbox(self.outputFrame, background="#15253F", selectmode="browse")
        self.outputCheckBoxes.grid(row=0, column=0, sticky="nsew")
        
        self.outputText = tk.Label(self.outputFrame, text="File Count : 0")
        self.outputText.grid(row=1, column=0, sticky="nsew")

        
        self.outputFileButton = tk.Button(self.outputFrame, text="Output File", command=self.output)
        self.outputFileButton.grid(row=2, column=0, sticky="nsew")
        
        self.outputCheckBoxes.bind("<<ListboxSelect>>", self.checkboxesSelect)
        
        
        self.runButton = tk.Button(self.parametersFrame, text="Run", command=self.run)
        self.runButton.grid(row=5, column=0, sticky="nsew")
        
        self.exportGraphButton = tk.Button(self.parametersFrame, text="Export Graph", command=self.exportGraph)
        self.exportGraphButton.grid(row=6, column=0, sticky="nsew")
        
    def listBoxSelect(self, event):
        #Used for debugging
        print(self.fileListBox.curselection())
        
    def exportGraph(self):
        fileName = tk.filedialog.asksaveasfilename(filetypes=[("PNG Files", "*.png")])
        
        #If the user cancels the file selection, return
        if(fileName == ""):
            return
        
        self.upsetPlot.exportGraph(fileName)
        
    def checkboxesSelect(self, event):
        
        if(len(self.outputCheckBoxes.curselection()) == 0):
            return
        
        index = self.outputCheckBoxes.curselection()[0]
        
        self.fileSelectionVars[index].set(not self.fileSelectionVars[index].get())
              
        self.outputCheckBoxes.itemconfig(index, bg="green" if self.fileSelectionVars[index].get() else "gray")
        
        #Update the textbox to display the insertion count
        self.key = 0
        for i, value in enumerate(self.fileSelectionVars[::-1]):
            if(value.get()):
                self.key += 2 ** i
                
        self.outputText["text"] = f"File Count : {self.upsetPlot.getInsertionCounts(self.key)}"
        
        
    
        
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
        filename = tk.filedialog.asksaveasfilename(filetypes=[("Csv Files", "*.csv")])
        
        if(filename == ""):
            return
                
        upsetPlot.fileOutput([self.fileNamesDict[i] for i in self.fileNamesDict], filename, self.key)

            
    def addFile(self):
        fileName = tk.filedialog.askopenfilename(initialdir = "/",
                                            title = "Select a File",
                                            filetypes = [("CSV Files", "*.csv")])
        
        
        
        if(fileName == ""):
            return
        
        #File paths are long, so shorten the name increases the readability
        shortName = f"{fileName.split('/')[-1]} - {os.path.getsize(fileName)}" 
        print(shortName)
        if(shortName not in self.fileNamesDict):
            self.fileNamesDict[shortName] = fileName
            self.fileListBox.insert(tk.END, shortName)
            
            self.fileSelectionVars.append(tk.BooleanVar(value=False))
            
            i = len(self.fileSelectionVars) - 1
            
            self.outputCheckBoxes.insert(tk.END, chr(65 + i))
            self.outputCheckBoxes.itemconfig(i, bg="gray")
        
        
            
    def removeFile(self):
        if(len(self.fileListBox.curselection()) == 0):
            return
        
        index = self.fileListBox.curselection()[0]
        
        #Get the value at the index
        fileName = self.fileListBox.get(index)
        self.fileListBox.delete(index)
        self.fileNamesDict.pop(fileName)
        
        print(self.fileNamesDict)
        
        i = len(self.fileSelectionVars) - 1
        
        self.fileSelectionVars.pop(i)
        
        self.outputCheckBoxes.delete(i)
        
        
    def run(self):
        fileNames = [self.fileNamesDict[i] for i in self.fileNamesDict]
        
        if(len(fileNames) <= 1 ):
            return
        
        runPopUp = runningPopUp(self)
        self.upsetPlot = upsetPlot(self, fileNames)
        self.upsetPlot.grid(row=0, column=0, sticky="nsew")
        runPopUp.finishedProgram()
    

        
#---------------------------------------------------------#
#Class: upsetPlot
#Description: Class that creates the upset plot
#---------------------------------------------------------#
class upsetPlot(tk.Frame):
    def __init__(self, controller, fileNames = []) -> None:
        super().__init__(controller)
        
        self.insertsectionDict = {}
        
        self.fig = Figure(figsize=(14, 6), dpi=100)
        # Above works, but causes issue of figure having a lot of white space. Still needs work
        self.figureCanvas = FigureCanvasTkAgg(self.fig, self)
        
        #Set the title for the figure
        self.fig.suptitle("Upset Plot")
        
        plt.style.use('fivethirtyeight')
        plt.tight_layout()
        
        #Creating subplots for the upset plot
        self.ax = self.fig.subplots(2,2, gridspec_kw={'height_ratios': [4, 1], 'width_ratios': [1, 4]})
        
        self.ax1 = self.ax[0,0]
        self.ax2 = self.ax[0,1]
        self.ax3 = self.ax[1,0]
        self.ax4 = self.ax[1,1]
        
        if(len(fileNames) > 0):
            self.chanceFileNames(fileNames)
            
        # self.fig.tight_layout()
        # self.fig.subplots_adjust(wspace=0, hspace=0)
        
        self.figureCanvas.get_tk_widget().pack()
        
    def getInsertionCounts(self, key):
        return self.insertsectionDict.get(key, "N/A")
    
    def exportGraph(self, fileName):
        self.fig.savefig(fileName)
        
    
    def createUpSetGraph(self):
        pass
    
    def barSubGraph(self):
        # ax.bar(tempAxis, data, align='center')
        
        # for i in range(len(data)):
        #     ax.bar(i, data[i],0.5, align='center', color="black")
        
        maxInsersection = self.totalUnionSize
       
        #Create the labels for the graph 
        t = [f"{round(x/maxInsersection * 100)}%" for x in self.sortedInsersectionCounts]
        
        for i in range(len(self.sortedInsersectionCounts)):
            rect = Rectangle((i,0), 0.8,self.sortedInsersectionCounts[i] , color='black', alpha=0.5)
            # modifedSizeLabel = np.format_float_scientific(self.sortedInsersectionCounts[i], precision = 2)
            label = self.ax2.text(i + 0.4, self.sortedInsersectionCounts[i], t[i] , ha="center", va="bottom", color="black", fontsize="x-small")
            
            self.ax2.add_patch(rect)
        self.ax2.set_xlim(0,len(self.sortedInsersectionCounts))

        self.ax2.set_ylim(0,max(self.sortedInsersectionCounts) * 1.2)
        self.ax2.set(ylabel='Number of intersections')
        self.ax2.yaxis.set_label_position("right")
        self.ax2.yaxis.tick_right()
        self.ax2.yaxis.label.set_size(10)
        self.ax2.tick_params(axis='y', labelsize=7)
        
        
        self.ax2.yaxis.major.formatter._useMathText = True
        self.ax2.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.1e}'))
        
        

        self.ax2.xaxis.set_visible(False)
    
    
    def setAssioationGraph(self):
        # self.ax4.get_xaxis().set_visible(False)
        self.ax4.set_xticklabels([])
        self.ax4.set_xticks([])
        
        self.ax4.set_xlim(0, len(self.sortedAxis))

        positives = []
        negatives = []
        connectingLines = []
        
        for i, key in enumerate(self.sortedAxis):
            t = np.arange(0, len(self.fileNames))
            t = np.power(2, t)
            
            #Flip the bits so that the array represents the sets
            t = np.flip(t)
            
            #Match the key to the sets
            t = np.bitwise_and(key, t)
            t = np.sign(t)
            
            connectingLines.append([])
            
            for j, value in enumerate(t):
                if(value == 1):
                    positives.append([i + 0.4,j])
                    connectingLines[i].append(j)
                else:
                    negatives.append([i + 0.4,j])
                    
        for i, line in enumerate(connectingLines):
            if(len(line) > 1):
                self.ax4.plot([i + 0.4, i + 0.4], [line[0], line[-1]], c="green")
                
        positives = np.array(positives)
        negatives = np.array(negatives)           
        

        
        background = np.array([x%2 ==  0 for x in range(len(self.fileNames))])
        background = background * len(self.sortedAxis)
        self.ax4.barh(self.fileNames, background, align='center', color="black", alpha=0.1)
        self.ax4.xaxis.label.set_size(10)
        self.ax4.get_yaxis().set_visible(False)
        self.ax4.set(xlabel='Set Association')
        
        self.ax4.scatter(positives[:,0], positives[:,1], c="green", s = 50)
        self.ax4.scatter(negatives[:,0], negatives[:,1], c="black", alpha = 0.5, s = 30, zorder=2)
        self.ax4.scatter(negatives[:,0], negatives[:,1], c="white", alpha = 1, s = 20, zorder=3)
    
    def fileSizeGraph(self):
        reference = "ABCDEFGHIJKLMNOP"
        
        shortenedFileNames = [reference[i] for i in range(len(self.fileNames))]
        
        self.ax3.barh(shortenedFileNames,self.fileLengths, color="black")
        self.ax3.invert_xaxis()
        self.ax3.yaxis.tick_right()
        self.ax3.yaxis.set_label_position("right")
        self.ax3.set(xlabel='File size')
        self.ax3.xaxis.label.set_size(10)
        self.ax3.tick_params(axis='x', labelsize=8)
        
        
        
        pass
    
    def namesGraph(self):
        #Display the names allong with the shortended names of the file names
        
        reference = "ABCDEFGHIJKLMNOP"
        shortenedFileNames = [reference[i] for i in range(len(self.fileNames))]

        self.ax1.set_xlim(0, len(self.fileNames)+1)
        
        trimmedFileNames = [os.path.basename(i) for i in self.fileNames]
        
        for i in range(len(self.fileNames)):
            self.ax1.text(i + 1, 0.5, f"{shortenedFileNames[i]} - {trimmedFileNames[i]}", ha="center", va="center", fontsize=8, rotation=90)
            
        self.ax1.set_yticklabels([])
        self.ax1.set_yticks([])
        
        self.ax1.set_xticklabels([])
        self.ax1.set_xticks([])
        
        self.ax1.set(ylabel='File Names')
        self.ax1.yaxis.label.set_size(10)
        
    
    @staticmethod
    def fileInsersection(key, fileSets):
        #Create a binary number where each bit represents if a set is in or not
        refrenceNumber = np.arange(0, len(fileSets))
        refrenceNumber = np.power(2, refrenceNumber)
        
        #Flip the bits so that the array represents the sets
        refrenceNumber = np.flip(refrenceNumber)
        
        #Match the key to the sets
        refrenceNumber = np.bitwise_and(key, refrenceNumber)
        refrenceNumber = np.sign(refrenceNumber)
        
        #Seperate the files to be unioned and the files to be differenced
        unionOrDifferencedToggle = np.argsort(refrenceNumber * -1, kind="stable")

        intersectionSet = pd.Index(fileSets[unionOrDifferencedToggle[0]])
        difference = pd.Index([])
        
   
        if(max(refrenceNumber) == 0):
            return []

        for i in range(1, len(unionOrDifferencedToggle)):

            comparision = pd.Index(fileSets[unionOrDifferencedToggle[i]])
            
            t = refrenceNumber[unionOrDifferencedToggle[i]]
 

            if(refrenceNumber[unionOrDifferencedToggle[i]] == 1):
                
                intersectionSet = intersectionSet.intersection(comparision)
            else:
                # s = np.setdiff1d(s, c)
                difference = difference.union(comparision)
                
        # print(len(unions), len(difference))
        return intersectionSet.difference(difference)
    
    @staticmethod
    def fileSetComparision(fileNames, intersectionDict = {}):
        maxNumber = 2 ** len(fileNames)
        insersectionConuts = []
        fileLengths = []
        
        totalFileUnion = pd.Index([])
        files = []
        for i in fileNames:
            pdData = pd.read_csv(i)
            fileLengths.append(len(pdData))
            files.append(sequenceValues := pdData.loc[:, "sequence"].values)
            totalFileUnion = totalFileUnion.union(sequenceValues)
            
        totalUnionSize = len(totalFileUnion)
        
        for i in range(maxNumber):
            
            insersectionConuts.append(t := len(upsetPlot.fileInsersection(i, files)))
            intersectionDict[i] = t
            
        return insersectionConuts, fileLengths, totalUnionSize
    
    @staticmethod
    def fileOutput(fileNames, output,key):
        
        #Create a binary number where each bit represents if a set is in or not
        refrenceNumber = np.arange(0, len(fileNames))
        refrenceNumber = np.power(2, refrenceNumber)
        
        #Flip the bits so that the array represents the sets
        refrenceNumber = np.flip(refrenceNumber)
        
        #Match the key to the sets
        refrenceNumber = np.bitwise_and(key, refrenceNumber)
        refrenceNumber = np.sign(refrenceNumber)
        
        files = []
        dfs = []
        for i, fileName in enumerate(fileNames):
            pdData = pd.read_csv(fileName)
            files.append(pdData.loc[:, "sequence"].values)
            if(refrenceNumber[i] == 1):
                pdData.set_index("sequence", inplace=True)
                dfs.append(pdData)
            
            
        t = upsetPlot.fileInsersection(key, files)
        
        #get the data from the files
        counts = None
        for df in dfs:
            if(counts is None):
                counts = df.loc[t]["mean"].values
            else:
                counts = np.add(counts, df.loc[t]["mean"].values)
            
        
        nedDF = pd.DataFrame({"sequence" : t, "mean" : counts})
        nedDF.set_index("sequence", inplace=True)
        nedDF.sort_values(by=["mean"], inplace=True, ascending=False)
        nedDF.insert(1, "std", 0)
        nedDF.to_csv(output)
        
        
    
    
        
    
    def chanceFileNames(self, fileNames):
        
        self.fileNames = fileNames
        
        self.insertsectionDict = {}
        
        self.insersectionConuts, self.fileLengths, self.totalUnionSize = upsetPlot.fileSetComparision(fileNames, self.insertsectionDict)
        
        self.sortedAxis = np.argsort([t * -1 for t in self.insersectionConuts], kind="stable")
        self.sortedInsersectionCounts = np.take(self.insersectionConuts, self.sortedAxis)
        
        self.barSubGraph()
        self.setAssioationGraph()
        self.fileSizeGraph()
        self.namesGraph()
        
        pass

        
    
    
if __name__ == "__main__":
    app(None)
    #print(upsetPlot.fileOutput(["testRunData/RunUnified/GAL_LB-Uni.csv", "testRunData/RunUnified/GAL-BSA-Uni.csv", "Amino_9.csv"], "test.csv", 7))
