
import pandas as pd
import numpy as np

import tkinter as tk
from tkinter import ttk
from configPopupV2 import parameterMenu
import os
import csv
import json

class app(tk.Tk):
    def __init__(self, controller):
        super().__init__()
        
        self.geometry("800x600")
        
        frame = runUnifierFrame(self)
        frame.pack(fill="both", expand=True)
        
        self.mainloop()
        
class runUnifierFrame(tk.Frame):
    def __init__(self, controller):
        super().__init__(controller)
        
        self.fileNamesDict = {}
    
        
        self.outPutFile = ""
        
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=3)
        self.grid_columnconfigure(1, weight=1)
    
        
        self.parametersFrame = tk.LabelFrame(self, text="Parameters", background="#DDDEE5")
        self.parametersFrame.grid(row=0, column=1, sticky="nsew")
        self.parametersFrame.grid_columnconfigure(0, weight=1)
        
        self.filesFrame = tk.LabelFrame(self.parametersFrame, text="Files", background="#C3C3C3")
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
        
        self.outputFileButton = tk.Button(self.outputFrame, text="Output File", command=self.output)
        self.outputFileButton.grid(row=0, column=0, sticky="nsew")
        
        self.tableView = ScrollableFrame(self, csvPath="")
        self.tableView.grid(row=0, column=0, sticky="nsew")
        
        
        self.runButton = tk.Button(self.parametersFrame, text="Run", command=self.run)
        self.runButton.grid(row=5, column=0, sticky="nsew")
        
        
    def listBoxSelect(self, event):
        
        
        print(self.fileListBox.curselection())   
        
        
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
        filename = tk.filedialog.asksaveasfilename(filetypes=[("Csv Files", "*.csv")])
        
        if(filename == ""):
            return
        else:
            self.outPutFile = filename
            self.outputFileButton["text"] = self.outPutFile.split("/")[-1]

            
    def addFile(self):
        fileName = tk.filedialog.askopenfilename(initialdir = "/",
                                            title = "Select a File",
                                            filetypes = [("CSV Files", "*.csv")])
        
        
        
        if(fileName == ""):
            return
        
        shortName = f"{fileName.split('/')[-1]} - {os.path.getsize(fileName)}" 
        print(shortName)
        if(shortName not in self.fileNamesDict):
            self.fileNamesDict[shortName] = fileName
            self.fileListBox.insert(tk.END, shortName)
            
    def removeFile(self):
        if(len(self.fileListBox.curselection()) == 0):
            return
        
        index = self.fileListBox.curselection()[0]
        
        #Get the value at the index
        fileName = self.fileListBox.get(index)
        self.fileListBox.delete(index)
        self.fileNamesDict.pop(fileName)
        
        print(self.fileNamesDict)
        
        
    def run(self):
        #load the parameters
        
        if(len(self.fileNamesDict) <= 1 or self.outPutFile == ""):
            return
        
        t = list(self.fileNamesDict.values())
        
        merge(list(self.fileNamesDict.values()), outputPath=self.outPutFile)
        
        #Update the table view with the new ouput file
        self.tableView.changeData(self.outPutFile)
        
class ScrollableFrame(ttk.Frame):
    def __init__(self, container, csvPath = ""):
        super().__init__(container)
        
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)
        
       
        self.canvasFrame = ttk.Frame(self)
        self.canvasFrame.grid(row=0, column=0, sticky="nsew")
        
        
        canvas = tk.Canvas(self.canvasFrame, background="#C3C3C3")
        
        scrollbar = ttk.Scrollbar(self.canvasFrame, orient="vertical", command=canvas.yview)
        scrollbar_h = ttk.Scrollbar(self.canvasFrame, orient="horizontal", command=canvas.xview)
        self.tableFrame = tableEntry(canvas, csvPath)
        

        self.tableFrame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )        

        canvas.create_window((0, 0), window=self.tableFrame, anchor="nw")

        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.configure(xscrollcommand=scrollbar_h.set)

        scrollbar_h.pack(side="bottom", fill="x")
        scrollbar.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)
        
        
        
    def changeData(self, data):
          
        print("Changing Data", data)
        self.tableFrame.changeData(data)
        
    def addData(self, data):
        self.tableFrame.addData(data)
        
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
        
def merge(filePaths, outputPath="merged.csv"):
    
    
    if(len(filePaths) < 1):
        return None
    
    DFs = []
    columnNames = ["m_index", "std"]
    for path in filePaths:
        DFs.append(pd.read_csv(path))
        columnNames.append(path.split("/")[-1])
        
    
        
    currentDF = DFs[0]
    currentDF.set_index('sequence', inplace=True)
    currentDF.rename(columns={'m_index': columnNames[2]}, inplace=True)
    currentDF.drop(['std'], axis=1, inplace=True)   
    
    for i in range(1, len(DFs)):
        DF = DFs[i]
        DF.set_index('sequence', inplace=True)
        DF.drop(['s_index'], axis=1, inplace=True)
        #Rename the mean col to count
        DF.rename(columns={'m_index': columnNames[i + 2]}, inplace=True)
        currentDF = currentDF.join(DF, how='outer')
        
      
    currentDF.fillna(0, inplace=True)
    print(currentDF)
    
    colSums = currentDF.sum(axis=0)
    currentDF = currentDF.div(colSums, axis=1)
    
    
    means = np.mean(currentDF.values, axis=1)
    stds = np.std(currentDF.values, axis=1)
    
    #Add the means and stds to the currentDF
    
    currentDF['m_index'] = means
    currentDF['s_index'] = stds
    
    total_row = [ 1/ (len(filePaths) * x) for x in colSums.values]
    t =  np.mean(total_row)
    total_row = ["NORMALIZED_ONE_COUNT", t, 0.0] + total_row
    
    dfColumns = columnNames
    dfColumns.insert(0, "sequence")
    
    header_DF = pd.DataFrame([total_row], columns=dfColumns)
    header_DF.set_index('sequence', inplace=True)
    
    
    currentDF.sort_values(by=['m_index'], inplace=True, ascending=False)
    
    totalDF = pd.concat([header_DF, currentDF])
    totalDF.to_csv(outputPath, index=True, columns = columnNames[1:])
    
        
    
    
if __name__ == "__main__":
    app(None)
