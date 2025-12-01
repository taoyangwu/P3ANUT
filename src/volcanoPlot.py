import tkinter as tk
from tkinter import ttk
import os

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from matplotlib.patches import Rectangle

import numpy as np
import pandas as pd
from scipy import stats as sp

from runningPopUP import runningPopUp

class app(tk.Tk):
    def __init__(self, controller):
        super().__init__()
        
        self.geometry("800x600")
        
        frame = volcanoPlotFrame(self)
        frame.pack(fill="both", expand=True)
        
        self.mainloop()

class volcanoPlotFrame(tk.Frame):
    def __init__(self, controller):
        super().__init__(controller )
    
        
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=3)
        self.grid_columnconfigure(1, weight=1)
        
        self.currentRatio = 5
        self.currentPvalue_log = 1.3
        self.currentPvalue = 0.05
        
        
        self.fileA = tk.StringVar(value="VolcanoPlot/mergedSorted copy.csv")
        self.fileB = tk.StringVar(value="VolcanoPlot/7 copy.csv")
        
        self.graphFrame = tk.LabelFrame(self, text="Graph", background="#DDDEE5")
        self.graphFrame.grid(row=0, column=0, sticky="nsew")
        
        self.parametersFrame = tk.LabelFrame(self, text="Parameters", background="#C3C3C3")
        self.parametersFrame.grid(row=0, column=1, sticky="nsew")
        self.parametersFrame.grid_columnconfigure(0, weight=1)
        


        # self.data = supportingLogic.csvComparision("VolcanoPlot/mergedSorted copy.csv", "VolcanoPlot/7 copy.csv")
        
        # print(self.data)
        
        self.maxRatio, self.minRatio = 1, 0
        self.maxPvalue, self.minPvalue = 1, 0
        
        self.ratioEntryFrame = tk.LabelFrame(self.parametersFrame, text="Ratio Settings", background="#7C8594")
        self.ratioEntryFrame.grid(row=0, column=0, sticky="nsew")
        self.ratioEntryFrame.grid_columnconfigure(0, weight=1)
        
        self.ratioEntry = tk.Entry(self.ratioEntryFrame, validate="key", validatecommand=(self.register(volcanoPlotFrame.validateRatio), '%P'))
        self.ratioEntry.grid(row=0, column=0, sticky="nsew")
        self.ratioEntry.insert(0, str(self.currentRatio))
        self.ratioEntry.bind("<Return>", lambda event: self.updateRatio(float(self.ratioEntry.get())))
        self.ratioEntry.bind("<FocusOut>", lambda event: self.updateRatio(float(self.ratioEntry.get())))
        
        #self.ratioExtInput =
        
        #update the ratio slider range
        
        
        self.pValueEntryFrame = tk.LabelFrame(self.parametersFrame, text="P-Value Settings", background="#6C8193")
        self.pValueEntryFrame.grid(row=1, column=0, sticky="nsew")
        self.pValueEntryFrame.grid_columnconfigure(0, weight=1)
        
        self.pValueEntry = tk.Entry(self.pValueEntryFrame, validate="key", validatecommand=(self.register(volcanoPlotFrame.validatePvalue), '%P'))
        self.pValueEntry.grid(row=0, column=0, sticky="nsew")
        self.pValueEntry.insert(0, str(self.currentPvalue))
        self.pValueEntry.bind("<Return>", lambda event: self.updatePvalue(float(self.pValueEntry.get())))
        self.pValueEntry.bind("<FocusOut>", lambda event: self.updatePvalue(float(self.pValueEntry.get())))
        
        #self.pValue
        
        self.graph = None
        
        self.fileAFrame = tk.LabelFrame(self.parametersFrame, text="1st File", background="#7C8594")
        self.fileAFrame.grid(row=2, column=0, sticky="nsew")
        self.fileAFrame.grid_columnconfigure(0, weight=1)
        
        self.fileAButton = tk.Button(self.fileAFrame, text="1st File", command=self.setfileA)
        self.fileAButton.grid(row=0, column=0, sticky="nsew")
        
        self.fileBFrame = tk.LabelFrame(self.parametersFrame, text="2nd File", background="#6C8193")
        self.fileBFrame.grid(row=3, column=0, sticky="nsew")    
        self.fileBFrame.grid_columnconfigure(0, weight=1)
        
        self.fileBButton = tk.Button(self.fileBFrame, text="2nd File", command=self.setfileB)
        self.fileBButton.grid(row=0, column=0, sticky="nsew")
        
        self.exportFrame = tk.LabelFrame(self.parametersFrame, text="Export", background="#7C8594")
        self.exportFrame.grid(row=4, column=0, sticky="nsew")
        self.exportFrame.grid_columnconfigure(0, weight=4)
        self.exportFrame.grid_columnconfigure(1, weight=1)
        
        self.quadrantList = ["Q1 - Red", "Q2 - Blue", "Q3 - Green", "Q4 - Yellow"]
        self.exportQuarant = tk.StringVar()
        self.exportQuarant.set(self.quadrantList[0])
        self.exportQuarantMenu = tk.OptionMenu(self.exportFrame, self.exportQuarant, *self.quadrantList)
        self.exportQuarantMenu.grid(row=0, column=0, sticky="nsew")
        
        
        self.exportBaseText = tk.Label(self.exportFrame, text="Export Base Name")
        self.exportBaseText.grid(row=1, column=0, sticky="nsew")
        self.exportBaseText.config(background="#7C8594")
        
        self.exportBaseName = tk.StringVar(value="A")
        self.exportBaseName.trace("w", lambda name, index, mode, sv=self.exportBaseName: self.updateExportBaseName(sv.get()))
        self.exportBaseNameText = tk.Entry(self.exportFrame, textvariable=self.exportBaseName)
        self.exportBaseNameText.grid(row=1, column=1, sticky="nsew")
        
        
        self.exportButton = tk.Button(self.exportFrame, text="Export", command=self.export)
        self.exportButton.grid(row=0, column=1)
        
        self.createGraphFrame = tk.LabelFrame(self.parametersFrame, text="Create Graph", background="#6C8193")
        self.createGraphFrame.grid(row=5, column=0, sticky="nsew")
        self.createGraphFrame.grid_columnconfigure(0, weight=1)
        
        self.createGraphButton = tk.Button(self.createGraphFrame, text="Create Graph", command=self.createGraph)
        self.createGraphButton.grid(row=0, column=0, sticky="nsew")
    
        
        self.exportGraphButton = tk.Button(self.createGraphFrame, text="Export Graph", command=self.exportGraph)
        self.exportGraphButton.grid(row=1, column=0, sticky="nsew")
        
        self.distrabutionFrame = tk.Button(self.createGraphFrame, text=f"Export {self.exportBaseNameText.get()}1v{self.exportBaseNameText.get()}2 Distrabution", command=self.createDistrabution)
        self.distrabutionFrame.grid(row=2, column=0, sticky="nsew")
      
    #Print out the value within the ratio entry box  
    def temp(self):
        print(self.ratioEntry.get())
        
    def updateRatio(self, new_val):
        print("new Ratio val: ", new_val)
        self.currentRatio = float(new_val)
        
        if(self.graph is not None):
            
            self.graph.setRatio(self.currentRatio)
        

    def exportGraph(self):
        fileName = tk.filedialog.asksaveasfilename(filetypes=[("PNG Files", "*.png")])
        
        if(fileName == ""):
            return
        
        if(self.graph is not None):
            self.graph.exportGraph(fileName)
            
    def createDistrabution(self):
        
        
        if(self.fileA.get() == "" or self.fileB.get() == ""):
            return
        
        fileName = tk.filedialog.asksaveasfilename(filetypes=[("Txt Files", "*.txt")])
        
        supportingLogic.createDistrabution(self.fileA.get(), self.fileB.get(), fileName)
        
        
        
        
        
        
    def updatePvalue(self, new_val):
        print("new pValue val: ", new_val)
        self.currentPvalue = new_val
        self.currentPvalue_log = float(-np.log10(float(self.currentPvalue)))
        
        if(self.graph is not None):
            
            self.graph.setPvalue(self.currentPvalue_log, new_val)
        
    def setfileA(self):
        filename = tk.filedialog.askopenfilename(initialdir = "/",
                                            title = "Select a File",
                                            filetypes = [("CSV Files", "*.csv")])
        
        print(filename)
        
        #Validate the file that it contains the correct columns
        if(filename == ""):
            return
        
        self.fileA.set(filename)
        self.fileAButton["text"] = self.fileA.get().split("/")[-1]
        
        
    def setfileB(self):
        filename = tk.filedialog.askopenfilename(initialdir = "/",
                                            title = "Select a File",
                                            filetypes = [("CSV Files", "*.csv")])
        
        if(filename == ""):
            return
        
        self.fileB.set(filename)
        self.fileBButton["text"] = self.fileB.get().split("/")[-1]
        
    def createGraph(self):
        
        if(self.fileA.get() == "" or self.fileB.get() == ""):
            return
        
        if(self.graph is not None):
            self.graph.destroy()
        
        runPopUp = runningPopUp(self)
        self.data = supportingLogic.csvComparision(self.fileA.get(), self.fileB.get())
        
        print(self.data)
        
        self.maxRatio = max(self.data["AvB_Ratio"])
        self.minRatio = min(self.data["AvB_Ratio"])
        self.maxPvalue = max(self.data["-log10(P-Value)"])
        self.minPvalue = min(self.data["-log10(P-Value)"])
        
        # self.currentPvalue = float(self.pValueEntry.get())
        # self.currentRatio = float(self.ratioEntry.get())
        
        self.updateRatio(self.ratioEntry.get())
        self.updatePvalue(self.pValueEntry.get())
        
        # self.ratioSlider.config(from_=self.minRatio, to=self.maxRatio, resolution=(self.maxRatio - self.minRatio)/100)
        # self.pValueSlider.config(from_=self.minPvalue, to=self.maxPvalue)
        
        # #Set the current slider values to the average of the max and min values
        # self.ratioSlider.set(self.currentRatio)
        # self.pValueSlider.set(self.currentPvalue)
        
        #Remove the old graph
        if(self.graph is not None):
            self.graph.destroy()
        
        self.graph = moveableLine(self.graphFrame, self.data, [self.fileA.get(), self.fileB.get()], self.exportBaseName.get(),self.currentPvalue_log, self.currentRatio)
        self.graph.pack()
        
        runPopUp.finishedProgram()
        
        self.updateRatio(self.ratioEntry.get())
        self.updatePvalue(self.pValueEntry.get())
        
    def export(self):
        
        filename = tk.filedialog.asksaveasfilename(filetypes=[("CSV Files", "*.csv")])
        
        if(filename == ""):
            return
        
        quadrant = self.quadrantList.index(self.exportQuarant.get()) + 1
        
        quadDf = supportingLogic.returnQuadrant(self.data, self.currentRatio, self.currentPvalue_log, quadrant)
        
        exportDF = supportingLogic.exportDF(self.fileA.get(), self.fileB.get(), quadDf)
       
        quadDf.drop([x for x in quadDf.columns.values if x not in ["sequence", "m_index", "s_index"]], axis=1, inplace=True)
        
        exportDF.to_csv(filename)
        
    def validateRatio(P):
        return P.isdigit() or P == "" or P.replace(".", "", 1).isdigit()
    
    def validatePvalue(P):
        return P.isdigit() or P == "" or P.replace(".", "", 1).isdigit()
    
    def updateExportBaseName(self, new_val):
        print("new Export Base Name val: ", new_val)
        self.exportBaseName = str(new_val)
        self.distrabutionFrame.config(text=f"Show {self.exportBaseName}1v{self.exportBaseName}2 Distrabution")
        
        
        
class moveableLine(tk.Frame):
    def __init__(self, controller, data, filenames, exportBaseName, pValue = 0.5, ratio = 0.5, displayResolution = 0.01, distabutionBins = 10, 
                 **kwargs):
        super().__init__(controller)
        
        self.data = data
        self.pValue = pValue
        self.pValue_linear = 0.05
        self.ratio = ratio
        self.displayResolution = displayResolution
        self.distabutionBins = distabutionBins
        self.fileNames = filenames
        self.exportBaseName = exportBaseName
        
        #Change this later to mathc the change in the data
        self.maxRatio = max(self.data["AvB_Ratio"])
        self.minRatio = min(self.data["AvB_Ratio"])
        self.maxPvalue = max(self.data["-log10(P-Value)"])
        self.minPvalue = min(self.data["-log10(P-Value)"])
        
        
        self.q1 = Rectangle((self.ratio, self.pValue), self.maxRatio - self.ratio, self.maxPvalue - self.pValue,
                            fill=True, linewidth=2, color = "red", alpha=0.5)
        self.q2 = Rectangle((self.ratio, self.minPvalue), self.maxRatio - self.ratio, self.pValue,
                            fill=True, linewidth=2, color = "yellow",  alpha=0.5)
        self.q3 = Rectangle((self.minRatio, self.pValue), self.ratio, self.maxPvalue - self.pValue,
                            fill=True, linewidth=2, color = "blue",  alpha=0.5)
        self.q4 = Rectangle((self.minRatio, self.minPvalue), self.ratio, self.pValue,
                            fill=True, linewidth=2, color = "green",  alpha=0.5)
        
        
        
        self.fig = Figure(figsize=(6, 4), dpi=100)
        self.figureCanvas = FigureCanvasTkAgg(self.fig, self)
        
        #Set the title for the figure
        self.fig.suptitle("Volcano Plot")
        
        subFigs = self.fig.subplots(2,2, gridspec_kw={'height_ratios': [5, 1], 'width_ratios': [1, 5]})
        self.ax1 = subFigs[0,0]
        self.ax2 = subFigs[0,1]
        self.ax3 = subFigs[1,0]
        self.ax4 = subFigs[1,1]
        
        self.ax1.set_xlim(0, len(self.fileNames)+1)
        
        trimmedFileNames = [os.path.basename(i) for i in self.fileNames]
        
        # for i in range(len(self.fileNames)):
        #     self.ax1.text(i + 1, 0.5, f"{chr(ord('A') + i )} - {trimmedFileNames[i]}", ha="center", va="center", fontsize=8, rotation=90)
        self.ax1.text(1, 0.5, f"{self.exportBaseName}1 - {trimmedFileNames[0]}", ha="center", va="center", fontsize=8, rotation=90)
        self.ax1.text(2, 0.5, f"{self.exportBaseName}2 - {trimmedFileNames[1]}", ha="center", va="center", fontsize=8, rotation=90)     
        
        self.ax1.set_yticklabels([])
        self.ax1.set_yticks([])
        
        self.ax1.set_xticklabels([])
        self.ax1.set_xticks([])
        
        self.ax1.set(ylabel='File Names')
        self.ax1.yaxis.label.set_size(10)
        
        
       #self.drawCurrentValues()
        
        #Main Volcano Plot
        roundedRatio = np.round(self.data["AvB_Ratio"], 1)
        self.ax2.scatter( roundedRatio, self.data["-log10(P-Value)"], c="black", linewidth=0.5, )
        self.yLine = self.ax2.axhline(y=0.5, color='r', linestyle='-')
        self.xLine = self.ax2.axvline(x=0.5, color='r', linestyle='-')
    
        self.ax2.set_xscale('log')
        

        
        
        self.ax2.add_patch(self.q1)
        self.ax2.add_patch(self.q2)
        self.ax2.add_patch(self.q3)
        self.ax2.add_patch(self.q4)
        
        self.ax2.set_xlim(self.minRatio,self.maxRatio)
        self.ax2.set_ylim(self.minPvalue,self.maxPvalue)
        self.ax2.set_xlabel(f"{self.exportBaseName}1v{self.exportBaseName}2_Ratio")
        self.ax2.set_ylabel("-log10(P-Value)")
        self.ax2.yaxis.tick_right()
        self.ax2.yaxis.set_label_position("right")
        
        
        
        #Count Distrabution
        
        self.countDistabution = supportingLogic.quarantCounts(self.data, self.ratio, self.pValue)
        
        self.ax3_T1 = self.ax3.text(0.75, 0.75, self.countDistabution[0], ha="center", va="center", color="black", fontsize=6)
        self.ax3_T2 = self.ax3.text(0.25, 0.75, self.countDistabution[1], ha="center", va="center", color="black", fontsize=6)
        self.ax3_T3 = self.ax3.text(0.25, 0.25, self.countDistabution[2], ha="center", va="center", color="black", fontsize=6)
        self.ax3_T4 = self.ax3.text(0.75, 0.25, self.countDistabution[3], ha="center", va="center", color="black", fontsize=6)
        
        countRect1 = Rectangle((0.5,0.5), 0.5,0.5 , color='red', alpha=0.5)
        countRect2 = Rectangle((0,0.5), 0.5,0.5 , color='blue', alpha=0.5)
        countRect3 = Rectangle((0.5,0), 0.5,0.5 , color='yellow', alpha=0.5)
        countRect4 = Rectangle((0,0), 0.5,0.5 , color='green', alpha=0.5)
        
        self.ax3.add_patch(countRect1)
        self.ax3.add_patch(countRect2)
        self.ax3.add_patch(countRect3)
        self.ax3.add_patch(countRect4)
        
        self.ax3.xaxis.set_ticklabels([])
        self.ax3.get_yaxis().set_visible(False)
    
        self.ax3.set_xlabel("Quadrant Counts")
        
        
        
        
        
        self.fig.tight_layout()
        self.fig.subplots_adjust(wspace=0, hspace=0)
        
        self.figureCanvas.get_tk_widget().pack()
        
    def drawCurrentValues(self):
        #Clear the ax4
        self.ax4.clear()
        self.ax4.axis('off')    
        self.ax4.set_xlim(0, 1)
        self.ax4.set_yticklabels([])
        self.ax4.set_yticks([])
        
        self.ax4.set_xticklabels([])
        self.ax4.set_xticks([])
        
        self.ax4.text(1, 0.1, f"P-Value {self.pValue_linear} ({round(self.pValue,2)})", ha="right", va="top", fontsize=10)
        self.ax4.text(1, 0.4, f"{self.exportBaseName}1v{self.exportBaseName}2_Ratio {round(self.ratio,2)}", ha="right", va="top", fontsize=10)
        
    def exportGraph(self, fileName):
        self.fig.savefig(fileName)
        
    def setData(self, newData, displayResolution = 0.01):
        self.ax2.clear()
        self.ax2.plot(newData, c="black", linewidth=0.5)
        self.figureCanvas.draw()
        
    def calcAxisDistrabution(self, data, bins = 10):
        self.ax2.hist(data, bins = bins)
        self.figureCanvas.draw()
        
    def setPvalue(self, new_val, textValue = 0.05):
        self.pValue = float(new_val)
        self.pValue_linear = textValue
        print("new pValue val: ", new_val)


        self.yLine.set_ydata([self.pValue])
        
        self.drawCurrentValues()
        self.drawQuadrantSquares()
        self.calculateQuadrant()
        self.figureCanvas.draw()
        
        
    def setRatio(self, new_val):
        self.ratio = float(new_val)
        print("new Ratio val: ", new_val)

        self.xLine.set_xdata([self.ratio])
        self.drawCurrentValues()
        self.drawQuadrantSquares()
        self.calculateQuadrant()
        self.figureCanvas.draw()
        
    def drawQuadrantSquares(self):
        self.q1.set_xy((self.ratio, self.pValue))
        self.q1.set_width(self.maxRatio - self.ratio)
        self.q1.set_height(self.maxPvalue - self.pValue)
        
        self.q2.set_xy((self.ratio, self.minPvalue))
        self.q2.set_width(self.maxRatio - self.ratio)
        self.q2.set_height(self.pValue - self.minPvalue)
        
        self.q3.set_xy((self.minRatio, self.pValue))
        self.q3.set_width(self.ratio - self.minRatio)
        self.q3.set_height(self.maxPvalue - self.pValue)
        
        self.q4.set_xy((self.minRatio, self.minPvalue))
        self.q4.set_width(self.ratio - self.minRatio)
        self.q4.set_height(self.pValue - self.minPvalue)
        
        
        
    def calculateQuadrant(self):
        self.countDistabution = supportingLogic.quarantCounts(self.data, self.ratio, self.pValue)
        
        self.ax3_T1.set_text(self.countDistabution[0])
        self.ax3_T2.set_text(self.countDistabution[1])
        self.ax3_T3.set_text(self.countDistabution[2])
        self.ax3_T4.set_text(self.countDistabution[3])
        
    def export(self, quardant = 1):
        pass
        
        
class supportingLogic:


    def  csvComparision(fileA, fileB):
        
        
        dataFrameA = pd.read_csv(fileA, index_col=0)
        dataFrameA.drop([x for x in dataFrameA.columns.values if x not in ["sequence", "m_index", "s_index"]], axis=1, inplace=True)
        
        #Check in the NORMALIZED_ONE_COUNT row is within the file
        if("NORMALIZED_ONE_COUNT" in dataFrameA.index):
            dfA_oneCount = dataFrameA.loc["NORMALIZED_ONE_COUNT"][0]
            dataFrameA.drop("NORMALIZED_ONE_COUNT", inplace=True)
        else:
            dfA_oneCount = dataFrameA["m_index"].min()
        
        dataFrameB = pd.read_csv(fileB, index_col=0)
        dataFrameB.drop([x for x in dataFrameB.columns.values if x not in ["sequence", "m_index", "s_index"]], axis=1, inplace=True)
        
        if("NORMALIZED_ONE_COUNT" in dataFrameB.index):
            dfB_oneCount = dataFrameB.loc["NORMALIZED_ONE_COUNT"][0]
            dataFrameB.drop("NORMALIZED_ONE_COUNT", inplace=True)
        else:
            dfB_oneCount = dataFrameB["m_index"].min()
        
        dfa_ColumnCount = 3
        dfb_ColumnCount = 3
        
        
        joined = dataFrameA.join(dataFrameB, how='outer', lsuffix='_a', rsuffix='_b')
        
        joined.fillna(0, inplace=True)
        
        sumA, sumB = joined.iloc[:, 0].sum(), joined.iloc[:, 2].sum()
        
        [statistic, pValue] = sp.ttest_ind_from_stats(joined['m_index_a'], joined['s_index_a'], dfa_ColumnCount, joined['m_index_b'], joined['s_index_b'], dfb_ColumnCount, equal_var=False, alternative='greater')
        
        
        joined['m_index_a'] = joined['m_index_a'].replace(0, dfA_oneCount)
        joined['m_index_b'] = joined['m_index_b'].replace(0, dfB_oneCount)
        
        joined['AvB_Ratio'] = (joined['m_index_a'] / sumA) / (joined['m_index_b'] / sumB)
        joined['-log10(P-Value)'] = -(np.log10(pValue + 1e-10))
        
        #Drop columns
        joined.drop(['s_index_a', 's_index_b', 'm_index_a', 'm_index_b'], axis=1, inplace=True)
        
        joined.to_csv("debug_volcano.csv")
        
        return joined

    def quarantCounts(df, ratio, pvalue):
        # df1T = df1.
        q1 = df[(df['AvB_Ratio'] >= ratio) & (df['-log10(P-Value)'] >= pvalue)]
        q2 = df[(df['AvB_Ratio'] < ratio) & (df['-log10(P-Value)'] >= pvalue)]
        q3 = df[(df['AvB_Ratio'] < ratio) & (df['-log10(P-Value)'] < pvalue)]
        q4 = df[(df['AvB_Ratio'] >= ratio) & (df['-log10(P-Value)'] < pvalue)]
        
        return [len(q1), len(q2), len(q3), len(q4)]

    def returnQuadrant(df, ratio, pvalue, quadrant = 1):
        if(quadrant == 1):
            return df[(df['AvB_Ratio'] >= ratio) & (df['-log10(P-Value)'] >= pvalue)]
        elif(quadrant == 2):
            return df[(df['AvB_Ratio'] < ratio) & (df['-log10(P-Value)'] >= pvalue)]
        elif(quadrant == 3):
            return df[(df['AvB_Ratio'] < ratio) & (df['-log10(P-Value)'] < pvalue)]
        elif(quadrant == 4):
            return df[(df['AvB_Ratio'] >= ratio) & (df['-log10(P-Value)'] < pvalue)]
        else:
            raise ValueError('Invalid Quadrant Number. Must be 1, 2, 3, or 4.')
        
    def trimmedDF(csv, quarantDF):
        dataFrameA = pd.read_csv(csv)
        dataFrameA.set_index('sequence', inplace=True)
        return dataFrameA.loc[quarantDF.index]

    def createDisto(list, maxValue = 1, minValue = 0, bins = 10):
        disto = np.zeros(bins)
        counts = np.zeros(bins)
        
        for i in list:
            t = int((i - minValue) / (maxValue - minValue) * bins)
            t = min(t, bins-1)
            disto[t] += 1
            counts[0:t] += 1
            
        maxDisto = np.max(disto)
        disto = disto / maxDisto
        counts = counts / counts[0] if counts[0] != 0 else counts
        
        axes = (np.arange(bins) * maxValue)/bins + 1/(2*bins)  
            
        return disto, counts, axes

    def exportDF(fileA, fileB,  quadrantDF):
        
        dataFrameA = pd.read_csv(fileA)
        dataFrameB = pd.read_csv(fileB)
        
        # print(dataFrameA)
        
        dataFrameA.set_index('sequence', inplace=True)
        dataFrameB.set_index('sequence', inplace=True)
        
        #Drop all columns except for m_index and s_index
        dataFrameA.drop([x for x in dataFrameA.columns.values if x not in ["m_index", "s_index"]], axis=1, inplace=True)
        dataFrameB.drop([x for x in dataFrameB.columns.values if x not in ["m_index", "s_index"]], axis=1, inplace=True)
        
        joined = dataFrameA.join(dataFrameB, how='outer', lsuffix='_a', rsuffix='_b')
        
        joined.fillna(0, inplace=True)
        
        joined['m_index'] = joined['m_index_a'] + joined['m_index_b']
        
        
        joined.drop(['s_index_a', 's_index_b', 'm_index_a', 'm_index_b'], axis=1, inplace=True)
        joined.insert(1, 's_index', 0)
        
        trimmed = joined.loc[quadrantDF.index]
        trimmed.sort_values(by='m_index', inplace=True, ascending=False)
        print(trimmed)
        
        
        return trimmed
    
    def createDistrabution(fileA, fileB, fileName):
        joined = supportingLogic.csvComparision(fileA, fileB)
        joined.to_csv("debug_volcano_distro.csv")
        
        ratioCount = {}
        
        start, end, step = 0, 10, 0.1
        stepDecimal = len(str(step).split(".")[-1])
        rolling_count, rolling_above, rolling_below = 0, 0, 0
        for i in np.linspace(start, end, int((end - start)/step) + 1):
            
            count_withinRange_belowPvalue = len(joined[(joined['AvB_Ratio'] >= i) & (joined['AvB_Ratio'] < i + step) & (joined['-log10(P-Value)'] >= -np.log10(0.05))])
            count_withinRange_abovePvalue = len(joined[(joined['AvB_Ratio'] >= i) & (joined['AvB_Ratio'] < i + step) & (joined['-log10(P-Value)'] < -np.log10(0.05))])
            
            count_withinRange = count_withinRange_belowPvalue + count_withinRange_abovePvalue
            ratioCount[float(i)] = {"belowPvalue": count_withinRange_belowPvalue, "abovePvalue": count_withinRange_abovePvalue, "total": count_withinRange}
            rolling_count += count_withinRange
            rolling_below += count_withinRange_belowPvalue
            rolling_above += count_withinRange_abovePvalue
        
        count_withinRange_belowPvalue = len(joined[(joined['AvB_Ratio'] > i + step) & (joined['-log10(P-Value)'] >= -np.log10(0.05))])
        count_withinRange_abovePvalue = len(joined[(joined['AvB_Ratio'] > i + step) & (joined['-log10(P-Value)'] < -np.log10(0.05))])
        
        count_withinRange = count_withinRange_belowPvalue + count_withinRange_abovePvalue
        
        ratioCount[float(end)] = {"belowPvalue": count_withinRange_belowPvalue, "abovePvalue": count_withinRange_abovePvalue, "total": count_withinRange}
        rolling_count += count_withinRange
        rolling_below += count_withinRange_belowPvalue
        rolling_above += count_withinRange_abovePvalue
            
            
        with open(fileName, 'w') as f:
            f.write("Ratio_Range, Total_Count, Rolling_Count, Count_Below_Pvalue, Count_Above_Pvalue\n")
             
            for key, value in ratioCount.items():
                roundedKey = round(key, stepDecimal)
                if(key == end):
                    f.write(f"{roundedKey}-{roundedKey}+, {value["total"]}, {rolling_count}, {value["belowPvalue"]}, {value["abovePvalue"]}\n")
                else:
                    f.write(f"{roundedKey}-{round(key + step, stepDecimal)}, {value["total"]}, {rolling_count}, {value["belowPvalue"]}, {value["abovePvalue"]}\n")
                rolling_count -= ratioCount[key]["total"]
            
        pass
    
    
    
    
    

    
    
        
        
        
    


def main():
    app(None)

if(__name__ == "__main__"):
    main()
