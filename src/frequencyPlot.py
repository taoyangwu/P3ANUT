import tkinter as tk
from tkinter import ttk
import os

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from matplotlib.patches import Rectangle
from matplotlib.pyplot import Polygon

import numpy as np
import pandas as pd
from scipy import stats as sp
import bisect

from runningPopUP import runningPopUp

class app(tk.Tk):
    def __init__(self, controller):
        super().__init__()
        
        self.geometry("800x600")
        
        frame = frequencyPlotFrame(self)
        frame.pack(fill="both", expand=True)
        
        self.mainloop()

class frequencyPlotFrame(tk.Frame):
    def __init__(self, controller):
        super().__init__(controller )
    
        
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=3)
        self.grid_columnconfigure(1, weight=1)
        
        self.currentSlope = 1
        self.currentB = 0
        self.currentCount = 100
        self.inlcudeFile1 = True
        self.inlcudeFile2 = False
        self.logAxis = False
        
        self.selectedColor = "#00AA00"
        self.unselectedColor = "#CCCCCC"
        
        self.fileA = tk.StringVar(value="VolcanoPlot/mergedSorted copy.csv")
        self.fileB = tk.StringVar(value="VolcanoPlot/7 copy.csv")
        
        self.graphFrame = tk.LabelFrame(self, text="Graph", background="#DDDEE5")
        self.graphFrame.grid(row=0, column=0, sticky="nsew")
        
        self.parametersFrame = tk.LabelFrame(self, text="Parameters", background="#C3C3C3")
        self.parametersFrame.grid(row=0, column=1, sticky="nsew")
        self.parametersFrame.grid_columnconfigure(0, weight=1)
        


        # self.data = supportingLogic.csvComparision("VolcanoPlot/mergedSorted copy.csv", "VolcanoPlot/7 copy.csv")
        
        # print(self.data)
        
        parm_row = 0
        self.slopeFrame = tk.LabelFrame(self.parametersFrame, text="Slope(M) of Line", background="#7C8594")
        self.slopeFrame.grid(row=parm_row, column=0, sticky="nsew")
        self.slopeFrame.grid_columnconfigure(0, weight=1)
        
        self.slopeEntry = tk.Entry(self.slopeFrame, validate="key", validatecommand=(self.register(frequencyPlotFrame.validateSlope), '%P'))
        self.slopeEntry.grid(row=0, column=0, sticky="nsew")
        self.slopeEntry.insert(0, str(self.currentSlope))
        self.slopeEntry.bind("<Return>", lambda event: self.updateSlope())
        self.slopeEntry.bind("<FocusOut>", lambda event: self.updateSlope())
        
        #self.ratioExtInput =
        
        #update the ratio slider range
        
        parm_row += 1
        # self.bValueEntryFrame = tk.LabelFrame(self.parametersFrame, text="", background="#6C8193")
        # self.bValueEntryFrame.grid(row=parm_row, column=0, sticky="nsew")
        # self.bValueEntryFrame.grid_columnconfigure(0, weight=1)
        
        # self.bValueEntry = tk.Entry(self.bValueEntryFrame, validate="key", validatecommand=(self.register(frequencyPlotFrame.validateB), '%P'))
        # self.bValueEntry.grid(row=0, column=0, sticky="nsew")
        # self.bValueEntry.insert(0, str(self.currentB))
        # self.bValueEntry.bind("<Return>", lambda event: self.updateB(float(self.bValueEntry.get())))
        # self.bValueEntry.bind("<FocusOut>", lambda event: self.updateB(float(self.bValueEntry.get())))
        
        self.numberOfPointsEntryFrame = tk.LabelFrame(self.parametersFrame, text="Number of Points", background="#6C8193")
        self.numberOfPointsEntryFrame.grid(row=1, column=0, sticky="nsew")
        self.numberOfPointsEntryFrame.grid_columnconfigure(0, weight=1)
        self.numberOfPointsEntryFrame.grid_columnconfigure(1, weight=5)
        
        self.percentOrCountLabel = tk.OptionMenu(self.numberOfPointsEntryFrame, tk.StringVar(value="%"), "%", "#")
        self.percentOrCountLabel.grid(row=0, column=0, sticky="nsew")
        
        self.numberOfPointsEntry = tk.Entry(self.numberOfPointsEntryFrame, validate="key", validatecommand=(self.register(frequencyPlotFrame.validateCutoff), '%P'))
        self.numberOfPointsEntry.grid(row=0, column=1, sticky="nsew")
        self.numberOfPointsEntry.insert(0, str(self.currentCount))
        
        
        #self.pValue
        
        self.graph = None
        
        parm_row += 1
        self.fileAFrame = tk.LabelFrame(self.parametersFrame, text="1st File", background="#7C8594")
        self.fileAFrame.grid(row=parm_row, column=0, sticky="nsew")
        self.fileAFrame.grid_columnconfigure(0, weight=1)
        
        self.fileAButton = tk.Button(self.fileAFrame, text="1st File", command=self.setfileA)
        self.fileAButton.grid(row=0, column=0, sticky="nsew")
        
        parm_row += 1
        self.fileBFrame = tk.LabelFrame(self.parametersFrame, text="2nd File", background="#6C8193")
        self.fileBFrame.grid(row=parm_row, column=0, sticky="nsew")    
        self.fileBFrame.grid_columnconfigure(0, weight=1)
        
        self.fileBButton = tk.Button(self.fileBFrame, text="2nd File", command=self.setfileB)
        self.fileBButton.grid(row=0, column=0, sticky="nsew")
        
        parm_row += 1
        self.primaryFileFrame = tk.LabelFrame(self.parametersFrame, text="Which Frame to Count", background="#7C8594")
        self.primaryFileFrame.grid(row=parm_row, column=0, sticky="nsew")
        self.primaryFileFrame.grid_columnconfigure(0, weight=1)
        self.primaryFileFrame.grid_columnconfigure(1, weight=1)
        
        self.file1Select = tk.Button(self.primaryFileFrame, text="Frame 1", command=lambda: self.setFileToCount(file1=True))
        self.file1Select.grid(row=0, column=0, sticky="nsew")
        self.file1Select.config(foreground=self.selectedColor, activeforeground=self.selectedColor)
        
        self.file2Select = tk.Button(self.primaryFileFrame, text="Frame 2", command=lambda: self.setFileToCount(file2=True))
        self.file2Select.grid(row=0, column=1, sticky="nsew")
        self.file2Select.config(foreground=self.unselectedColor, activeforeground=self.unselectedColor)
        
        parm_row  += 1
        self.xAxisFrame = tk.LabelFrame(self.parametersFrame, text="X Axis Scale", background="#6C8193")
        self.xAxisFrame.grid(row=parm_row, column=0, sticky="nsew")
        self.xAxisFrame.grid_columnconfigure(0, weight=1)
        self.xAxisFrame.grid_columnconfigure(1, weight=1)
        
        self.linearScaleButton = tk.Button(self.xAxisFrame, text="Linear", command=lambda: self.setAxis(False))
        self.linearScaleButton.grid(row=0, column=0, sticky="nsew")
        #Set the background color of the button to reflect it is selected
        self.linearScaleButton.config(foreground=self.selectedColor, activeforeground=self.selectedColor)
        
        self.logScaleButton = tk.Button(self.xAxisFrame, text="Log", command=lambda: self.setAxis(True))
        self.logScaleButton.grid(row=0, column=1, sticky="nsew")
        self.logScaleButton.config(foreground= self.unselectedColor, activeforeground=self.unselectedColor)
        
        
        
        
        parm_row += 1
        self.exportFrame = tk.LabelFrame(self.parametersFrame, text="Export", background="#7C8594")
        self.exportFrame.grid(row=parm_row, column=0, sticky="nsew")
        self.exportFrame.grid_columnconfigure(0, weight=4)
        self.exportFrame.grid_columnconfigure(0, weight=1)
        
        self.aboveBelowList = ["Bellow - Red", "Above - Blue"]
        self.exportAboveBelow = tk.StringVar()
        self.exportAboveBelow.set(self.aboveBelowList[0])
        self.exportAboveBelowMenu = tk.OptionMenu(self.exportFrame, self.exportAboveBelow, *self.aboveBelowList)
        self.exportAboveBelowMenu.grid(row=0, column=0, sticky="nsew")
        
        
        self.exportBaseText = tk.Label(self.exportFrame, text="Export Base Name")
        self.exportBaseText.grid(row=1, column=0, sticky="nsew")
        self.exportBaseText.config(background="#7C8594")
        
        self.exportBaseName = tk.StringVar(value="A")
        self.exportBaseNameText = tk.Entry(self.exportFrame, textvariable=self.exportBaseName)
        self.exportBaseNameText.grid(row=1, column=1, sticky="nsew")
        
        
        self.exportButton = tk.Button(self.exportFrame, text="Export", command=self.export)
        self.exportButton.grid(row=0, column=1)
        
        parm_row += 1
        self.createGraphFrame = tk.LabelFrame(self.parametersFrame, text="Create Graph", background="#6C8193")
        self.createGraphFrame.grid(row=parm_row, column=0, sticky="nsew")
        self.createGraphFrame.grid_columnconfigure(0, weight=1)
        
        self.createGraphButton = tk.Button(self.createGraphFrame, text="Create Graph", command=self.createGraph)
        self.createGraphButton.grid(row=0, column=0, sticky="nsew")
    
        
        self.exportGraphButton = tk.Button(self.createGraphFrame, text="Export Graph", command=self.exportGraph)
        self.exportGraphButton.grid(row=1, column=0, sticky="nsew")


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
        
    def updateSlope(self):
        self.currentSlope = float(self.slopeEntry.get())
        
        if(self.graph is not None):
            
            self.graph.slope = self.currentSlope
            self.graph.drawAX2(self.inlcudeFile1, self.inlcudeFile2, self.percentOrCount, self.currentCount, self.logAxis)
            self.graph.drawAX3()
        
        
        
    def updateB(self, new_val):
        print("new pValue val: ", new_val)
        self.currentB = float(new_val)
        
        if(self.graph is not None):
            
            self.graph.moveCutoffLine(self.currentSlope, new_val)
        
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
        
    def setAxis(self, logScale):
        self.logAxis = logScale
        print("Log Scale: ", logScale)
        
        #Set the button backgrounds
        if(logScale):
            print("Setting log scale")
            self.logScaleButton.config(foreground=self.selectedColor, activeforeground=self.selectedColor)
            self.linearScaleButton.config(foreground= self.unselectedColor, activeforeground=self.unselectedColor)
        else:
            print("Setting linear scale")
            self.linearScaleButton.config(foreground=self.selectedColor, activeforeground=self.selectedColor)
            self.logScaleButton.config(foreground= self.unselectedColor, activeforeground=self.unselectedColor)
        
        if(self.graph is not None):
            self.graph.drawAX2(self.inlcudeFile1, self.inlcudeFile2, self.percentOrCount, self.currentCount, self.logAxis)
            
    def setFileToCount(self, file1 = None, file2 = None):
        if(file1 is not None):
            print("Toggling File 1")
            self.inlcudeFile1 = not self.inlcudeFile1 or not self.inlcudeFile2
            
            if(self.inlcudeFile1):
                self.file1Select.config(foreground=self.selectedColor, activeforeground=self.selectedColor)
            else:
                self.file1Select.config(foreground=self.unselectedColor, activeforeground=self.unselectedColor)
            
        if(file2 is not None):
            print("Toggling File 2")
            self.inlcudeFile2 = not self.inlcudeFile2 or not self.inlcudeFile1
            
            if(self.inlcudeFile2):
                self.file2Select.config(foreground=self.selectedColor, activeforeground=self.selectedColor)
            else:
                self.file2Select.config(foreground=self.unselectedColor, activeforeground=self.unselectedColor)
                
        if(self.graph is not None):
            print("Updating Graph include file1:", self.inlcudeFile1, " include file2:", self.inlcudeFile2)
            self.graph.drawAX2(self.inlcudeFile1, self.inlcudeFile2, self.percentOrCount, self.currentCount, self.logAxis)
            self.graph.drawAX3()
        
    def createGraph(self):
        
        self.currentCount = float(self.numberOfPointsEntry.get())
        self.percentOrCount = self.percentOrCountLabel.cget("text")
        self.currentSlope = float(self.slopeEntry.get())
        
        if(self.fileA.get() == "" or self.fileB.get() == ""):
            return
        
        if(self.graph is not None):
            self.graph.destroy()
        
        self.data = supportingLogic.csvComparision(self.fileA.get(), self.fileB.get())
        
        self.percentOrCount = self.percentOrCountLabel.cget("text")
        
        #Remove the old graph
        if(self.graph is not None):
            self.graph.destroy()
        
        self.graph = graphFrame(self.graphFrame, self.data, [self.fileA.get(), self.fileB.get()], self.exportBaseName.get(),
                                  self.currentSlope, self.currentB, self.currentCount, self.percentOrCount, self.inlcudeFile1, self.inlcudeFile2)
        self.graph.pack()
        # self.updateB(self.bValueEntry.get())
        
    def export(self):
        
        filename = tk.filedialog.asksaveasfilename(filetypes=[("CSV Files", "*.csv")])
        
        if(filename == ""):
            return
        
        quadrant = self.aboveBelowList.index(self.exportAboveBelow.get()) + 1
        
        quadDf = supportingLogic.returnQuadrant(self.data, self.currentSlope, self.currentB, quadrant)
        
        exportDF = supportingLogic.exportDF(self.fileA.get(), self.fileB.get(), quadDf)
       
        quadDf.drop([x for x in quadDf.columns.values if x not in ["sequence", "mean", "std"]], axis=1, inplace=True)
        
        exportDF.to_csv(filename)
        
    def validateSlope(P):
        return P.isdigit() or P == "" or P.replace(".", "", 1).isdigit()
    
    def validateCutoff(P):
        return P.isdigit() or P == "" or P.replace(".", "", 1).isdigit()
    
        
        
        
class graphFrame(tk.Frame):
    def __init__(self, controller, data, filenames, exportBaseName, slope = 0.5, b = 0.5,  points = 100, percentOrCount = "%",
                 inludeFrameA = True, includeFrameB = False, **kwargs):
        super().__init__(controller)
        
        self.data = data
        self.slope = slope
        self.b = b
        self.fileNames = filenames
        self.exportBaseName = exportBaseName
        self.points = points
        self.percentOrCount = percentOrCount
        
        
        
        self.fig = Figure(figsize=(6, 4), dpi=100)
        self.figureCanvas = FigureCanvasTkAgg(self.fig, self)
        
        #Set the title for the figure
        self.fig.suptitle("Volcano Plot")
        
        subFigs = self.fig.subplots(2,2, gridspec_kw={'height_ratios': [5, 1], 'width_ratios': [1, 5]})
        self.ax1 = subFigs[0,0]
        self.ax2 = subFigs[0,1]
        self.ax3 = subFigs[1,0]
        self.ax4 = subFigs[1,1]
        
        self.drawAX1()
        self.drawAX2()
        self.drawAX3()
        self.drawAX4()
        
        self.fig.tight_layout()
        self.fig.subplots_adjust(wspace=0, hspace=0)
        
        self.figureCanvas.get_tk_widget().pack()
        
        
        self.initialized = True
        
    def redraw(self, logAxis = True, includeFrameA = True, includeFrameB = False, numberOfPoints = 100, slope = 1, b = 0.5, countOrPercent = "%"):
        pass
        
    def drawAX1(self):
        
        self.ax1.clear()
        
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
        
    def drawAX2(self, includeFile1 = True, includeFile2 = False, percentOrCount = "%", points = 100, x_axis_log = False):
        
        self.ax2.clear()
        
        #Main Volcano Plot
        self.x,self.y = supportingLogic.gatherScatterData(self.data, includeFile1, includeFile2, percentOrCount, points)
        self.sc = self.ax2.scatter( self.x,self.y, c="black", linewidth=0.5, )
        #self.ax2.plot([0, self.points], [self.b, self.points * self.slope + self.b], c="green", linewidth=1)
        # self.ax2.plot([0, self.points], [self.b, (self.points * self.slope) + self.b], c="green", linewidth=1)
        x_line = np.arange(0, np.round(self.points), 1)
        y_line = self.slope * x_line + self.b
        self.ax2.plot(x_line, y_line, c="green", linewidth=1)
        
        self.ax2.set_xlabel(f"{self.exportBaseName}1v{self.exportBaseName}2_Ratio")
        self.ax2.set_ylabel("-log10(P-Value)")
        self.ax2.yaxis.tick_right()
        self.ax2.yaxis.set_label_position("right")
        
        if(x_axis_log):
            self.ax2.set_xscale('log')
            
            
        self.annot = self.ax2.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                    bbox=dict(boxstyle="round", fc="w"),
                    arrowprops=dict(arrowstyle="->"))
        self.annot.set_visible(False)
        
        self.fig.canvas.mpl_connect("motion_notify_event", self.hover)
        
        self.figureCanvas.draw()
        
        # if(self.b == 0):
            
        #     corner1 = Polygon([[0, 0], [self.points, 2 * self.points], [0, self.points]], closed=True, fill=True, color="red", alpha=0.5)
        #     corner2 = Polygon([[0, 0], [self.points, self.points], [self.points, 0]], closed=True, fill=True, color="blue", alpha=0.5)
            
        #     self.ax2.add_patch(corner1)
        #     self.ax2.add_patch(corner2)
        # else:
            
        #     corner1 = Polygon([[0, 0], [self.points, self.points * self.slope + self.b], [0, self.points]], closed=True, fill=True, color="red", alpha=0.5)
        #     corner2 = Polygon([[0, 0], [self.points, self.points * self.slope + self.b], [self.points, 0]], closed=True, fill=True, color="blue", alpha=0.5)
            
        #     self.ax2.add_patch(corner1)
        #     self.ax2.add_patch(corner2)
        
    def hover(self, event):
        vis = self.annot.get_visible()
        if event.inaxes == self.ax2:
            cont, ind = self.sc.contains(event)
            if cont:
                self.update_annot(ind)
                self.annot.set_visible(True)
                self.fig.canvas.draw_idle()
            else:
                if vis:
                    self.annot.set_visible(False)
                    self.fig.canvas.draw_idle()
        
    def update_annot(self, ind):
        
        point_index = ind["ind"][0]
    
        pos = self.sc.get_offsets()[point_index]
        self.annot.xy = pos
        
        text = f"{self.data['fileA_index'][point_index][0]}"
        self.annot.set_text(text)
        # annot.get_bbox_patch().set_facecolor(cmap(norm(c[ind["ind"][0]])))
        self.annot.get_bbox_patch().set_alpha(0.4)
        
    def drawAX3(self):
        
        self.ax3.clear()
        
        countT1, countT2 = supportingLogic.aboveBelowCounts(self.x, self.y, self.slope, self.b)
        print("Counts T1:", countT1, " Counts T2:", countT2)
        
        self.ax3_T1 = self.ax3.text(0.25, 0.75, countT1, ha="center", va="center", color="black", fontsize=6)
        self.ax3_T2 = self.ax3.text(0.75, 0.25, countT2, ha="center", va="center", color="black", fontsize=6)
        
        corner1 = Polygon([[0, 0], [1, 1], [0, 1]], closed=True, fill=True, color="red", alpha=0.5)
        corner2 = Polygon([[0, 0], [1, 1], [1, 0]], closed=True, fill=True, color="blue", alpha=0.5)
        
        self.ax3.add_patch(corner1)
        self.ax3.add_patch(corner2)
        
        self.ax3.xaxis.set_ticklabels([])
        self.ax3.get_yaxis().set_visible(False)
    
        self.ax3.set_xlabel("Counts")
        
        self.figureCanvas.draw()
        
    
    def drawAX4(self):
        #Clear the ax4
        self.ax4.clear()
        self.ax4.axis('off')    
        self.ax4.set_xlim(0, 1)
        self.ax4.set_yticklabels([])
        self.ax4.set_yticks([])
        
        self.ax4.set_xticklabels([])
        self.ax4.set_xticks([])
        
        self.ax4.text(1, 0.1, f"P-Value {self.slope}", ha="right", va="top", fontsize=10)
        self.ax4.text(1, 0.4, f"X-Intercept {self.b}", ha="right", va="top", fontsize=10)
        
        
    
    def exportGraph(self, fileName):
        self.fig.savefig(fileName)
        
        
      
        
    def drawQuadrantSquares(self):
        pass
        
        

        
    def export(self, quardant = 1):
        pass
        
        
class supportingLogic:


    def  csvComparision(fileA, fileB):
        
        
        fileA_index = []
    
        with open(fileA, 'r') as f:
            for i, line in enumerate(f):
                line_split = line.strip().split(',')
                if(line_split[0] in ["sequence" or "NORMALIZED_ONE_COUNT"]  ):
                    continue
                fileA_index.append((line_split[0], float(line_split[1])))
                
        fileA_index.sort(key=lambda x: x[1], reverse=True)
        fileA_dict = {seq: [rank + 1, freq] for rank, (seq, freq) in enumerate(fileA_index)}
        
                
        fileB_index = []
    
        with open(fileB, 'r') as f:
            for i, line in enumerate(f):
                line_split = line.strip().split(',')
                if(line_split[0] in ["sequence" or "NORMALIZED_ONE_COUNT"]  ):
                    continue
                fileB_index.append((line_split[0], float(line_split[1])))
                
        fileB_index.sort(key=lambda x: x[1], reverse=True)
        fileB_Dictionary = {seq: [rank + 1, freq] for rank, (seq, freq) in enumerate(fileB_index)}
        
        return {"fileA_index": fileA_index, "fileB_index": fileB_index,
                "fileA_Dictionary" : fileA_dict, "fileB_Dictionary": fileB_Dictionary}
        
    
    def gatherScatterData(data, fileA = True, fileB = False, percentOrCount = "%", maskValues = 100):
        
        points = [[], []]
        
        if(percentOrCount == "%" and False):
            if(fileA and not fileB):
                x = data["fileA_index"]
                 
                min_i = min(len(x), len(data["fileB_index"]))
                for i, (seq, freq) in enumerate(x):
                    if(i >= min_i):
                        break
                    elif(freq < maskValues / 100):
                        min_i = i
                        break
                
                for i, (seq, _) in enumerate(data["fileA_index"][:min_i]):
                    y = data["fileB_Dictionary"].get(seq, -1)
                    if(y == -1):
                        continue
                    
                    points[0].append(i + 1)
                    points[1].append(y[0])
            elif(not fileA and fileB):
                y = data["fileB_index"]
                freqs = [data["fileB_index"][seq][0] for seq in y]
                y_index = supportingLogic.smallest_greater_index(freqs, maskValues / 100)
                
                for i, (seq, _) in enumerate(data["fileB_index"][:y_index]):
                    x = data["fileA_Dictionary"].get(seq, -1)
                    if(x == -1):
                        continue
                    
                    points[0].append(x[0])
                    points[1].append(i + 1)
            elif(fileA and fileB):
                x = data["fileA_index"]
                freqs_x = [data["fileA_index"][seq][0] for seq in x]
                x_index = supportingLogic.smallest_greater_index(freqs_x, maskValues / 100)
                
                y = data["fileB_index"]
                freqs_y = [data["fileB_index"][seq][0] for seq in y]
                y_index = supportingLogic.smallest_greater_index(freqs_y, maskValues / 100)
                
                top_x = data["fileA_index"][:x_index]
                top_y = data["fileB_index"][:y_index]
                
                top_x_sequences = set([seq for seq, _ in top_x])
                top_y_sequences = set([seq for seq, _ in top_y])
                
                combined_sequences = top_x_sequences.union(top_y_sequences)
                for seq in combined_sequences:
                    x = data["fileA_Dictionary"].get(seq, -1)
                    y = data["fileB_Dictionary"].get(seq, -1)
                    if(x == -1 or y == -1):
                        continue
                    points[0].append(x[0])
                    points[1].append(y[0])
            else:
                raise ValueError('Invalid Option must include at least one file')
                
        # elif(percentOrCount == "#"):
        else:
            if(fileA and not fileB):
                
                #Get the top N sequences from file A
                print(maskValues)
                x = data["fileA_index"][:int(maskValues)]
                
                #Loop through and get the corresponding value from file B
                for i, (seq, _) in enumerate(x):
                    
                    #Safe guard against missing sequences
                    y = data["fileB_Dictionary"].get(seq, -1)
                    if(y == -1):
                        continue
                    
                    #Append the values to the points list
                    points[0].append(i + 1)
                    points[1].append(y[0])
                
            elif(not fileA and fileB):
                
                #Get the top N sequences from file B
                y = data["fileB_index"][:int(maskValues)]
                
                #Loop through and get the corresponding value from file A
                for i, (seq, _) in enumerate(y):
                    
                    #Safe guard against missing sequences
                    x = data["fileA_Dictionary"].get(seq, -1)
                    if(x == -1):
                        continue
                    
                    #Append the values to the points list
                    points[0].append(x[0])
                    points[1].append(i + 1)
                    
            elif(fileA and fileB):
                
                
                
                top_x = data["fileA_index"][:int(maskValues)]
                top_y = data["fileB_index"][:int(maskValues)]
                
                top_x_sequences = set([seq for seq, _ in top_x])
                top_y_sequences = set([seq for seq, _ in top_y])
                
                combined_sequences = top_x_sequences.union(top_y_sequences)
                for seq in combined_sequences:
                    x = data["fileA_Dictionary"].get(seq, -1)
                    y = data["fileB_Dictionary"].get(seq, -1)
                    if(x == -1 or y == -1):
                        continue
                    points[0].append(x[0])
                    points[1].append(y[0])
                
            else:
                raise ValueError('Invalid Option must include at least one file')
        
        print("Points gathered: ", len(points[0]))
        return points

    def aboveBelowCounts(x, y, slope, b):
        above, below = 0, 0
        
        x_arr = np.array(x)
        y_arr = np.array(y)
        
        y_line = slope * x_arr + b
        
        y_sign = np.sign(y_arr-y_line)
        
        above = np.sum(y_sign > 0) + np.count_nonzero(y_sign == 0)
        below = np.sum(y_sign < 0)
        
        return above, below
        

    def returnQuadrant(df, slope, b, above = True):
        
        x = np.array(df['FileA_Index'])
        y = x * slope + b
        
        if(above):
            return 0
        elif(not above):
            return 1
        else:
            raise ValueError('Invalid Option must be above (True) or below (False)')


    def exportDF(fileA, fileB,  quadrantDF):
        
        dataFrameA = pd.read_csv(fileA)
        dataFrameB = pd.read_csv(fileB)
        
        # print(dataFrameA)
        
        dataFrameA.set_index('sequence', inplace=True)
        dataFrameB.set_index('sequence', inplace=True)
        
        joined = dataFrameA.join(dataFrameB, how='outer', lsuffix='_a', rsuffix='_b')
        
        joined.fillna(0, inplace=True)
        
        joined['mean'] = joined['mean_a'] + joined['mean_b']
        joined.drop(['std_a', 'std_b', 'mean_a', 'mean_b'], axis=1, inplace=True)
        joined.insert(1, 'std', 0)
        
        trimmed = joined.loc[quadrantDF.index]
        trimmed.sort_values(by='mean', inplace=True, ascending=False)
        print(trimmed)
        
        
        return trimmed
    
    def smallest_greater_index(data, value, ascending=True):
        """
        Return index of the smallest element in sorted 'data' that is > value.
        If no element is greater, returns len(data).
        Works with list/tuple/numpy array/pandas Series.
        ascending=True for ascending-sorted data (default).
        """
        # normalize to a Python sequence for bisect
        seq = data if isinstance(data, (list, tuple)) else list(data)

        if not seq:
            return 0

        if ascending:
            # bisect_right gives insertion point after any equals -> first element > value
            return bisect.bisect_right(seq, value)
        else:
            # for descending order, flip sign to reuse bisect on ascending data
            neg_seq = [-x for x in seq]
            return bisect.bisect_right(neg_seq, -value)
    
    
    
    
    
    
    

    
    
        
        
        
    


def main():
    app(None)

if(__name__ == "__main__"):
    main()
