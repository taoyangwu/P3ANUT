import tkinter as tk

#Importing all of the frames
from GUI_sequenceCount import sequenceCountFrame
from GUI_pairedAssembler import pairedAssemblerFrame
from runUnifier import runUnifierFrame
from volcanoPlot import volcanoPlotFrame
from upsetPlot import upsetPlotFrame
from rankingPlot import rankingPlotFrame

class app(tk.Tk):
    def __init__(self, controller):
        super().__init__()
        
        
        
        
        self.geometry("800x800")
        
        self.grid_rowconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=5)
        
        self.grid_columnconfigure(0, weight=1)
        
        self.buttonFrame = tk.Frame(self)
        self.buttonFrame.grid(row=0, column=0, sticky="nsew")
        
        self.mainFrame = tk.Frame(self)
        
        self.title("P3ANUT : A Python Pipeline for Phage Analysis through a Normative Unified Toolset")
        
        
        self.pairedAssemblerButton = tk.Button(self.buttonFrame, text="Paired Assembler", command=lambda: self.changeMainFrame("pairedAssembler"))
        self.pairedAssemblerButton.grid(row=0, column=0, sticky="nsew")
        
        
        
        self.sequenceCountButton = tk.Button(self.buttonFrame, text="Sequence Count", command=lambda: self.changeMainFrame("sequenceCount"))
        self.sequenceCountButton.grid(row=0, column=1, sticky="nsew")
        
    
        
        self.runUnifierButton = tk.Button(self.buttonFrame, text="Run Unifier", command=lambda: self.changeMainFrame("runUnifier"))
        self.runUnifierButton.grid(row=0, column=2, sticky="nsew")
        
        self.volcanoPlotButton = tk.Button(self.buttonFrame, text="Volcano Plot", command=lambda: self.changeMainFrame("volcanoPlot"))
        self.volcanoPlotButton.grid(row=0, column=3, sticky="nsew")
        
        self.upsetPlotButton = tk.Button(self.buttonFrame, text="Upset Plot", command=lambda: self.changeMainFrame("upsetPlot"))
        self.upsetPlotButton.grid(row=0, column=4, sticky="nsew")
        
        self.upsetPlotButton = tk.Button(self.buttonFrame, text="Ranking Plot", command=lambda: self.changeMainFrame("rankingPlot"))
        self.upsetPlotButton.grid(row=0, column=5, sticky="nsew")
        
        
        
        self.mainloop()
        
    def changeMainFrame(self, frame):
        self.mainFrame.destroy()
        
        if(frame == "pairedAssembler"):
            self.mainFrame = pairedAssemblerFrame(self)
        elif(frame == "sequenceCount"):
            self.mainFrame = sequenceCountFrame(self)
        elif(frame == "runUnifier"):
            self.mainFrame = runUnifierFrame(self)
        elif(frame == "volcanoPlot"):
            self.mainFrame = volcanoPlotFrame(self)
        elif(frame == "upsetPlot"):
            self.mainFrame = upsetPlotFrame(self)
        elif(frame == "rankingPlot"):
            self.mainFrame = rankingPlotFrame(self)
            
        self.mainFrame.grid(row=1, column=0, sticky="nsew")
        
        
if(__name__ == "__main__"):
    app = app(None)
    
    # app.pairedAssemblerButton.bind("<Button-1>", lambda x: pairedAssemblerFrame(app))
    # app.sequenceCountButton.bind("<Button-1>", lambda x: sequenceCountFrame(app))
    # app.runUnifierButton.bind("<Button-1>", lambda x: runUnifierFrame(app))
    # app.volcanoPlotButton.bind("<Button-1>", lambda x: volcanoPlotFrame(app))
    # app.upsetPlotButton.bind("<Button-1>", lambda x: upsetPlotFrame(app)
