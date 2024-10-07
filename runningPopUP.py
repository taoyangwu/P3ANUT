import tkinter as tk


class runningPopUp(tk.Toplevel):
    def __init__(self, controller):
        super().__init__()
        
        self.controller = controller
        
        # self.geometry("400x400")
        
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)
        
        self.title("Running")
        
        self.text = tk.Text(self)
        self.text.grid(row=0, column=0, sticky="nsew")
        
        self.text.insert(tk.END, "The program is running. Please do not close this window. It will automatically update when the program is finished.")
        
        # self.protocol("WM_DELETE_WINDOW", self.on_closing)
        
    def finishedProgram(self):
        self.text.insert(tk.END, "\n\nThe program has finished. You may now close this window. This popup will close in 5 second.")
        self.text.see(tk.END)
        
        self.after(5000, self.destroy)
        
        

class app(tk.Tk):
    def __init__(self, controller):
        super().__init__()
        
        self.geometry("500x500")
        self.parms = {}
        
        self.parameterButton = tk.Button(self, text="Parameters", command=self.showParameters)
        self.parameterButton.pack()
        
        self.testButton = tk.Button(self, text="Test", command=self.test)
        self.testButton.pack()

    def showParameters(self):
        dataPath = ""
        print(hex(id(self.parms)))
        self.running = runningPopUp(self)
        
    def test(self):
        self.running.finishedProgram()
        
    
if __name__ == "__main__":
    app(None).mainloop()