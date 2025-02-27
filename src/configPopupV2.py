import tkinter as tk
from tkinter import ttk
import yaml
from tkinter.filedialog import askopenfilename,asksaveasfilename
from datetime import datetime
from copy import deepcopy

#---------------------------------------------------------#
#Script: configPopupV2.py
#Description: Script to create a popup window that allows the user to change the configuration settings for a YAML file
#---------------------------------------------------------#

#---------------------------------------------------------#
#Class: app
#Description: Class used to test the popup window
#---------------------------------------------------------#
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
        parameterMenu(self, dataPath, self.parms)
        
    def test(self):
        print("testing")
        print(hex(id(self.parms)))
        print(self.parms)
 
#---------------------------------------------------------#
#Class: parameterMenu
#Description: Class used to create the popup window that allows the user to change the configuration settings for a YAML file
#---------------------------------------------------------#       
class parameterMenu(tk.Toplevel):
    #---------------------------------------------------------#
    #Function: __init__
    #Description: Function to initialize the parameterMenu class
    #Inputs: parent - the parent window
    #        dataPath - the path to the yaml file
    #        parameters - the parameters to be changed
    #        currentScript - the current script, used to improve the user experience by automatically selecting the current script
    #Outputs: None
    def __init__(self, parent, dataPath = "", parameters = {}, currentScript = "pariedAssembler"):
        super().__init__(parent)
        
        #Storing the parameters in the class
        self.dataPath = dataPath
        self.data = {}
        self.scriptSelction = self.data.keys()
        self.params = parameters
        self.currentScript = tk.StringVar(value = currentScript if currentScript != "" else "None")
        
        
        #Statement to validate that the parameters are not being changed
        # print(hex(id(self.params)))
        
        
        # self.grid_rowconfigure(0, weight=1)
        
        
        self.title("Config Change Menu")
        self.titleLabel = tk.Label(self, text="Yaml Config", font=("Arial", 20), justify=tk.CENTER)
        self.titleLabel.grid(column=0, row=0, sticky=tk.W + tk.E)
        self.grid_columnconfigure(0, weight=1)
        
        #-------------#
        
        self.scriptDropdownFrame = tk.Frame(self, background="red")
        self.scriptDropdownFrame.grid(column=0, row=1, sticky=tk.W + tk.E)
        self.scriptDropdownFrame.columnconfigure(0, weight=1)
        self.scriptDropdownFrame.columnconfigure(1, weight=1)
        
        #-------------#
        
        self.scriptDropsdownLabel = tk.Label(self.scriptDropdownFrame, text="Script")
        self.scriptDropsdownLabel.grid(column=0, row=0, sticky=tk.W + tk.E)
        
        #-------------#
        
        self.dataEntriesFrame = ScrollableFrame(self, data = {})
        self.dataEntriesFrame.grid(column=0, row=2, sticky=tk.W + tk.E)
        
        #-------------#
        self.scriptDropdown = tk.OptionMenu(self.scriptDropdownFrame, self.currentScript, [], command=self.optionChanged)
        self.scriptDropdown.grid(column=1, row=0, sticky=tk.W + tk.E)
        
        #-------------#
        self.controllFrame = tk.Frame(self)
        self.controllFrame.grid(column=0, row=3, sticky=tk.W + tk.E)
        
        
        #-------------#
        
        self.loadFrame = tk.LabelFrame(self.controllFrame, text="Load", background="red")
        self.loadFrame.grid(column=0, row=3, sticky=tk.W + tk.E)
        
        self.loadYamlButton = tk.Button(self.loadFrame, text="Load Yaml", command=self.loadYaml)
        self.loadYamlButton.grid(column=0, row=0, sticky=tk.W + tk.E)
        
        #-------------#
        self.saveFrame = tk.LabelFrame(self.controllFrame, text="Save", background="green")
        self.saveFrame.grid(column=0, row=4, sticky=tk.W + tk.E)
        
        self.saveYamlButton = tk.Button(self.saveFrame, text="Save Changes", command=self.save)
        self.saveYamlButton.grid(column=0, row=0, sticky=tk.W + tk.E)
        
        self.saveYamlButton = tk.Button(self.saveFrame, text="Save As", command=self.saveAs)
        self.saveYamlButton.grid(column=1, row=0, sticky=tk.W + tk.E)
        
        #-------------#
        
        self.exitFrame = tk.LabelFrame(self.controllFrame, text="Exit", background="blue")
        self.exitFrame.grid(column=0, row=5, sticky=tk.W + tk.E)        
        
        self.applyButton = tk.Button(self.exitFrame, text="Apply Changes", command=self.applyChanges)
        self.applyButton.grid(column=0, row=0, sticky=tk.W + tk.E)
        
        self.closeButton = tk.Button(self.exitFrame, text="Close Window", command=self.destroy)
        self.closeButton.grid(column=1, row=0, sticky=tk.W + tk.E)
        
        if(self.dataPath != ""):
            self.loadYaml(self.dataPath)
        
        
        # self.pack()
    
    #---------------------------------------------------------#
    #Function: optionChanged
    #Description: Function to change the data entries when the option is changed
    #Inputs: None
    #Outputs: None
    #---------------------------------------------------------#  
    def optionChanged(self, *args):
        print("Option Changed")
        print(self.currentScript.get())
        self.dataEntriesFrame.changeData(self.data.get(self.currentScript.get(), {}))
      
    #---------------------------------------------------------#
    #Function: optionChangedCurrentScript
    #Description: Function to change the data entries when the option is changed. Slight variation of the optionChanged
    #               function is needed to update currentScript
    #Inputs: value - the value of the option
    #Outputs: None
    #---------------------------------------------------------#
    def optionChangedCurrentScript(self, value, *args):
        print("Option Changed")
        self.currentScript.set(value)
        print(self.currentScript.get())
        self.dataEntriesFrame.changeData(self.data.get(self.currentScript.get(), {}))
        
    #---------------------------------------------------------#
    #Function: saveAs
    #Description: Function to save the yaml file as a new file
    #Inputs: None
    #Outputs: None
    #---------------------------------------------------------#
    def saveAs(self):
        filePath = asksaveasfilename(defaultextension=".yaml", filetypes=[("Yaml Files", "*.yaml")])
        
        #If the user cancels the save, return with out creating a new file
        if(filePath == ""):
            return
        
        for key, value in self.dataEntriesFrame.getValues().items():
            
            #Loop through each script or major section within the yaml file
            for i in self.data.keys():
                #Loop through each of the parameters within the script
                if(key in self.data[i].keys()):
                    print("Change", key, value)
                    self.data[i][key]["Value"] = value
                    continue
        
        #Save the new file
        with open(filePath, 'w') as file:
            yaml.dump(self.data, file)
    
    #---------------------------------------------------------#
    #Function: save
    #Description: Function to save the yaml file. Overwrites the current file
    #Inputs: None
    #Outputs: None
    #---------------------------------------------------------#      
    def save(self):
        
        # print("Save", self.currentScript.get())
        # print(self.dataEntriesFrame.getValues())    
        
        for key, value in self.dataEntriesFrame.getValues().items():
            
            for i in self.data.keys():
                if(key in self.data[i].keys()):
                    print("Change", key, value)
                    self.data[i][key]["Value"] = value
                    continue
            # print(key, value)
            # self.data[self.currentScript.get()][key]["Value"] = value
        
        print("Save", self.dataPath)
        with open(self.dataPath, 'w') as file:
            yaml.dump(self.data, file)
      
    #---------------------------------------------------------#
    #Function: loadYaml
    #Description: Function to load the yaml file
    #Inputs: file - the file path to the yaml file
    #Outputs: None
    #---------------------------------------------------------#  
    def loadYaml(self, file = ""):
        
        #clear all options from the dropdown
        self.scriptDropdown['menu'].delete(0, 'end')
        
        #If the file is not provided, ask the user for the file
        if(file == ""):
            file = askopenfilename(filetypes=[("Yaml Files", "*.yaml")])
         
        #Save the file path to the class   
        self.dataPath = file

        #Attempt to load the yaml file
        with open(file, 'r') as stream:
            try:
                self.data = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)
          
        #Load each of the supported scripts within the yaml file      
        self.scriptSelction = list(self.data.keys())
        print("scriptSelection", self.scriptSelction)
        
        #Load the current script, if the current script is not in the list, set the current script to the first script
        self.currentScript = tk.StringVar(value = self.currentScript.get() if self.currentScript.get() == "None" else list(self.scriptSelction)[0]) 
        print("Current Script", self.currentScript.get())
        
        #Add each of the possible scripts to the dropdown
        for script in self.scriptSelction:
            self.scriptDropdown['menu'].add_command(label=script, command=tk._setit(self.currentScript, script, self.optionChangedCurrentScript))
        self.scriptDropdown.configure(textvariable=self.currentScript)

        
        #Update the data entries
        self.dataEntriesFrame.changeData(self.data.get(self.currentScript.get(), {}))
     
     
    #---------------------------------------------------------#
    #Function: applyChanges
    #Description: Function to apply the changes to the local parameters, does not update the orignla yaml file.
    #           Used to test the changes before saving
    #Inputs: None
    #Outputs: None
    #---------------------------------------------------------#
    def applyChanges(self):
        
        #Apply the changes to the local parameters
        for key, value in self.data.items():
            self.params[key] = value
        
        #Update the parameters within the parent window, can not pass a referecene due to memory ownership  
        for key, value in self.dataEntriesFrame.getValues().items():
            
            for i in self.data.keys():
                if(key in self.data[i].keys()):
                    print("Change", key, value)
                    self.data[i][key]["Value"] = value
                    continue

        self.closeWindow()
        
    #---------------------------------------------------------#
    #Function: loadChanges
    #Description: Function to load the changes from the data entries
    #Inputs: None
    #Outputs: None
    #---------------------------------------------------------#
    def loadChanges(self):
        print("Load Changes")
        
        #Attempt to determine the config chage file
        configChangeFile = self.data["global"]["configChangeHistory"].get("Value", "")
        
        #A common set of operations to write the changes to a file, so this is a lambda function
        changesToString = lambda configName, key, newValue, oldValue: f"{configName}-{datetime.now()}-{key}\nNew Value: {newValue}\nOld Value: {oldValue}\n"
        changes = self.dataEntriesFrame.getValues()
        
        changesTxt = ""
        
        #Create a string that contains the changes for each of the parameters that have been changed
        for key, value in changes.items():
            if((newValue := self.data[self.currentScript.get()][key]["Value"]) != value):
                print("Change", key, value)
                self.data[self.currentScript.get()][key]["Value"] = value
                changesTxt += changesToString(self.currentScript.get(), key, value, newValue)
                
        print(changesTxt)
        
        #Write the changes to the config change file
        try:
            configChangeFile = self.data["global"]["configChangeHistory"].get("Value", "")
            if(configChangeFile != ""):
                with open(configChangeFile, 'a') as file:
                    file.write(changesTxt + "\n")
        except Exception as e:
            print("Error Writing to Config Change File", e)
        # print(self.dataEntriesFrame.getValues())
        
        pass
    
    #---------------------------------------------------------#
    #Function: closeWindow
    #Description: Function to close the window
    #Inputs: None
    #Outputs: None
    #---------------------------------------------------------#
    def closeWindow(self):
        print("Close Window")
        self.destroy()
        
    
    
#---------------------------------------------------------#
#Class: errorWindow
#Description: Class used to create a popup window that displays an error message
#---------------------------------------------------------#    
class errorWindow(tk.Toplevel):
    def __init__(self, master = None, message = None):    
        super().__init__(master = master)
        self.title("Error")
        self.geometry("200x200")
        label = tk.Label(self, text = message)
        label.pack()
 
#---------------------------------------------------------#
#Class: helpWindow
#Description: Class used to create a popup window that displays help information
#---------------------------------------------------------#   
class helpWindow(tk.Toplevel):
    #Create a refernce dictionary that contains infomration about each of the scripts and the functions
    def __init__(self, master = None, function = None, script = None):    
        super().__init__(master = master)
        self.title("Help")
        self.geometry("200x200")
        label = tk.Label(self, text = function)
        label.pack()

#---------------------------------------------------------#
#Class: ScrollableFrame
#Description: Class used to create a scrollable frame. Used to contain the data entries and 
#           allow the user to scroll to see all of the data entries
#---------------------------------------------------------#
class ScrollableFrame(ttk.Frame):
    #---------------------------------------------------------#
    #Function: __init__
    #Description: Function to initialize the ScrollableFrame class
    #Inputs: container - the parent window
    #        data - the data to be displayed. A reference to the data is pass to the dataEntries class
    #Outputs: None
    #---------------------------------------------------------#
    def __init__(self, container, data = {}):
        super().__init__(container)
        canvas = tk.Canvas(self)
        scrollbar = ttk.Scrollbar(self, orient="vertical", command=canvas.yview)
        self.scrollable_frame = dataEntries(canvas, data = data)

        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(
                scrollregion=canvas.bbox("all")
            )
        )

        canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")

        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
       
    #---------------------------------------------------------#
    #Function: getValues
    #Description: A wrapper function to get the Values in the dataEntries class
    #Inputs: None
    #Outputs: None
    #---------------------------------------------------------#
    def getValues(self):
        return self.scrollable_frame.getValues()
    
    #---------------------------------------------------------#
    #Function: changeData
    #Description: A wrapper function to change the data in the dataEntries class
    #Inputs: data - the new data to be displayed
    #Outputs: None
    #---------------------------------------------------------#
    def changeData(self, data):
        self.scrollable_frame.changeData(data)
        
        

#---------------------------------------------------------#
#Class: dataEntries
#Description: Class used to create the data entries for the yaml file
#---------------------------------------------------------#
class dataEntries(ttk.Frame):
    #---------------------------------------------------------#
    #Function: __init__
    #Description: Function to initialize the dataEntries class
    #Inputs: parent - the parent window
    #        data - the data to be displayed, passed as a reference to the shared data dictionary
    #Outputs: None
    def __init__(self, parent, data = {}):
        super().__init__(parent)
        
        #Maintain a reference to the data, so that the data can be changed
        self.data = data
        
        self.entries = {}
        
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(0, weight=1)
        
        self.columnconfigure(0, weight=1)
        self.columnconfigure(1, weight=1)
        self.rowconfigure(0, weight=1)
        
        self.createChildFrames()
        
    #---------------------------------------------------------#
    #Function: changeData
    #Description: Function to change the data in the dataEntries class, removes all of the child frames and creates new ones
    #Inputs: data - the new data to be displayed
    #Outputs: None
    #---------------------------------------------------------#
    def changeData(self, data):
        self.data = data
        

        for l in self.grid_slaves():
            l.destroy()
        
        # print("Change Data", self.data)
        self.createChildFrames()
        
    #---------------------------------------------------------#
    #Function: createChildFrames
    #Description: Creates child frames that display the data entries and help buttons
    #Inputs: None
    #Outputs: None
    #---------------------------------------------------------#
    def createChildFrames(self):
        #Function to create the child frames for each of the data entries
        
        for i, (key, value) in enumerate(self.data.items()):
            t = dataEntry(self, key, value)
            self.entries[key] = t
            t.grid(column=0, row=i + 1, sticky=tk.W + tk.E)
         
    #---------------------------------------------------------#
    #Function: getValues
    #Description: Function to get the values from the data entries
    #Inputs: None
    #Outputs: dictionary - the values from the data entries
    #---------------------------------------------------------#   
    def getValues(self):
        #Function to get the values from the data entries
        return {k:v.getValue() for k, v in self.entries.items()}
        
#---------------------------------------------------------#
#Class: dataEntry
#Description: Class used to create the data entry for the yaml file
#---------------------------------------------------------#
class dataEntry(tk.Frame):
    #---------------------------------------------------------#
    #Function: __init__
    #Description: Function to initialize the dataEntry class
    #Inputs: parent - the parent window
    #        key - the key for the data entry
    #        data - the data to be displayed
    #Outputs: None
    #---------------------------------------------------------#
    def __init__(self, parent, key, data):
        
    #Sample Data 
    #   Value : "errorData.json"
    #   Type : 'String'
    #   Description : 'Name or location of the file to append the error results to'
    #   UsedFunctions : 'pairedAssembler'
        super().__init__(parent)
        
        #Save the key and data to the class
        self.key = key
        self.data = data
        self.itemType = data["Type"]
        
        #Setting up the label frame
        self.labelFrame = tk.LabelFrame(self, text=key)
        self.labelFrame.grid_columnconfigure(0, weight=4)
        self.labelFrame.grid_columnconfigure(1, weight=1)
        
        #Establishing the domain of the data entry
        if(self.itemType.lower() == 'boolean'):
            self.value = tk.BooleanVar(value = data["Value"])
            self.valueEntry = booleanEntry(self.value, self.labelFrame)
        elif(self.itemType.lower() == 'int'):
            self.value = tk.IntVar(value = data["Value"])
            self.valueEntry = valueEntry(self.value, self.labelFrame, validChars = '0123456789')
        elif(self.itemType.lower() == 'float'):
            self.value = tk.DoubleVar(value = data["Value"])
            self.valueEntry = valueEntry(self.value, self.labelFrame, validChars = '0123456789.')
        elif(self.itemType.lower() == 'selection'):
            self.value = tk.StringVar(value = data["Value"])
            self.valueEntry = selectionEntry(self.value, data["Options"], self)
        else:
            self.value = tk.StringVar(value = data["Value"])
            self.valueEntry = valueEntry(self.value, self.labelFrame, validChars = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789.-+')
    
            
        self.valueEntry.grid(column=0, row=0, sticky=tk.W + tk.E)
        
        self.helpButton = tk.Button(self.labelFrame, text="Help", command=self.showHelp)
        self.helpButton.grid(column=1, row=0, sticky=tk.W + tk.E)
        
        self.labelFrame.pack()
        
    def showHelp(self):
        helpWindow(self, self.key, self.data["UsedFunctions"])
        
    def getValue(self):
        return self.value.get()
        
 
        
class booleanEntry(tk.Frame):
    def __init__(self, tkVar, parent = None):
        super().__init__(parent)
        self.options = ["True", "False"]
        
        self.tkVar = tkVar
        self.localValue = tk.StringVar(value = "True" if self.tkVar.get() else "False")
        
        self.updateValue = lambda x: self.tkVar.set(True) if x == "True" else self.tkVar.set(False)
        
        
        self.dropdown = tk.OptionMenu(self, self.localValue, *self.options, command=self.updateValue)
        self.dropdown.pack()
        #self.dropdown.configure(command=lambda x: print(x))

        
    def updateValue(self, *args):
        print("update")
        self.dropdown.configure(text=self.localValue.get())
        self.tkVar.set(self.localValue.get())
        
        # self.pack()
        
class selectionEntry(tk.Frame):
    def __init__(self, tkVar, options,  parent = None):
        super().__init__(parent)
        
        self.value = tkVar
        
        self.dropdown = tk.OptionMenu(self, self.value, *options)
        self.dropdown.pack()
        
        # self.pack()
        
    def changeOptions(self, options):
        self.dropdown = tk.OptionMenu(self, self.value, *options)
        self.dropdown.pack()
        
        # self.pack()
    
        
class valueEntry(tk.Frame):
    def __init__(self, tkVar, parent = None, validChars = '0123456789'):
        super().__init__(parent)
        
        self.value = tkVar
        self.validChars = set(list(validChars))
        
        self.entry = tk.Entry(self, textvariable=self.value)
        self.entry.pack()
        
        self.vcmd = (self.register(self.validate),
                '%d', '%i', '%P', '%s', '%S', '%v', '%V', '%W')
        self.entry.configure(validate = 'key', validatecommand = self.vcmd)
        
        # self.pack()
    
    def validate(self, action, index, value_if_allowed, prior_value, text, validation_type, trigger_type, widget_name):
        if(action=='1'):
            if text in self.validChars:
                print(text, self.validChars)
                try:
                    return True
                except ValueError:
                    return False
            else:
                return False
        else:
            return True

    
if(__name__ == "__main__"):
    app(None).mainloop()
