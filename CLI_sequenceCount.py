import yaml
import argparse
import sequenceCounter
import os
import json

#---------------------------------------------------------#
#Script: CLI_sequenceCount.py
#Description: Command line interface for the sequence counter. Implemented to allow increased usage options
#---------------------------------------------------------#

#---------------------------------------------------------#
#Function: main
#Description: Main function for the sequence counter command line interface. This function will parse the command line
#           arguments and call the sequence counter with the appropriate arguments.
#Inputs: None
#Outputs: None
#---------------------------------------------------------#
def main():
    args = parser.parse_args()
    
    #Load the yaml file
    if(args.yamlFile):
        with open(args.yamlFile, 'r') as stream:
            try:
                yamlSettings = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)
                return
        
        #Trim the yaml file to only the settings that are needed
        trimedYAML = {k:v.get("Value") for k,v in yamlSettings.get("sequenceCount", {}).items() if v.get("Value") is not None}
    else:
        yamlSettings = {}
        
    
    
    #Overwrite the default settings with the yaml file
    if(args.method):
        trimedYAML["method"] = args.method
        
    if(args.encoding):
        trimedYAML["encoding"] = args.encoding

    print(trimedYAML)
    
    if(args.file.split(".")[-1] in ["fastq", "fq"]):
        sequenceCounter.parseFastq(args.file, args.parseDNA, args.parseAmino,
                                   basedirectory = args.ouptut, **trimedYAML)
        
    elif(args.file.split(".")[-1] in ["csv"]):
        sequenceCounter.parseCSV(args.file, args.parseDNA, args.parseAmino,
                                  basedirectory = args.ouptut, **trimedYAML)
        
    elif(args.file.split(".")[-1] in ["json"]):
        sequenceCounter.parseJSON(args.file, args.parseDNA, args.parseAmino,
                                   basedirectory = args.ouptut, **trimedYAML)
    


    
    
#---------------------------------------------------------#
#Function: fileInput
#Description: Function to validate that the input is a valid file format
#Inputs: string - the string to validate
#Outputs: string - the string if it is a valid file
#---------------------------------------------------------#
def fileInput(string):
    #Validate that the input is a valid file format
    if not os.path.isfile(string):
        raise FileNotFoundError(string)
    
    if(string.split(".")[-1] not in ["fastq", "fq", "csv"]):
        raise ValueError("File must be a fastq or csv file")
    
    return string

#---------------------------------------------------------#
#Function: dir_path
#Description: Function to check that the input is a valid file path. Raises a FileNotFoundError if the file does not exist
#Inputs: string - the string to validate
#Outputs: string - the string if it is a valid file
#---------------------------------------------------------#
def dir_path(string):
    if os.path.isdir(string):
        return string
    else:
        raise FileNotFoundError(string)
    
#---------------------------------------------------------#
#Function: yaml_path
#Description: Function to check that the input is a valid YAML file and can be loaded.
#             Raises a FileNotFoundError if the file does not exist
#             Raises a yaml.YAMLError if the file cannot be loaded as a yaml file
#Inputs: string - the string to validate
#Outputs: string - the string if it is a valid file
#---------------------------------------------------------#

def yaml_path(string):
    #Check that the file exists and is a valid yaml file
    if os.path.isfile(string):
        with open(string, 'r') as stream:
            try:
                yamlSettings = yaml.safe_load(stream)
                return string
            except yaml.YAMLError as exc:
                print(exc)
                raise exc
    else:
        raise FileNotFoundError(string)

#---------------------------------------------------------#
#If the file is being run as the main file, create a parser and parse the arguments
#---------------------------------------------------------#
if(__name__ == "__main__"):
    parser = argparse.ArgumentParser(
                    prog='P3ANUT Sequence Counter CLI interface',
                    description='Subsection of the P3ANUT Sequence Counter that allows for command line interface',)
    #Input arguments
    #forward File - required
    parser.add_argument('-f',  metavar='file', type=fileInput,
                        required=True, help='The fastq file')
    
    #Output directory to save the output file
    parser.add_argument('-o', metavar='ouptut', type=dir_path, required=True,
                        help='Name of the output Folder that the output file will be saved to')
    
    #Boolean Toggles to determine weither to encode DNA and Amino Acids
    parser.add_argument('--parseDNA', action='store_false',
                        help='Count the DNA sequences with in the file. The default acition is to not count the DNA sequences.')
    parser.add_argument('--parseAmino', action='store_True',
                        help='Count the Amino Acid sequences with in the file. The default acition is to count the Amino Acid sequences.')
    
    #Controlling the encoding, and method
    parser.add_argument('--encoding', metavar='encoding', choices=sequenceCounter.getEncodings(),
                        help='The encoding to use for the sequences. The default is to use the default encoding')
    
    parser.add_argument('--method', metavar='method', choices=sequenceCounter.getMethods(),
                        help='The method to use for the sequences. The default is to use the default method.')
    
    #Loading the advanced settings - YAML File
    parser.add_argument('-y','--yamlFile', metavar='yamlFile', type=yaml_path,
                        help="name of the yaml file that contains the advance settings")
    
    
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    
    main()