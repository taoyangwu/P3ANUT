import yaml
import argparse
import multiprocessedPairAssembler as pairedassembler
import os
import json

#---------------------------------------------------------#
#File: CLI_pairedAssembler.py
#Description: Command line interface for the paired assembler. Implemented to allow increased usage options. Major optinos
#           are included as cli commands. The remains are accesable through the addition of a yaml file.
#---------------------------------------------------------#


#---------------------------------------------------------#
#Function: main
#Description: Main function for the paired assembler command line interface. This function will parse the command line
#           arguments and call the paired assembler with the appropriate arguments.
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
        trimmedYaml = {k:v.get("Value") for k,v in yamlSettings.get("pairedAssembler", {}).items() if v.get("Value") is not None}
    else:
        yamlSettings = {}
        
    
    
    #Overwrite the default settings with the yaml file
    if(args.radius):
        trimmedYaml["radius"] = args.radius
        
    parsedData = {}
    print(trimmedYaml)
    
    
    pairedassembler.parse(args.forwardFile, args.reverseFile, parsedData, **trimmedYaml)

    
    #Support multiple output formats
    if(args.outputFile):
        with open(args.outputFile, "w") as file:
            json.dump(parsedData, file, indent=4)
    
    
    
        pass


#---------------------------------------------------------#
#Function: dir_path
#Description: Function to check that the input is a valid file path. Raises a FileNotFoundError if the file does not exist
#Inputs: string - the string to check
#Outputs: string - the string if it is a valid file path
#---------------------------------------------------------#
def dir_path(string):
    if os.path.isfile(string):
        return string
    else:
        raise FileNotFoundError(string)
    
#---------------------------------------------------------#
#Function: yaml_path
#Description: Function to check that the input is a valid YAML file and can be loaded. 
#             Raises a FileNotFoundError if the file does not exist
#             Raises a yaml.YAMLError if the file cannot be loaded as a yaml file
#Inputs: string - the string to check
#Outputs: string - the string if it is a valid file path
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
#then call the main function passing the parsed arguments
#---------------------------------------------------------#

if(__name__ == "__main__"):
    parser = argparse.ArgumentParser(
                    prog='Fastq Assembler Pipeline',
                    description='Assemble fastq files into a single file')
    #Input arguments
    
    #forward File - required
    parser.add_argument('-f', '--forwardFile', metavar='forwardFile', type=dir_path,
                        required=True, help='The forward fastq file')
    #Reverse File - optional, support the single assembler
    parser.add_argument('-r', '--reverseFile', metavar='reverseFile', type=dir_path,
                        required=True, help='The reverse fastq file')
    #Output File - optional, default to pairedAssembleDate.json
    parser.add_argument('-o','--outputFile', metavar='outputFile', type=str,
                        default='pairedAssembleOutput.json', help='Name of the output file')
    
    
    
    #Check if the output file already exists, 
    
    
    #Yaml file - optional, automatically create one? or just use defaults?
    parser.add_argument('-y','--yamlFile', metavar='yamlFile', type=yaml_path,
                        help="name of the yaml file that contains the advance settings")


    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    #Overload the default settings with the yaml file
    parser.add_argument('--radius', metavar='radius', type=int,
                        help='The radius of k-mers to use for the assembler')
    parser.add_argument('--scoreOffset', metavar='scoreOffset', type=int,
                        help='The scoreOffset of k-mers to use for the assembler')
    

    
    main()