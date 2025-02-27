import yaml
import argparse
from runUnifier import merge
import os
import json

#---------------------------------------------------------#
#Script: CLI_runUnifier.py
#Description: Command line interface for the run unifier. Implemented to allow increased usage options
#---------------------------------------------------------#


#---------------------------------------------------------#
#Function: main
#Description: Main function for the run unifier command line interface. This function will parse the command line
#           arguments and call the run unifier with the appropriate arguments.
#Inputs: None
#Outputs: None
#---------------------------------------------------------#
def main():
    args = parser.parse_args()
    
    merge(args.f, args.o)
    

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
    
    if(string.split(".")[-1] not in ["csv"]):
        raise ValueError("File must be a csv file")
    
    return string

#---------------------------------------------------------#
#Function: fileOutput
#Description: Function to validate that the input is a valid file format. If the file is not a csv file, a ValueError is raised
#Inputs: string - the string to validate
#Outputs: string - the string if it is a valid file
#---------------------------------------------------------#
def fileOutput(string):
    #Validate that the input is a valid file format
    
    if(string.split(".")[-1] not in ["csv"]):
        raise ValueError("Ouput File must be a csv file")
    
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
#If the file is being run as the main file, create a parser and parse the arguments
#then call the main function passing the parsed arguments
#---------------------------------------------------------#  
if(__name__ == "__main__"):
    parser = argparse.ArgumentParser(
                    prog='P3ANUT Sequence Counter CLI interface',
                    description='Subsection of the P3ANUT Sequence Counter that allows for command line interface',)

    #forward Files - required, support multiple files
    parser.add_argument('-f',  metavar='files', type=fileInput, nargs='+',
                        required=True, help='The fastq file')
    
    #Output directory to save the output file - required
    parser.add_argument('-o', metavar='ouptut', type=fileOutput, required=True,
                        help='Name of the output File that the output file will be saved to')
    
    
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    
    #python CLI_runUnifier.py -o test.csv -f /Users/ethankoland/Desktop/3rd Year Project/code/testRunData/t/Amino_5.csv /Users/ethankoland/Desktop/3rd Year Project/code/testRunData/t/Amino_5.csv
    main()
