import argparse
from upsetPlot import upsetPlot
import os


#---------------------------------------------------------#
#Function: main
#Description: Main function for the upset plot command line interface. This function will parse the command line
#           arguments and call the upset plot with the appropriate arguments.
#Inputs: parser - the parser object
#Outputs: None
#---------------------------------------------------------#
def main(parser):
    
    args = parser.parse_args()
    
    upsetPlot.fileOutput(set(args.files), args.key, args.output)
    
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
#If the file is being run as the main file, create a parser and parse the arguments
#---------------------------------------------------------#
if(__name__ == "__main__"):
    parser = argparse.ArgumentParser(
                    prog='P3ANUT Upset Plot CLI interface',
                    description='Subsection of the P3ANUT Supset Plot that allows for command line interface',)
    
    parser.add_argument('-k', "--key", metavar='files', type=int,
                        required=True, help='The int representation of the binary key. If 3 files are being compared, the key could be 101 (enter 5). This would find the insertsection of the first and thrid file and excluuding any sequences found in the second file')
    parser.add_argument('-f', "--files", metavar='files', type=fileInput, nargs='+',
                        required=True, help='Files to be compared')
    
    parser.add_argument('-o', "--output", metavar='output', type=str, default="upsetPlot.csv",
                        help='The output file')
    
    main(parser)
    
    