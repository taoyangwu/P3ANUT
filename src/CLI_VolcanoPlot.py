import argparse
from volcanoPlot import supportingLogic
import os

#---------------------------------------------------------#
#Subsection of the P3ANUT Volcano Plot that allows for command line interface
#Script Name: CLI_VolcanoPlot.py
#Description: Command line interface for the volcano plot. 
# Limited options are included as cli commands. The CLI only supports outputing a csv of possible sequences
#---------------------------------------------------------#


#---------------------------------------------------------#
#Function: main
#Description: Main function for the volcano plot command line interface. This function will parse the command line
#           arguments and call the volcano plot with the appropriate arguments.
#Inputs: parser - the parser object
#Outputs: None
#---------------------------------------------------------#
def main(parser):
    args = parser.parse_args()
    
    data = supportingLogic.csvComparision(args.InputA, args.InputB)

    
        
    quadDf = supportingLogic.returnQuadrant(data, args.ratio, args.pvalue, args.quadrant)
  
        
    exportDF = supportingLogic.exportDF(args.InputA, args.InputB, quadDf)

    exportDF.to_csv(args.output)    

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
                    prog='P3ANUT Volcano Plot CLI interface',
                    description='Subsection of the P3ANUT Volcano Plot that allows for command line interface',)
    
    parser.add_argument('-a', "--InputA", metavar='files', type=fileInput,
                        required=True, help='The A file to be used in the comparision')
    parser.add_argument('-b', "--InputB", metavar='files', type=fileInput,
                        required=True, help='The B file to be used in the comparision')
    
    parser.add_argument('-r', "--ratio", metavar='ratio', type=float, default=0.5,
                         help='The ratio to be used in the comparision')
    parser.add_argument('-p', "--pvalue", metavar='pvalue', type=float, default=0.05,
                         help='The pvalue to be used in the comparision')
    
    parser.add_argument('-q', "--quadrant", metavar='quadrant', type=int, default=1, choices=[1,2,3,4])
    
    parser.add_argument('-o', "--output", metavar='output', type=str, default="volcanoPlot.csv",)
    
    main(parser)
    
    #python3.12 CLI_VolcanoPlot.py -a '/Users/ethankoland/Desktop/3rd Year Project/code/testRunData/RunUnified/GAL_LB-Uni.csv' -b '/Users/ethankoland/Desktop/3rd Year Project/code/testRunData/RunUnified/GAL-BSA-Uni.csv' -r 470 -p 0.5 -q 2 -o '/Users/ethankoland/Desktop/3rd Year Project/code/testRunData/RunUnified/volcanoPlot.csv'
