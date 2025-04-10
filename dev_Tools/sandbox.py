import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats as sp

fileA = "/Users/ethankoland/Desktop/FOR ETHAN/temp.csv"
fileB = "/Users/ethankoland/Desktop/FOR ETHAN/m7man-cona-vs-m7man-bsa.csv"

dataFrameA = pd.read_csv(fileA)
dataFrameB = pd.read_csv(fileB)

print(dataFrameA.columns.values)

print([x for x in dataFrameA.columns.values if x not in ["sequence", "mean", "std"]])
