import pandas as pd

df = pd.read_csv("debug_volcano_distro.csv")
print(max(df['AvB_Ratio']))