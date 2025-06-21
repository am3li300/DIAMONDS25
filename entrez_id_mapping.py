import pandas as pd

filePath = input("Enter file path for gene info: ")
f = open(filePath)

headers = f.readline().split()
mapp = {header:[] for header in headers}







# [entrez gene id 1][entrez gene id 2][posterior prob.]