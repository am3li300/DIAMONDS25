import pandas as pd

PATH="../data/matrix_small.csv"
CHUNK_SIZE=10**6

with pd.read_csv(PATH, iterator=True, chunksize=CHUNK_SIZE) as reader:
    print(reader.read())
    exit(1)

