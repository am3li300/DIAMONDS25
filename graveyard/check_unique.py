import os

# cross_validation/validation_rankings/SZ_STRING_adaptive_k
directory_path = input("Enter directory path: ")
k = int(input("Enter k value: "))
flag = int(input("Include scores? 1: Yes, 0: No "))
seen = set()

for fileName in os.listdir(directory_path):
    file = open(os.path.join(directory_path, fileName))
    top_k = tuple([(line.split()[0] if not flag else f"{line.split()[0]} {line.split()[1]}") for line in file][-k:])
    print(top_k)
    if top_k not in seen:
        seen.add(top_k)

    else:
        print(f"File {i} is not unique")

"""
flag = int(input("See distribution? 1: Yes, 0: No "))
if flag:
    file = open(f"data/20_data_drug_{"rheumatoid_arthritis" if disease == "RA" else "schizophrenia"}.txt")
    included = {node[:-1]: 0 for node in file}

    filePath = "cross_validation/" + ("rheumatoid_arthritis" if disease == "RA" else "schizophrenia") \
                + f"/{disease if disease == "RA" else "schizophrenia"}_new_seeds_"

    for i in range(25):
        file = open(f"{filePath}{str(i)}.txt")
        for line in file:
            node = line[:-1]
            included[node] += 1

    print()
    print(included)
    print()

    print("Min:", min(included.values()))
    print("Max:", max(included.values()))
"""