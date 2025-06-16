disease = input("Check SZ or RA output: ")
filePath = "output/validation_rankings/" + ("rheumatoid_arthritis" if disease == "RA" else "schizophrenia") \
            + f"/{disease}_cross_validation_"

k = int(input("Enter k-value: "))
seen = set()
for i in range(25):
    file = open(f"{filePath}{str(i)}.out")
    top_k = tuple([line.split()[0] for line in file][-k:])
    if top_k not in seen:
        seen.add(top_k)

    else:
        print(f"File {i} is not unique")