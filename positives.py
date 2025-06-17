import os
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np

# output/validation_output_labels/schizophrenia
# output/validation_output_labels/rheumatoid_arthritis

NUM_GENES = 11882 # for STRING network

def count_positives(file, threshold):
    trueP = falseP = 0
    for line in file:
        if not threshold:
            break
            
        try:
            label = int(line.split()[2])
            if label == 1:
                trueP += 1

            else:
                falseP += 1

            threshold -= 1

        except:
            continue

    return trueP, falseP

def plot_auroc(averages):
    x = [arr[1] for arr in averages if arr != [0,0]]  # skip 0th index
    y = [arr[0] for arr in averages if arr != [0,0]]
    plt.plot(x, y, label="Mean ROC curve")
    plt.plot([0, 1], [0, 1], 'k--', label='Random')
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("Average ROC Curve Across Folds")
    plt.legend()
    plt.show()

    # Trapezoidal Rule to get area under curve
    auc_score = np.trapz(y, x)
    print(f"AUC score: {auc_score:.4f}")

def main():
    """
    -numPos represents the number of true positives; how many genes in the ranking file trying to be recovered (not seeds)
    -numNeg represents the number of true negatives; includes both genes that were never disease genes and disease genes that were used as seeds
    -trueP represents the number of true positives found so far given a threshold; any gene trying to be recovered within the threshold
    -falseP represents the number of false positives found so far given a threshold; any gene not trying to be recovered within the threshold
    """

    directory = input("Enter the folder path: ")

    numDisease = 163 if "schizophrenia" in directory else 95
    numPos = numDisease // 2
    numNeg = NUM_GENES-numPos

    # [0] True Positives, [1] False Positives
    positives = [[0.0]*2 for _ in range(NUM_GENES + 1)]

    numFiles = 0
    for file in os.listdir(directory):
        numFiles += 1
        # sum the TPR and FPR for each threshold across each data set
        f = open(directory + '/' + file, 'r')
        for threshold in range(1, NUM_GENES + 1):
            f.seek(0)
            trueP, falseP = count_positives(f, threshold)
            positives[threshold][0] += trueP*1.0 / numPos
            positives[threshold][1] += falseP*1.0 / numNeg

        f.close()   
            
    print("Top-K value when threshold == 100: ", positives[100][0])
    print("Top-K value when threshold == 250: ", positives[250][0])

    # get the average for each threshold       
    averages = [[0]*2 for _ in range(NUM_GENES + 1)]
    for i in range(1, NUM_GENES + 1):
        averages[i][0] = positives[i][0]/numFiles
        averages[i][1] = positives[i][1]/numFiles

    plot_auroc(averages)

if __name__ == "__main__":
    main()