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

def plot_auroc(avgFPR, avgRecall):
    plt.plot(avgFPR, avgRecall, label="Mean ROC curve")
    plt.plot([0, 1], [0, 1], 'k--', label='Random')
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("Average ROC Curve Across Folds")
    plt.legend()
    plt.show()

    # Trapezoidal Rule to get area under curve
    auc_score = np.trapz(avgRecall, avgFPR)
    print(f"AUROC score: {auc_score:.4f}")

def plot_auprc(avgRecall, avgPrecision):
    plt.plot(avgRecall, avgPrecision, label="Mean PRC curve")
    plt.plot([0, 1], [0, 1], 'k--', label='Random')
    plt.xlabel("Recall")
    plt.ylabel("Precision/TPR")
    plt.title("Average PRC Curve Across Folds")
    plt.legend()
    plt.show()

    # Trapezoidal Rule to get area under curve
    auc_score = np.trapz(avgPrecision, avgRecall)
    print(f"AUPRC score: {auc_score:.4f}")


def main():
    """
    -numPos represents the number of true positives; how many genes in the ranking file trying to be recovered (not seeds)
    -numNeg represents the number of true negatives; includes both genes that were never disease genes and disease genes that were used as seeds
    -trueP represents the number of true positives found so far given a threshold; any gene trying to be recovered within the threshold
    -falseP represents the number of false positives found so far given a threshold; any gene not trying to be recovered within the threshold
    """

    directory = input("Enter the folder path: ")
    partition = 3 if "3-fold" in directory else 2

    numDisease = 163 if "schizophrenia" in directory else 95
    numPos = numDisease // partition
    numNeg = NUM_GENES - numPos

    recall = [0]*(NUM_GENES + 1) # same as TPR
    FPR = [0]*(NUM_GENES + 1)
    precision = [0]*(NUM_GENES + 1)

    numFiles = 0
    for file in os.listdir(directory):
        numFiles += 1
        # sum the TPR and FPR for each threshold across each data set
        f = open(directory + '/' + file, 'r')
        for threshold in range(1, NUM_GENES + 1):
            f.seek(0)
            if threshold == 100:
                print(trueP*1.0/numPos)

            trueP, falseP = count_positives(f, threshold)
            recall[threshold] += trueP*1.0 / numPos
            FPR[threshold] += falseP*1.0 / numNeg
            precision[threshold] += trueP*1.0 / threshold

        f.close()   

    # get the average for each threshold       
    avgRecall = [0]*(NUM_GENES + 1)
    avgFPR = [0]*(NUM_GENES + 1)
    avgPrecision = [0]*(NUM_GENES + 1)
    for i in range(1, NUM_GENES + 1):
        avgRecall[i] = recall[i]/numFiles
        avgFPR[i] = FPR[i]/numFiles
        avgPrecision[i] = precision[i]/numFiles
    

    plot_auroc(avgFPR, avgRecall)

    plot_auprc(avgRecall, avgPrecision)

    print("Average Top-K value when threshold == 100: ", avgRecall[100])
    print("Average Top-K value when threshold == 250: ", avgRecall[250])

if __name__ == "__main__":
    main()