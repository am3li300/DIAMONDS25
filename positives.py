import os
import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np

NUM_GENES = int(input("Enter number of genes in network: ")) # 11882 for STRING network, 10317 for humanbase

def compute_truncated_auroc(fpr, tpr, fpr_max=0.10):
    """
    fpr, tpr : lists of equal length with values in [0,1]
               They must start at index 0 = threshold 0.
    fpr_max  : right-hand truncation bound (ADAGIO uses 0.10).

    Returns the truncated AUROC rescaled to [0,1].
    """
    # strip off the index-0 padding and make numpy arrays
    fpr = np.asarray(fpr[1:], dtype=float)
    tpr = np.asarray(tpr[1:], dtype=float)

    # only keep the segment fpr <= fpr_max
    mask = fpr <= fpr_max
    fpr_seg = fpr[mask]
    tpr_seg = tpr[mask]

    if not fpr_seg.size or fpr_seg[-1] < fpr_max:
        # find first point where fpr > fpr_max
        idx_hi = np.argmax(fpr > fpr_max)

        fpr_lo, fpr_hi = fpr[idx_hi - 1], fpr[idx_hi]
        tpr_lo, tpr_hi = tpr[idx_hi - 1], tpr[idx_hi]
        slope = (tpr_hi - tpr_lo) / (fpr_hi - fpr_lo)

        tpr_interp = tpr_lo + slope * (fpr_max - fpr_lo)
        fpr_seg = np.concatenate((fpr_seg, [fpr_max]))
        tpr_seg = np.concatenate((tpr_seg, [tpr_interp]))
    
    # prepend the (0,0) start so trapezoid() runs correctly
    fpr_seg = np.concatenate(([0.0], fpr_seg))
    tpr_seg = np.concatenate(([0.0], tpr_seg))

    # rescale
    partial_auc = np.trapz(tpr_seg, fpr_seg)
    return partial_auc / fpr_max

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
    plt.axvline(0.10, ls='--', color='grey', alpha=0.6)   # t-AUROC cutoff
    plt.plot([0, 1], [0, 1], 'k--', label='Random')
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("Average ROC Curve Across Folds")
    plt.legend()
    plt.show()

    # Trapezoidal Rule to get area under curve
    auc_score = np.trapezoid(avgRecall, avgFPR)
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
    auc_score = np.trapezoid(avgPrecision, avgRecall)
    print(f"AUPRC score: {auc_score:.4f}")


def main():
    """
    -numPos represents the number of true positives; how many genes in the ranking file trying to be recovered (not seeds)
    -numNeg represents the number of true negatives; includes both genes that were never disease genes and disease genes that were used as seeds
    -trueP represents the number of true positives found so far given a threshold; any gene trying to be recovered within the threshold
    -falseP represents the number of false positives found so far given a threshold; any gene not trying to be recovered within the threshold
    """

    # cross_validation/validation_output_labels/
    directory = input("Enter the folder path to validation labels: ")

    # 2-fold, 3-fold, etc.
    partition = int(input("Enter number of folds used for validation: "))

    # 163 for schizophrenia, 95 for RA, 114 for humanbase SZ
    numDisease = int(input("Enter total number of seed genes: "))

    # numPos = numDisease // partition
    # clustering louvain: 51; clustering markov: 50; walktrap: 49
    numPos = int(input("Enter number of total positives (genes trying to be recovered): "))
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

    t_auroc = compute_truncated_auroc(avgFPR, avgRecall, fpr_max=0.10)
    print(f"Truncated AUROC (FPR ≤ 0.10): {t_auroc}")


if __name__ == "__main__":
    main()