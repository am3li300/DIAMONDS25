import os
import matplotlib.pyplot as plt
import numpy as np

PRC_RANDOM_Y = float("inf")
DISEASE = ""

############################ Internal Helpers ##################################
def plot_auroc(avgFPR, avgRecall, save_dir=None):
    plt.plot(avgFPR, avgRecall, label="Mean ROC curve")
    plt.axvline(0.10, ls='--', color='grey', alpha=0.6)   # t-AUROC cutoff
    plt.plot([0, 1], [0, 1], 'k--', label='Random')
    plt.xlabel("False Positive Rate")
    plt.ylabel("True Positive Rate")
    plt.title("Average ROC Curve Across Folds")
    plt.legend()
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        plt.savefig(os.path.join(save_dir, "{0}_ROC.png".format(DISEASE)))
        plt.close()

    else:
        plt.show()

    # Trapezoidal Rule to get area under curve
    auc_score = np.trapezoid(avgRecall, avgFPR)
    print(f"AUROC score: {auc_score:.4f}")

def plot_auprc(avgRecall, avgPrecision, save_dir=None):
    plt.plot(avgRecall, avgPrecision, label="Mean PRC curve")
    plt.axhline(y=PRC_RANDOM_Y, linestyle='--', color='grey', alpha=0.6, label="Random")
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.title("Average PRC Curve Across Folds")
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.legend()
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        plt.savefig(os.path.join(save_dir, "{0}_PRC.png"))
        plt.close()

    else:
        plt.show()

    # Trapezoidal Rule to get area under curve
    auc_score = np.trapezoid(avgPrecision, avgRecall)
    print(f"AUPRC score: {auc_score:.4f}")

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
    partial_auc = np.trapezoid(tpr_seg, fpr_seg)
    return partial_auc / fpr_max

################################ actually use ##################################

def positives(label_directory, numPos, numGenes, disease, source, method):
    global PRC_RANDOM_Y
    global DISEASE
    """
    -numPos represents the number of true positives; how many genes in the ranking file trying to be recovered (not seeds)
    -numNeg represents the number of true negatives; includes both genes that were never disease genes and disease genes that were used as seeds
    -trueP represents the number of true positives found so far given a threshold; any gene trying to be recovered within the threshold
    -falseP represents the number of false positives found so far given a threshold; any gene not trying to be recovered within the threshold
    """

    numNeg = numGenes - numPos

    PRC_RANDOM_Y = numPos/(numPos+numNeg)
    DISEASE = disease

    recall = [0]*(numGenes + 1) # same as TPR
    FPR = [0]*(numGenes + 1)
    precision = [0]*(numGenes + 1)


    def process_file(file):
        f = open(label_directory + '/' + file, 'r')
        trueP = falseP = 0.0
        threshold = 1
        for line in f:
            if threshold > numGenes:
                break

            try:
                label = int(line.split()[2])
                if label == 1:
                    trueP += 1.0

                else:
                    falseP += 1.0

                recall[threshold] += trueP / numPos
                FPR[threshold] += falseP / numNeg
                precision[threshold] += trueP / threshold
                threshold += 1

            except:
                continue

        p = trueP / numPos
        n = falseP / numNeg
        while threshold <= numGenes:
            recall[threshold] += p
            FPR[threshold] += n
            precision[threshold] += trueP / threshold
            threshold += 1

        f.close()
    
    numFiles = 0
    for file in os.listdir(label_directory):
        numFiles += 1
        process_file(file)

    # get the average for each threshold       
    avgRecall = [0]*(numGenes + 1)
    avgFPR = [0]*(numGenes + 1)
    avgPrecision = [0]*(numGenes + 1)

    for i in range(1, numGenes + 1):
        avgRecall[i] = recall[i]/numFiles
        avgFPR[i] = FPR[i]/numFiles
        avgPrecision[i] = precision[i]/numFiles

    graph_dir = "cross_validation/{0}/graphs/{1}/{2}".format(source, method, DISEASE)
    plot_auroc(avgFPR, avgRecall, graph_dir)
    plot_auprc(avgRecall, avgPrecision, graph_dir)

    print("Average Top-K value when threshold == 100: ", avgRecall[100])
    print("Average Top-K value when threshold == 250: ", avgRecall[250])

    t_auroc = compute_truncated_auroc(avgFPR, avgRecall, fpr_max=0.10)
    print(f"Truncated AUROC (FPR ≤ 0.10): {t_auroc}")

def count_lines(file):
    with open(file, "r") as f:
        numLines = sum(1 for line in f if line.strip())
        return numLines