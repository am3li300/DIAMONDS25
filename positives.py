import os
import matplotlib as plt

NUM_GENES = 11882 # for STRING network

def count_positives(file, threshold):
    trueP = 0
    for line in file:
        if not threshold:
            break
            
        try:
            label = int(line.split()[2][:-1])
            trueP += label
            threshold -= 1
        except:
            continue

    return trueP, threshold-trueP

def plot_auroc(averages):
    x = [arr[1] for arr in averages]
    y = [arr[0] for arr in averages]
    plt.plot(x, y)
    plt.show()

def main():
    directory = input("Enter the folder path: ")

    numDisease = 163 if "schizophrenia" in directory else 95
    numPos = numDisease // 2

    # [0] True Positives, [1] False Positives
    positives = [[0.0]*2 for _ in range(NUM_GENES + 1)]

    numFiles = 0
    for file in os.listdir(directory):
        numFiles += 1
        f = open(directory + '/' + file, 'r')
        # sum the TPR and FPR for each threshold across each data set
        for threshold in range(1, NUM_GENES + 1):
            trueP, falseP = count_positives(f, threshold)
            positives[threshold][0] += trueP / numPos
            positives[threshold][1] += falseP / (NUM_GENES - numPos)

        f.close()   
            
    # get the average for each threshold       
    averages = [[0]*2 for _ in range(NUM_GENES + 1)]
    for i in range(1, NUM_GENES + 1):
        averages[i][0] = positives[i][0]/numFiles
        averages[i][1] = positives[i][1]/numFiles

    plot_auroc(averages)

if __name__ == "__main__":
    main()