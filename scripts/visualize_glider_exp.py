# visualize_glider_exp.py 

#!/cluster/tufts/cowenlab/.envs/tmap/bin/python
import glob
import argparse
import string
import matplotlib.pyplot as plt

#/cluster/tufts/cowenlab/mgelem01/E-MTAB-5214/tissue_files
parser = argparse.ArgumentParser(description='Visualize ROC results')
parser.add_argument('-r', '--result-path', type=str, required=True,
                    help='Path to folder that contains true/false positive data')

def main(args: argparse.Namespace):
    # Results files will be all files in path ending in ".result"
    result_files = glob.glob(args.result_path + '/*.result')
    
    # Plot data in each file
    for result_file in result_files:
        
        # First line: true positives
        # Second line: false positives
        with open(result_file) as file:
           yvals = file.readline()
           xvals = file.readline()
           file.close()

        #reformat strings
        yvals = format_roc_vals(yvals)
        xvals = format_roc_vals(xvals)

        #plot them
        plot_roc(yvals, xvals)               


def format_roc_vals(roc_vals):
    roc_vals = roc_vals[5:-2]
    roc_vals = roc_vals.split(",")

    for i in range(len(roc_vals)):
        roc_vals[i] = float(roc_vals[i])

    return roc_vals


def plot_roc(yvals, xvals):
    plt.plot(xvals, yvals)
    plt.xlabel("False Positive")
    plt.ylabel("True Positive")
    plt.savefig('/cluster/tufts/cowenlab/mgelem01/E-MTAB-5214/GLIDER_Params_Experiment/glider_exp_results.jpg')


if __name__ == '__main__':
    args = parser.parse_args()
    main(args)