from scipy.stats import chi2
import re
import os
import sys
from pathlib import Path

#import positional arguments
if len(sys.argv) < 3:
        print("Usage: python LRT.py model1 model2")
        sys.exit(1)

model1 = sys.argv[1]
model2 = sys.argv[2]

#crete basename

ortho = os.path.basename(model1)[:9]

#from the full model result, extract the likelihood, the numbero fo degrees of freedom, and the three values of omega inferred
def full_likelihood(full_model):
        with open(full_model, "r") as file:
                for line in file:
                        if "lnL" in line:
                                line_lnL = line.rstrip()
                        elif "(dN/dS) for branches" in line:
                                line_w = line.rstrip()

        pattern = r"[-+]?\d*\.\d+|\d+"
        #set the following variables as global ones, so other functions can access to them
        global likefull, paramfull, w0, w1, w2
        likefull = float(re.findall(pattern, line_lnL)[2])
        paramfull = int(re.findall(pattern, line_lnL)[1])
        w0 = float(re.findall(pattern, line_w)[0])
        w1 = float(re.findall(pattern, line_w)[1])
        w2 = float(re.findall(pattern, line_w)[2])

#from the reduced model result extract the likelihood
def reduced_likelihood(reduced_model):
        with open(reduced_model, "r") as file:
                for line in file:
                        if "lnL" in line:
                                line_lnL = line.rstrip()
                                break

        pattern = r"[-+]?\d*\.\d+|\d+"
        global likered, paramredu
        likered = float(re.findall(pattern, line_lnL)[2])
        paramredu = int(re.findall(pattern, line_lnL)[1])

#perform the LRT writing the result in a specific file with all the values worthy to be successively analysed 
def LRT(likefull, likered, doffull, dofred, outfile, ortho, w0, w1, w2):
        p_val = "{:.2e}".format(chi2.sf(2*(likefull-likered),(doffull-dofred)))
        significance = "Significant" if float(p_val) <= 0.05 else "Not significant"
        #trait = "TRAIT1" if "short" in directory else  "TRAIT2"
        #traccer = "accelerated" if "acc" in directory else "constrained"
        deltaw = round(w1 - w2, 5)
        with open(outfile, "a") as file:
                file.write(f"{ortho}\t{likefull}\t{likered}\t{p_val}\t{significance}\t{w0}\t{w1}\t{w2}\t{deltaw}\n")

#create the header if the file doesn't exist or is empty
if not Path("LRT_results.tsv").exists() or Path("LRT_results.tsv").stat().st_size == 0:
    with open("LRT_results.tsv", "w") as f:
        f.write("Orthogroup\tLikelihood_full\tLikelihood_reduced\tp_val\tSignificance\tw0\tw1\tw2\tDeltaw\n")

full_likelihood(model1)
reduced_likelihood(model2)
LRT(likefull, likered, paramfull, paramredu, "LRT_results.tsv", ortho, w0, w1, w2)