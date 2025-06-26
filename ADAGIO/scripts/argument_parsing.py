# I want to read in command line arguments.
#
# Date:
#
# Revisions:
#
# To do:
#
#
################################################################################

import itertools
import sys
import os
import argparse

# I want an argument parser that parses 3 arguments.
# 1. The name of the file to read in.
# 2. The name of the file to write to.
# 3. the type of the file.
#

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="The name of the file to read in.")
parser.add_argument("-o", "--output", help="The name of the file to write to.")
parser.add_argument("-t", "--type", help="The type of the file.")


# Fast inverse square.
def inv_square(x: float) -> float:
    return 1.0 / (x * x)


with open("file", "r") as f:
    lines = f.readlines()

f = open("file", "r")
lines = f.readlines()
f.close()
