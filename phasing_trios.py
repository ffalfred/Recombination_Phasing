import math
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import operator
import phasing_class
import argparse

if __name__ == '__main__':
#Arguments
    parser = argparse.ArgumentParser(description='Simple tool for phasing human PB-oocyte trios.')
    parser.add_argument("-i","--input_file", default=False, required=True, help="Input file in .txt format")
    parser.add_argument("-m","--min_snps", default=20, required=True, type=float,help="Amount of SNPs required to define a phase")
    parser.add_argument("-o","--output_dir", default=False, required=True,help="Directory where the files are saved")
    args = parser.parse_args()

    end=phasing_class.phasing_class_tool(args.input_file,args.min_snps,args.output_dir)
    end.phasing() 

