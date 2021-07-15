#!/usr/bin/env python

from argparse import ArgumentParser
import os
from Bio import SeqIO
import pandas as pd
import numpy as np
from joblib import load

import time

def string_two_mers():
    """Returns list of strings of all combinations of possible letters. """
    LETTERS = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]
    tuples = [(letter1, letter2) for letter1 in LETTERS for letter2 in LETTERS]
    return [''.join(tupl) for tupl in tuples]

def read_fastas_from_file(file):
    """Returns list of fasta records from a file. """
    fastas = []
    try:
        for record in SeqIO.parse(file, 'fasta'):
            fastas.append(record)
        print(f"Number of sequences in '{file}': {len(fastas)}")
        return fastas
    except:
        print(f"Failed to load file '{file}'.")
        return []

def read_fastas_from_directory(directory):
    """Returns list of fasta records from a directory. """
    print(f"Reading from directory '{directory}'.")
    fastas = []
    for filename in os.listdir(directory):
        fastas += read_fastas_from_file(os.path.join(directory, filename))
    return fastas

def count_two_mers(sequence, two_mers):
    """Returns dictionary, keys are possible couplesl of letters, 
    values are number of their occurances in a given sequence. """
    dictionary = {}
    for mer in two_mers:
        dictionary[mer] = sequence.count(mer)
    return dictionary

def construct_dataframe(two_mers, fastas):
    """Returns dataframe: rows are fasta records, 
    cols are 2mers in alphabetical order."""
    dictionary = {}
    print("Constructing 2mer dataframe:")
    for fasta in fastas:
        two_mers_count = count_two_mers(fasta.seq, two_mers)
        dictionary[fasta.id] = two_mers_count
    df = pd.DataFrame.from_dict(dictionary, orient='index')
    print("Dataframe constructed!")
    return df

def make_prediction(df):
    """Returns np array of predictions. """
    print("Predicting...")
    model_source = os.path.dirname(os.path.abspath(__file__)) + "/sr_a.joblib"
    model = load(model_source)
    prediction = model.predict(df)
    print("Finished prediction!")
    return prediction

def get_fastas(source):
    """Loads fasta records either from file or dictionary. """
    fastas = []
    if os.path.isfile(source):
        fastas = read_fastas_from_file(source)
    else:
        fastas = read_fastas_from_directory(source)
    return fastas
    

def start(source, output):
    """Starts the process... """
    fastas = get_fastas(source)
    if(len(fastas) == 0):
        print("No fasta records found.")
        return
    else:
        print(f"Found {len(fastas)} fasta records.")
    
    two_mers = string_two_mers()
    
    two_mer_df = construct_dataframe(two_mers, fastas)
    print("DF shape:", two_mer_df.shape)
    print("Head of DF:")
    print(two_mer_df.head(10))
    
    prediction = make_prediction(two_mer_df)
    
    print(prediction)
    print(type(prediction), len(prediction))
    print(f"Number of positively predicted: {np.count_nonzero(prediction)}")
    
    two_mer_df['prediction'] = prediction
    result = pd.DataFrame(two_mer_df['prediction'])
    result.to_csv(output, index=True, index_label="Protein ID")
    print(f"Predictions saved to {output}.")
    
def now():
    return time.time()

def main():
    """Takes the arguments from parser and calls
    the start function."""
    
    parser = ArgumentParser()
    parser.add_argument('-f', "--fasta",  help='path to fasta file or directory of fasta files', required=True)
    parser.add_argument('-o', "--output", help='output file', required=True)
    args = parser.parse_args()
    fasta_file = args.fasta
    output_file = args.output

    start_time = now()
    print("Welcome to phibers.")
    start(fasta_file, output_file)
    end_time = now()
    print(f"Finished! It took {end_time - start_time} seconds.")

if __name__ == '__main__':
    main()
