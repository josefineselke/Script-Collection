#!/usr/bin/env python3

# this script counted positive cells from an alignment.sam file
import sys
import glob
import csv
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from Bio import SeqIO
import time


def read_metrics(file):
    linecount_filtered = 0
    barcodes = {} # Barcode: ID
    total = 0
    fsize = Path(file).stat().st_size
    sample = file.split('/')[1]

    fastq = SeqIO.index(glob.glob(params[0] + sample + '*R1*fastq')[0], 'fastq')
    with open(file, 'r') as filtered:
        with tqdm(total=fsize, desc=file) as pbar:
            start = time.time()
            for line in filtered:
                ID = line.split('\t')[0]
                linecount_filtered += 1
                if str(fastq[ID].seq)[0:16] not in barcodes.keys():
                    barcodes[str(fastq[ID].seq)[0:16]] = ID
                total += len(line)
                pbar.update(total - pbar.n)
        end = time.time()
    tsv = str(params[1] + '/metrics.txt')
    log = str(params[1] + '/metrics.log')
    with open(tsv, 'a') as metrics:
        tsv_writer = csv.writer(metrics, delimiter='\t', lineterminator='\n')
        tsv_writer.writerow([file, str(linecount_filtered), str(len(barcodes))])
    with open(log, 'a') as log:
        log_writer = csv.writer(log, delimiter='\t', lineterminator='\n')
        log_writer.writerow([file, round(end-start, 2)])
        for record in barcodes:
            log_writer.writerow([record, barcodes[record]])

if __name__ == '__main__':
    print("START!")
    global params
    # params include path to fastq, directory with results, tool
    params = sys.argv[1:]
    files = glob.glob(params[2] + '/*.filtered.sam', recursive = True)

    # parallelization if possible
    with ThreadPoolExecutor() as executor:
        result = executor.map(read_metrics, files)
    # if not 
    for file in files:
        read_metrics(file)
