"""Importing packages"""
from Bio import SeqIO
from Bio import pairwise2
import argparse
import os, sys, tarfile, gzip
import time
from datetime import date
from time import sleep
from tqdm import tqdm

"""Number of samples function"""
def numos(record, sample, k, start, end):
    result = ''
    a = pairwise2.align.globalms(record[start:end], sample, 1, -1, -1, -1, penalize_end_gaps = False , score_only = True, one_alignment_only = True)
    if a:
        if a >= len(sample)-k:
            result = a
    return result

"""Fasta and fastq handle"""
def fastq_handle(file, destpath, samples, k, start, end, i, b):
    i = 0
    a = {}
    td = today_is()
    flag = False
    header = ""
    c = 0
    outfiles={}
    flag={}
    for sample in samples:
        outfiles[sample] = open(destpath+'_'+td+'_'+'0'+str(i)+'_'+sample+'.fastq', "w")
        print('Processing sample {0}'.format(sample))
        flag[sample]=False
        a[sample]=0
    for num_line in tqdm(range(b)):
        line = file.readline()
        i+=1
        if i%4 == 1:
            header = line
        elif i%4 == 2:
            for sample in samples:
                res = numos(line, sample, k, start, end)
                if res:
                    outfiles[sample].write(header)
                    outfiles[sample].write(line)
                    outfiles[sample].write('+\n')
                    flag[sample] = True
                    a[sample]+=1
        elif i%4 == 0:
            for sample in samples:
                if flag[sample]==True:
                    outfiles[sample].write(line)
                flag[sample] = False
    for sample in samples:
        outfiles[sample].close()
    return a

def fastq_gz_handle(file, destpath, samples, k, start, end, i, b):
    i = 0
    a = {}
    td = today_is()
    flag = False
    header = ""
    c = 0
    outfiles={}
    flag={}
    for sample in samples:
        outfiles[sample] = open(destpath+'_'+td+'_'+'0'+str(i)+'_'+sample+'.fastq', "w")
        print('Processing sample {0}'.format(sample))
        a[sample]=0
        flag[sample]=False
    for num_line in tqdm(range(b)):
        line = file.readline()
        i+=1
        if i%4 == 1:
            header = line
        elif i%4 == 2:
            for sample in samples:
                res = numos(line, sample, k, start, end)
                if res:
                    outfiles[sample].write(header.decode())
                    outfiles[sample].write(line.decode())
                    outfiles[sample].write('+\n')
                    flag[sample] = True
                    a[sample]+=1
        elif i%4 == 0:
            for sample in samples:
                if flag[sample]==True:
                    outfiles[sample].write(line.decode())
                flag[sample] = False
    for sample in samples:
        outfiles[sample].close()
    return a

def fasta_handle(file, destpath, sample, k, start, end, i):
    a = 0
    td = today_is()
    with open(destpath+'_'+td+'_'+'0'+str(i)+'_'+sample+'.fasta', "w") as output:
        for record in SeqIO.parse(file, "fasta"):
            text = str(record.seq)
            result = numos(text, sample, k, start, end)
            if result:
                record.id += "_" + result
                SeqIO.write(record, output, 'fasta')
                a += 1
    return a

def file_type(file):
    line = file.readline()
    ftype = ""
    if line.startswith('>'):
        ftype = "a"
    elif line.startswith('@'):
        ftype = "q"
    else:
        print ("This file is neither fasta nor fastq. Consider another file for input") 
    return ftype

def archived_file_type(file):
    line = file.readline()
    ftype = ""
    if line.startswith(b'>'):
        ftype = "a"
    elif line.startswith(b'@'):
        ftype = "q"
    else:
        print ("This file is neither fasta nor fastq. Consider another file for input") 
    return ftype

def today_is():
    today = date.fromtimestamp(time.time())
    t = today.timetuple()
    td = str(t[0])+'-'+str(t[1])+'-'+str(t[2])
    return td

"""Open file"""
def process_file(path, destpath, samples, k, start, end, i):
    count = 0
    b = 0
    name = os.path.basename(path)
    outpath = os.path.join(destpath, name)
    if path.endswith('.gz'):
        with gzip.open(path) as file:
            ftype = archived_file_type(file)
        with gzip.open(path) as file:
            for line in file:
                b+=1
        if ftype == "a":
            with gzip.open(path) as file:
                count = fasta_handle(file, outpath, samples, k, start, end, i)
                for key, value in count.items():
                    print("The number of sequences {0} in file is {1}".format(key, value))
        elif ftype == "q":
            with gzip.open(path) as file:
                count = fastq_gz_handle(file, outpath, samples, k, start, end, i, b)
                for key, value in count.items():
                    print("The number of sequences {0} in file is {1}".format(key, value))
    else:
        with open(path) as file:
            ftype = file_type(file)
        with open(path) as file:
            for line in file:
                b+=1
        if ftype == "a":
            with open(path) as file:
                count = fasta_handle(file, outpath, samples, k, start, end, i)
                for key, value in count.items():
                    print("The number of sequences {0} in file is {1}".format(key, value))
        elif ftype == "q":
            with open(path) as file:
                count = fastq_handle(file, outpath, samples, k, start, end, i, b)
                for key, value in count.items():
                    print("The number of sequences {0} in file is {1}".format(key, value))
    return count

if __name__ == "__main__":
    "Input"
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', default = 2, type=int, help='number of allowed errors for levenshtein distance (default=2)')
    parser.add_argument('--input','-i', default = 'AAAAAAAA', help='a pattern or a comma-separated list of patterns (default=AAAAAAAAAA)')
    parser.add_argument('--start', '-s', default = 0, type=int, help='0-based starting position of the search area (default=0)')
    parser.add_argument('--end', '-e', default = 30, type=int, help='0-based ending position of search area (default=30)')
    parser.add_argument('--path', '-p', required=True, help='path to the input file or a folder with input files')
    parser.add_argument('--output', '-o', default = './', help='path to the output files or folder')
    args = parser.parse_args()


    k = args.k
    samples = args.input.split (',')
    start = args.start
    end = args.end
    path = args.path
    destpath = args.output
    files = []
    i = 0

    if os.path.isfile(path):
        files.append(path)
        print(files)
        if destpath:
            pass
        else:
            destpath = os.path.dirname(path)
    elif os.path.isdir(path):
        if destpath:
            pass
        else:
            destpath = path + '/'
        for entry in os.listdir(path):
            entry = os.path.join(path, entry)
            if os.path.isfile(entry):
                files.append(entry)
            else:
                print(entry)
                print('File not found')

        
    for file in files:    
        print("File {0} is {1}".format(i+1, file))
        t1 = time.clock()
        count = process_file(file, destpath, samples, k, start, end, i) 
        t2 = time.clock()
        print (t2-t1)
        i += 1
            






