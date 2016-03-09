# useful-scripts
Various useful scripts
--

### Search patterns in FASTQ files

Made by Vasily Cvetkov [@seaneuron](https://github.com/seaneuron).

#### Dependencies:
 
 - Python3

 - Biopython

 - tqdm

#### Usage

Help message:
```
usage: vsearch.py [-h] [-k K] [--input INPUT] [--start START] [--end END]
                  --path PATH [--output OUTPUT]

optional arguments:
  -h, --help            show this help message and exit
  -k K                  number of allowed errors for levenshtein distance (default=2)
  --input INPUT, -i INPUT
                        a pattern or a comma-separated list of patterns (default=AAAAAAAAAA)
  --start START, -s START
                        0-based starting position of the search area (default=0)
  --end END, -e END     0-based ending position of search area (default=30)
  --path PATH, -p PATH  path to the input file or a folder with input files
  --output OUTPUT, -o OUTPUT
                        path to the output files or folder
```

Example:

```
$python3 vsearch.py -k 3 -i ACGATTCA,ACACG,AAAA -s 10 -e 50
```