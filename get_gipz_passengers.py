import argparse
import sys
import pandas as pd

def get_passengers(gipz, fasta):
    common5s = []
    passengers = []
    loops = []
    dd = pd.read_table(gipz)
    titles = dd['pSM2c Oligo ID']
    hairpins = dd['Hairpin Sequence']
    for hairpin in hairpins:
        common5s.append(hairpin[0:18])
        passengers.append(hairpin[18:40])
        loops.append(hairpin[40:58])
    for i in range(0, len(titles)):
        fasta.write('>'+titles[i]+'\n')
        fasta.write(passengers[i]+loops[i]+'\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Make FASTA sequence of passengers from GIPZ library spreadsheet')
    parser.add_argument('gipz_library', type=argparse.FileType('r'), help='GIPZ library file (tab-delimited text export of vendor provided Excel)')
    parser.add_argument('outfile', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help="Optional output file for fasta sequences (or write to STDOUT)")
    args = parser.parse_args()
    get_passengers(args.gipz_library, args.outfile)
