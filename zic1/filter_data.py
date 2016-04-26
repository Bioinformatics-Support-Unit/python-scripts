import sys

def main(f):
    zic1 = '206373_at'
    probesets = get_probesets()
    fh = open(f, 'r')
    for line in fh:
        if line.startswith("Scan"):
            print line.rstrip()
        if line.startswith(zic1):
            print line.rstrip()
        for probeset in probesets:
            if line.startswith(probeset):
                print line.rstrip()

def get_probesets():
    fh = open('probeset_list.txt', 'r')
    ps = []
    for line in fh:
        ps.append(line.rstrip())
    return ps

if __name__ == '__main__':
    main(sys.argv[1])
