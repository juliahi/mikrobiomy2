import sys

with open(sys.argv[1]) as infile:
    with open(sys.argv[2], 'w+') as outfile:
        counter = 0
        for line in infile:
            if line.startswith('>'):
                outfile.write(">contig-%d\n"%counter)
                counter += 1
            else:
                outfile.write(line)


