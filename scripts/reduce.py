import re

ROWS = 10

with open('../data/hyperleda.tsv', 'r') as fr:
    with open('../data/hyperleda_reduced.tsv', 'w') as fw:

        comments = True
        rows = 0

        for line in fr.readlines():
            if line[0] != '#':
                if comments:
                    line = re.sub(r' ', r'\t', line)
                    fw.write(line)

                    comments = False
                else:
                    fw.write(line)

                    rows += 1

            if ROWS is not None and rows >= ROWS:
                break
