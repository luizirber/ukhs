MAP = dict(zip("ACTG", "TGAC"))

def revc(seq):
    return "".join(MAP[c] for c in reversed(seq))

with open("data/res_7_20_4_0.txt", "r") as f:
    hashes = [l.strip() for l in f.readlines()]

query = "ACACCGTAGCCTCCAGATGC"
k = 7

matches = []
for i in range(len(query) - k + 1):
    if query[i: i+k] in hashes:
        matches.append(query[i:i+k])
print(*matches, sep="\n")

for match in matches:
    rc = revc(match)
    if rc in matches:
        print('kmer: {}, revc: {}'.format(match, rc))
