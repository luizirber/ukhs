with open("data/res_7_20_4_0.txt", "r") as f:
    hashes = [l.strip() for l in f.readlines()]

query = "ACACCGTAGCCTCCAGATGC"
k = 7

for i in range(len(query) - k + 1):
    if query[i: i+k] in hashes:
        print(query[i:i+k])
