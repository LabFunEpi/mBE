import sys
def find_gene(names, ids, gene):
    if gene[0:4] == "ENSG":
        if gene in ids.keys():
            return(ids[gene])
    for key in names.keys():
        if gene == key:
            return(key)
        else:
            for value in names[key]:
                if gene in value:
                    return(key)
    return(None)

f = open("aliases")
lines = f.readlines()
f.close()

names = {}
for line in lines:
    words = line.strip().split("\t")
    names[words[0]] = [words[3].split(", "), words[2].split(", "), words[1].split(", ")]

f = open("gencode.v19.genes.id")
lines = f.readlines()
f.close()
 
ids = {}
for line in lines:
    words = line.strip().split("\t")
    ids[words[0]] = words[1]
    
f = open("GeneHancer_Bed.txt")
lines = f.readlines()
f.close()

foundnames = {}
for line in lines:
    words = line.strip().split("\t")
    gene = words[4]
    name = ""
    if gene in foundnames.keys():
        name = foundnames[gene]
    else:
        name = find_gene(names, ids, gene)
        foundnames[gene] = name
    if name == None:
        name = gene
    print("\t".join(words[0:4]) + "\t" + name + "\t" + words[5])
    sys.stdout.flush()
