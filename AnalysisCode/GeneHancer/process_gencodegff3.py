

f = open("/references/gffs/gencode.v19.annotation.gff3")
lines = f.readlines()
f.close()

o1 = open("gencode.v19.genes.bed", "w")
o2 = open("gencode.v19.genes.id", "w")
for line in lines:
    if line[0] == "#":
        continue
    words = line.strip().split("\t")
    if words[2] == "gene":
        gene_id = ""
        gene_name = ""
        details = words[8].split(";")
        for detail in details:
            if detail[:8] == "gene_id=":
                gene_id = detail[8:]
            if detail[:10] == "gene_name=":
                gene_name = detail[10:]
        o1.write("\t".join([words[0], words[3], words[4], words[6], gene_name]) + "\n")
        o2.write("\t".join([gene_id[:15], gene_name]) + "\n")
o1.close()
o2.close()

