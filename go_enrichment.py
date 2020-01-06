import sys
from scipy.stats import hypergeom
from collections import defaultdict

gene_go = {}
go_dict = {}

association_file = sys.argv[1]

try:
   fassoc = open(association_file, "r")
         
   for line in fassoc:
     line = line.rstrip()
     id = line.split("\t")
     go_id = id[1]
     num_fields = len(id)
     go_dict[go_id] = num_fields - 2

     for i in range(3, num_fields):
         gene_id = id[i]

         if gene_id not in gene_go:
            go_value = []
            go_value.append(go_id) 
            gene_go[gene_id] = go_value
         else:
            gene_go[gene_id].append(go_id)
   fassoc.close()

except IOError:
       print ('cannot open', association_file)
       fassoc.close()


N = len(gene_go.keys())

n = 0
gofreq = {}

gene_file = sys.argv[2]

try:
   fgene = open(gene_file, "r")
   for gline in fgene:
     gline = gline.rstrip()
     n += 1
     geneids = gline.split(",")

     if geneids[0] in gene_go:
        go_list = gene_go[geneids[0]]

        for go in go_list: 
            if go in gofreq:
               gofreq[go] += 1 
            else:
               gofreq[go] = 1


   for go_key, frequency in gofreq.items():
       k = frequency
       K = go_dict[go_key]
       prb = hypergeom.sf(k, N, n, K)

       print ("GO ID = " + go_key + "\tN = " + str(N) + "\tK = " + str(K) + "\tn = " + str(n) + "\tk = " + str(k) + "\tSignificance = " +str(prb))
   fgene.close()
    
except IOError:
        print ('cannot open', gene_file)
        fgene.close()

