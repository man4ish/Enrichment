from scipy.stats import hypergeom
from collections import defaultdict

gene_go={}
go_count={}
go_dict = {}
f=open("167_go.gmt", "r")
for x in f:
  x=x.rstrip()
  id=x.split("\t")
  key = id[1]
  length = len(id)
  go_dict[key] = length-2;
  for i in range(3,length):
      go_count[id[i]]=1 
      gene_id = id[i]
      go_id = id[1]
      if gene_id not in gene_go:
         go_value = []
         go_value.append(go_id) 
         gene_go[gene_id] = go_value
         #print(id[i]+"-->"+key)
      else:
         gene_go[gene_id].append(go_id)
  #print(key)   
  #goannotation_dict[key].append(id[i])

N=len(go_count.keys())

#for k, v in gene_go.items():
#    print(k, v)

n = 0
gofreq = {}
fgene=open("genelist", "r")
for line in fgene:
  line=line.rstrip()
  n += 1
  geneids=line.split(",")
  if geneids[0] in gene_go:
     #print(gene_go[geneids[0]])
     go_list = gene_go[geneids[0]]
     for i in go_list: 
         if i in gofreq:
            gofreq[i] += 1 
         else:
            gofreq[i] = 1


for go_key, frequency in gofreq.items():
    #print(go_key, frequency)
    k = frequency
    K = go_dict[go_key]
    prb = hypergeom.sf(k, N, n, K)
    print ("GO ID = " + go_key + "\tN = " + str(N) + "\tK = " + str(K) + "\tn = " + str(n) + "\tk = " + str(k) + "\tSignificance = " +str(prb))
    


