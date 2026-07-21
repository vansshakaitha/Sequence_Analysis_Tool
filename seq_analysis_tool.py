from Bio import SeqIO, Entrez
from Bio.SeqUtils import gc_fraction
import pandas as pd

def find_orf(seq):
    orf_seqs=[]
    for frame in range(3):
        for i in range(frame,len(seq)-2,3):
            codon=seq[i:i+3]
            if codon=="ATG":
                j=i
                while j<len(seq)-2:
                    stop_codon=seq[j:j+3]
                    if stop_codon in ["TAA","TAG","TGA"]:
                        length=j+3-i
                        if length>100:
                            orf_seqs.append(seq[i:j+3]) 
                        break    
                    j+=3
    return orf_seqs

Entrez.email="vanssha@gmail.com"
handle=Entrez.esearch(db="nuccore", term="TP53 [Gene] AND Homo sapiens[Organism] AND CDS[Feature Key]", retmax=10)
read_rec=Entrez.read(handle)
handle.close()

id_list=read_rec["IdList"]
fetch_handle=Entrez.efetch(
    db="nuccore",
    id=id_list,
    rettype="fasta",
    retmode="text"
)
recs=list(SeqIO.parse(fetch_handle,"fasta"))
fetch_handle.close()

data=[]
for record in recs:
    orfs=find_orf(record.seq)
    if orfs:
        longest_orf=len(max(orfs,key=len))
    else:
        longest_orf=0
    data.append({
        "ID"       : record.id,
        "Length"   : len(record.seq),
        "GC %"     : round(gc_fraction(record.seq)*100,2),
        "ORF count": len(orfs),
        "Longest"  : longest_orf
    })

df=pd.DataFrame(data)

df.to_csv("tp53_cds_orf-results.csv",index=False)


    


    
