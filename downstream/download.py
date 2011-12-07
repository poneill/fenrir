import os, pickle,sys
from Bio import Entrez, SeqIO
sys.path.append("..")
from read_matrix import *
Entrez.email = "pon2@umbc.edu"
print os.getcwd()
print "parsing upstream regions"
if not os.path.isfile("upstream_regions.pkl"):
    upstream_regions = parse_urs("../upstream5000.fa")
    with open("upstream_regions.pkl",'w') as ur_handle:
        pickle.dump(upstream_regions,ur_handle)
else:
    with open("upstream_regions.pkl") as ur_handle:
        upstream_regions = pickle.load(ur_handle)
print "computing mapping"
mapping = open("gene_mrna_mapping.txt").readlines()
print "computing refseqs"
refseqs = [filter(lambda x: x,line.strip().split('\t')) for line in mapping if line.strip().startswith('NM')]
print "computing accession numbers"
nms = [nm[:nm.index('_',3)] for nm in upstream_regions.keys()]
print "finding relevant refseqs"
if not os.path.isfile("desired_refseqs.pkl"):
    desired_refseqs = [rs for rs in refseqs if rs[0] in nms]
    with open("desired_refseqs.pkl",'w') as pickle_handle:# heh
        pickle.dump(desired_refseqs,pickle_handle)
else:
    with open("desired_refseqs.pkl") as pickle_handle:
        desired_refseqs = pickle.load(pickle_handle)
wrong_lengths = []
bad_apples = []
retries = []
genes = {}
print "starting lookup"
for dr in desired_refseqs:
    if len(dr) == 3:
        idnum = [entry for entry in dr
                 if (entry != '1' and entry != '-1' and not entry.startswith('NM_'))][0]
        print "looking up ", idnum
        handle = Entrez.esearch(db="nucleotide", term=idnum)
        record = Entrez.read(handle)
        print "found", len(record['IdList']), "entries"
        genes[idnum] = {}
        for il in record['IdList']:
            if os.path.isfile(il + ".gbk"):
                print "already looked up", il
                with open(il + ".gbk") as in_handle:
                    try:
                        record = SeqIO.read(in_handle, "genbank")
                        genes[idnum][il] = record
                    except:
                        print "record for",il, "was empty...bailing!"
                        retries.append(il)
            else:
                print "fetching", il
                net_handle = Entrez.efetch(db="nucleotide",id=il,rettype="gb")
                try:
                    record = SeqIO.read(net_handle,"genbank")
                except:
                    print "couldn't read",il, "from net_handle...bailing!"
                    bad_apples.append(il)
                    continue
                print "writing to disk"
                try:
                    with open(il + ".gbk",'w') as out_handle:
                        out_handle.write(record.format("gb"))
                    print "storing in dictionary"
                    genes[idnum][il] = record
                except:
                    print "couldn't write", il, "to disk... bailing!"
                    bad_apples.append(record)
                    continue
    else:
        wrong_lengths.append(dr)

print "Summary"
print "_" * 80
print "scanned", len(desired_refseqs), "refseqs"
print "recovered", len([il for gene in genes for il in gene]), "genes"
print len(wrong_lengths), "were of wrong length"
print "encountered problems with", len(bad_apples)
print "should retry", len(retries)
