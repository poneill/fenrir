#!/bin/python
import os, pickle,sys, time
from Bio import Entrez, SeqIO
sys.path.append("..")
from read_matrix import *

def log(text):
    print text
    with open("download.log","a") as logfile:
        logfile.write(text + "\n")

Entrez.email = "pon2@umbc.edu"
log( os.getcwd())
log("starting at" + time.asctime(time.localtime()))
log( "parsing upstream regions")
if not os.path.isfile("upstream_regions.pkl"):
    upstream_regions = parse_urs("../upstream/upstream5000.fa")
    with open("../upstream/upstream_regions.pkl",'w') as ur_handle:
        pickle.dump(upstream_regions,ur_handle)
else:
    with open("../upstream/upstream_regions.pkl") as ur_handle:
        upstream_regions = pickle.load(ur_handle)
log("computing mapping")
mapping = open("gene_mrna_mapping.txt").readlines()
log("computing refseqs")
refseqs = [filter(lambda x: x,line.strip().split('\t'))
           for line in mapping if line.strip().startswith('NM')]
log("computing accession numbers")
nms = [nm[:nm.index('_',3)] for nm in upstream_regions.keys()]
log("finding relevant refseqs")
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
log("starting lookup")
for dr in desired_refseqs:
    if len(dr) == 3:
        idnum = [entry for entry in dr
                 if (entry != '1' and entry != '-1'
                     and not entry.startswith('NM_'))][0]
        log("looking up ", idnum)
        handle = Entrez.esearch(db="nucleotide", term=idnum)
        record = Entrez.read(handle)
        log("found", len(record['IdList']), "entries")
        genes[idnum] = {}
        for il in record['IdList']:
            if os.path.isfile(il + ".gbk"):
                log("already looked up", il)
                with open(il + ".gbk") as in_handle:
                    try:
                        record = SeqIO.read(in_handle, "genbank")
                        genes[idnum][il] = record
                    except:
                        log("record for",il, "was empty...bailing!")
                        retries.append(il)
            else:
                log("fetching", il)
                net_handle = Entrez.efetch(db="nucleotide",id=il,rettype="gb")
                try:
                    record = SeqIO.read(net_handle,"genbank")
                except:
                    log("couldn't read",il, "from net_handle...bailing!")
                    bad_apples.append(il)
                    continue
                log("writing to disk")
                try:
                    with open(il + ".gbk",'w') as out_handle:
                        out_handle.write(record.format("gb"))
                    log("storing in dictionary")
                    genes[idnum][il] = record
                except:
                    log("couldn't write", il, "to disk... bailing!")
                    bad_apples.append(record)
                    continue
    else:
        wrong_lengths.append(dr)

log("Summary")
log("_" * 80)
log("scanned", len(desired_refseqs), "refseqs")
log("recovered", len([il for gene in genes for il in gene]), "genes")
log(len(wrong_lengths), "were of wrong length")
log("encountered problems with", len(bad_apples))
log("should retry", len(retries))

if not os.path.isfile("genes.pkl"):
    with open("genes.pkl",'w') as gene_handle:# heh
        pickle.dump(genes,gene_handle)
else:
    log("warning: genes.pkl already exists")
    log("saving to backup")
    with open("genes.pkl.bak",'w') as gene_handle:# heh
        pickle.dump(genes,gene_handle)
