#!/usr/bin/env python
import os, pickle
# Script for converting a pickle directory to a flatfile (!)
cwd_contents = os.listdir('pickles')
pickle_filenames = [os.path.join("pickles",f)
                    for f in cwd_contents if "unique" in f]
with open("tf_hits.csv","w") as f:
    for pf in pickle_filenames:
        print pf
        with open(pf) as pf_handle:
            d = pickle.load(pf_handle)
        print "loaded dictionary"
        for gene in d:
            for tf in d[gene]:
                for (pos, score) in d[gene][tf]:
                    f.write("%s, %s, %s, %s\n" % (gene, tf, pos, score))
