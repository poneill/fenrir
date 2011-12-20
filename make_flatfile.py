import os, pickle
cwd_contents = os.listdir('.')
pickle_filenames = filter(lambda f: "unique" in f, cwd_contents)
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
