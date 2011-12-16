#!/usr/bin/python
print "starting"
import sys, re, math, random, itertools, time, pickle
import matrix_parser, analysis
default_sequence_file = "upstream/upstream5000.fa"
base_pair_ordering = "acgt"

def lexicographic_cmp(xs,ys):
    if not xs or not ys:
        return cmp(ys,xs)
    elif xs[0] != ys[0]:
        return cmp(xs,ys)
    else:
        return lexicographic_cmp(xs[1:],ys[1:])

def pick_winner(pssm):
    """Given a pssm (or count matrix), return the maximal scoring
    string"""
    indices = [col.index(max(col)) for col in pssm]
    return "".join([base_pair_ordering[index] for (col,index) in zip(pssm,indices)])

def lexicographic_cmp2(xs,ys):
    i = 0
    for i, (x,y) in enumerate(itertools.izip(xs,ys)):
        if x != y:
            return cmp(xs[i:],ys[i:])
    return cmp(ys[i:],xs[i:])

def suffixes_of(word):
    suffixes = sorted([word[i:] for i in range(len(word))],
                      cmp=lexicographic_cmp2)
    return [suf + "$" for suf in (suffixes + [''])]

def pprint(x):
    for row in x:
        print row
        
def nmers(n):
    if n == 1:
        return ["A","C","G","T"]
    else:
        return sum([map(lambda(b):b+c,nmers(n-1)) for c in base_pair_ordering],[])
    
def safe_log2(x):
    """Implements log2, but defines log2(0) = 0"""
    return math.log(x,2) if x > 0 else 0

def complement(base):
    return {"A":"T","T":"A","G":"C","C":"G"}[base]
    
def wc(word):
    return map(complement, word[::-1])

def split_on(xs, pred):
    """Split xs into a list of lists each beginning with the next x
    satisfying pred, except possibly the first"""
    indices = [i for (i,v) in enumerate(xs) if pred(v)]
    return [xs[i:j] for (i,j) in zip([0]+indices,indices+[len(xs)]) if i != j]
    
def matches_accession_number(line):
    """Return an re.match object for the accession number pattern """
    return re.search(r'^AC\s+([A-Z0-9]+)', line)

def matches_column(line):
    """Return an re.match object for the column pattern"""
    regexp = """
^[0-9]+     #Begins with column number
\s+         #Followed by some whitespace
            #and now numbers corresponding to the base count at that column:
([.0-9]+)\s+ #number of As
([.0-9]+)\s+ #Cs
([.0-9]+)\s+ #Gs
([.0-9]+)    #Ts
"""
    return re.search(regexp,line,re.VERBOSE)

def matrix_from_lines(lines):
    """Convert raw column lines into matrix.  Assumes all lines are
    column lines"""
    return [map(float,matches_column(line).groups()) for line in lines]

def parse_lines_for_matrices(lines):
    accession_chunks = split_on(lines,matches_accession_number)[1:]
    #list containing list of lines beginning with given accession
    #number.  First chunk is filler.
    count_chunks = [[line for line in chunk
                     if matches_accession_number(line)
                     or matches_column(line)]
                    for chunk in accession_chunks]
    matrices = {}
    for count_chunk in count_chunks:
        accession_line = count_chunk[0]
        accession_name = matches_accession_number(accession_line).groups(0)[0]
        matrix = matrix_from_lines(count_chunk[1:])
        matrices[accession_name] = matrix
    return matrices
    
def normalize(xs):
    return map(lambda(x): x/float(sum(xs)),xs)

def pssm_from_matrix(matrix,background_probs = (0.25,)*4):
    """Accept count matrix (as nested list) and return pssm (as nested list)"""
    def convert_column(col):
        return [safe_log2(c/p) for (c,p) in zip(normalize(col),background_probs)]
    return [convert_column(col) for col in matrix]

def score(pssm,word):
    return sum([col[base_pair_ordering.index(base)] for col, base in zip(pssm,word)])

def tail(xs):
    if xs:
        return xs[1:]
    else:
        return []
    
def partial_sums(xs,acc = 0):
    return [sum(xs[:i+1]) for i in range(len(xs))]

def partial_thresholds(pssm,theta):
    """Returns a list of partial thresholds to be interpreted as
    follows: After having read position i, you must have scored at
    least pt[i], hence pt[len(pssm)] >= theta if the word is a
    positive"""
    rev_maxes = map(max,pssm[::-1])
    rev_thetas = reduce(lambda ths,x: ths + [ths[-1] - x],rev_maxes,[theta])
    return rev_thetas[::-1][1:]
    
def pairs(xs):
    return zip([None] + xs, xs + [None])[1:-1]

def lcp_functional(x,y):
    if not x or not y or x[0] != y[0]:
        return 0
    else:
        return 1 + lcp(x[1:],y[1:])

def lcp(x,y):
    i = 0
    length = min(len(x),len(y))
    for i, (x,y) in enumerate(zip(x,y)):
        if x != y:
            return i
    return length
    
class ESA(object):
    def __init__(self,word):
        self.word = word
        self.suffixes = suffixes_of(word)
        self.suf = [len(word) - len(self.suffixes[i]) + 1
                    for i in range(len(word)+1)]
        self.lcp = [0] + [lcp(x,y) for (x,y) in pairs(self.suffixes)]
        self.skp = [min([j for j in range(i+1, len(word) + 1)
                         if self.lcp[j] < self.lcp[i]] + [len(word) + 1])
                         for i in range(len(word) + 1)]

test_pssm = [[1,3],[3,2]]

def algorithm1(esa, pssm, theta):
    matches = []
    C = {}
    thetas = partial_thresholds(pssm, theta)
    suf, lcp, skp, suffixes = esa.suf, esa.lcp, esa.skp, esa.suffixes
    depth = 0
    i = 0
    n = len(esa.word)
    if 'n' in esa.word:
        return [(-1,-1)] #if there exist ns in the string, deal with
                         #it later
    m = len(pssm)
    M = lambda d,char: pssm[d][base_pair_ordering.index(char)]
    while (i < n):
        if n - m < suf[i]: #if too far in to match
            while(n - m < suf[i] and (i < n)):
                i += 1
                depth = min(depth,lcp[i])
            if i >= n:
                return matches
        if depth == 0:
            score = 0
        else:
            score = C[depth - 1]
        d = depth - 1
        sentinel = True #hacks around the do-while loop
        while(sentinel or (d < m -1 and score >= thetas[d])):
            sentinel = False
            d = d + 1
            score = score + M(d,suffixes[i][d])
            C[d] = score
        if(d == m - 1 and score >= theta): 
            matches.append((suf[i],score))
            while(i < n):
                i += 1
                if lcp[i] >= m:
                    matches.append((suf[i],score))
                else:
                    break
        else:
            i = skipchain(lcp,skp,n,i,d)
        depth = lcp[i]
    return matches
    
def skipchain(lcp,skp,n,i,d):
    j = i + 1
    if i < n:
        while((j <= n) and (lcp[j] > d)):
            j = skp[j]
    else:
        j = n
    return j

def skipchain_old(lcp,skp,n,i,d):
    print "calling skipchain with n = {0},i  = {1}, d = {2}".format(n,i,d)
    if i < n:
        j = i + 1
        while((j <= n) and (lcp[j] > d)):
            j = skp[j] + 1
    else:
        j = n
    return j

def naive_match(esa,pssm,theta):
    word = esa.word
    window = len(pssm)
    scores = [score(pssm,word[i:i+window]) for i in range(len(word))]
    matches = [(i,s) for (i,s) in enumerate(scores) if s >= theta]
    return matches

def list_all_scores(pssm):
    if not pssm:
        return [0]
    else:
        return [w + x for w in pssm[0] for x in list_all_scores(pssm[1:])]

def sample_scores(pssm,n):
    return [sum(random.choice(w) for w in pssm) for i in range(n)]

cutoffs = {}
def return_cutoff(pssm, alpha, n):
    "Return an estimated cutoff value theta such that P(score(w,pssm) > theta) < alpha"
    if (str(pssm),alpha,n) in cutoffs:
        return cutoffs[(str(pssm),alpha,n)]
    else:
        samples = sorted(sample_scores(pssm,n),reverse=True)
        cutoffs[(str(pssm),alpha,n)] = samples[int(alpha*n)]
        return cutoffs[(str(pssm),alpha,n)]

def return_max(pssm, n):
    "Return sample max"
    samples = sorted(sample_scores(pssm,n),reverse=True)
    return samples[0]

def parse_urs(filename):
    "Parse filename (e.g. upstream5000.fa) and return a list of upstream regions"
    urs = {}
    current = ""
    next_name = ""
    for i, line in enumerate(open(filename)):
        match = re.match("[atgcn]+",line)
        if match:
            current += match.group()
        else:
            name_match = re.match(">(.*)",line)
            if name_match:
                current_name = next_name
                next_name = name_match.groups(0)[0]
                if current:
                    urs[current_name] = current
                    current = ""
            else:
                print "anomaly on line: ",i
    urs[next_name] = current
    return urs
            
def search_urs_for_pssms(urs,tfs,alpha,n=None):
    """Return a nested dictionary containing, for every / the first n
    ur and every tf, a list of locations and scores at which each tf
    binds"""
    print "starting at " + time.asctime()
    print "searching %s promoter regions for %s transcription factors" % (len(urs), len(tfs))
    if not n:
        n = len(urs)
    sites = {}
    for i, ur in enumerate(urs.keys()[:n]):
        print i, ur
        sites[ur] = {}
        print "Constructing ESA..."
        esa = ESA(urs[ur])
        print "scanning"
        for tf in tfs:
            cutoff = return_cutoff(tf.pssm,alpha,10000)
            sites[ur][tf.ID] = algorithm1(esa, tf.pssm, cutoff)
    print "finishing at " + time.asctime()
    return sites

def random_word(n):
    return "".join([random.choice("ATGC") for i in range(n)])

def test_lcps(n):
    for i in range(n):
        a = random_word(random.randrange(1000))
    	b = random_word(random.randrange(1000))
    	if lcp2(a,b) == lcp(a,b):
    		print True
    	else:
    		print False, a, b
print __name__
if __name__ == "__main__":
    gene_file = sys.argv[1] if len(sys.argv) > 1 else default_sequence_file
    print "making transfac table"
    tt = matrix_parser.TransfacTable("dat/matrix.dat")
    toi = [line.strip() for line in open("terms_of_interest.txt").readlines()]
    more_refined = [tf for term in toi
                    for tf in tt.entries
                    if (term.upper() in tf.ID.upper()
                        or term.upper() in tf.NA.upper())]
    print "parsing upstream regions"
    upstream_regions = parse_urs(gene_file)
    #eliminate redundancies:
    uniques = set(upstream_regions.values())
    urs = {}
    for ur in upstream_regions:
        if upstream_regions[ur] in uniques:
            urs[ur] = upstream_regions[ur]
    print len(urs)
    results = search_urs_for_pssms(urs,more_refined,.001)
    with open(filename + ".pickle",'w') as f:
        pickle.dump(results,f)
    print "done"
