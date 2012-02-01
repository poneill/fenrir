from __future__ import division
from scipy.stats import gaussian_kde
import numpy as np
from matplotlib import pyplot as plt
import random, os, pickle

def score_distribution(tf, score_dictionary):
    """given a dictionary of scores such that
    score_dictioanry[promoter][tf] returns a list of (position, score)
    pairs, return a list of scores for transcription factor tf"""
    return [score for promoter in score_dictionary
            for (position, score) in score_dictionary[promoter][tf]]


def total_score_distribution(score_dictionary):
    """given a dictionary of scores such that
    score_dictioanry[promoter][tf] returns a list of (position, score)
    pairs, return all scores"""
    return [score for promoter in score_dictionary
            for tf in score_dictionary[promoter]
            for (position, score) in score_dictionary[promoter][tf]]

def total_index_distribution(score_dictionary,cutoff=0):
    """given a dictionary of scores such that
    score_dictioanry[promoter][tf] returns a list of (position, score)
    pairs, return all scores"""
    return [position for promoter in score_dictionary
            for tf in score_dictionary[promoter]
            for (position, score) in score_dictionary[promoter][tf]
            if score > cutoff]

def index_distribution(promoter, score_dictionary, cutoff = 0,
                       entropy=True, bp_min=0,bp_max=5000):
    """Search a score_dictionary for a promoter and return the
    positions of hits falling above the cutoff (PSSM entropy by default)"""
    return [position for tf in score_dictionary[promoter]
            for (position, score) in score_dictionary[promoter][tf]
            if ((bp_min <= position <= bp_max) and
                 score > (tf_from_id(tf).entropy if entropy else cutoff))]

def tf_distribution(tf, score_dictionary, cutoff = 0,
                    entropy=True, bp_min=0,bp_max=5000):
    """Find the positional distribution of hits for a given TF over
    all promoters, i.e. where does the TF bind?"""
    return [position for promoter in score_dictionary
            for (position, score) in score_dictionary[promoter][tf]
            if ((bp_min <= position <= bp_max) and
                score > (tf_from_id(tf).entropy if entropy else cutoff))]

def tf_from_id(tf):
    return [mr for mr in more_refined if mr.ID == tf][0]

def plot_density(xs):
    density = gaussian_kde(map(float,xs))
    xs = np.arange(min(xs),max(xs),1)
    plt.plot(xs,density(xs))

def find_peaks(xs):
    return [i for i in range(len(xs)) if xs[i-1] < xs[i] > xs [i+1]]

def bin_indices(indices,n):
    """Partition indices into n bins"""
    bins = [[index for index in indices if i * (5000/n) <= index < (i + 1)*(5000/n)]  for i in range(n)]
    return bins

def len_bin_indices(indices,n):
    return map(len,bin_indices(indices,n))

def centroid(xs):
    return map(lambda coord: sum(coord)/len(coord), zip(*xs))

def distance(xs,centroid):
    return sum(map(lambda (x,c): (x - c)**2, zip(xs,centroid)))

def index_of_closest_centroid(xs,centroids):
    return [i for (i,d) in sorted([(i,distance(xs,centroid)) for i, centroid in enumerate(centroids)],
                key = lambda (i,d): d)][0]
    
def kmeans(xss, k):
    r = 0
    initial_assignment = [random.randint(0,k-1) for i in range(len(xss))]
    clusters = [[xs for i, xs in enumerate(xss) if initial_assignment[i] == j] for j in range(k)]
    centroids = map(centroid,clusters)
    updated_assignment = [index_of_closest_centroid(xs, centroids) for xs in xss]
    while initial_assignment != updated_assignment:
        r += 1
        initial_assignment = updated_assignment
        clusters = [[xs for i, xs in enumerate(xss) if initial_assignment[i] == j] for j in range(k)]
        centroids = map(centroid,clusters)
        updated_assignment = [index_of_closest_centroid(xs, centroids) for xs in xss]
        for i,cluster in enumerate(clusters):
            cent = centroids[i]
            distances = map(lambda x: distance(x,cent), cluster)
    for i,cluster in enumerate(clusters):
        cent = centroids[i]
        distances = map(lambda x: distance(x,cent), cluster)
    return assign(xss,updated_assignment)

def cluster_fitness(cluster):
    cent = centroid(cluster)
    return sum(map(lambda x: distance(x,cent)**2, cluster))
        
def fitness(clusters):
    return sum(map(cluster_fitness,clusters))
    
def assign(data,assignment):
    """group data by assignment"""
    return [[datum for i, datum in enumerate(data)
             if assignment[i] == j]
            for j in range(min(assignment),max(assignment) + 1)]

def select_from_population(pop,tournament_factor,xss):
    contender_assignments = [random.choice(pop) for i in range(tournament_factor)]
    contender_clusters = [assign(xss,ca) for ca in contender_assignments]
    (winner, fitness_winner) = sorted([(ca,fitness(assign(xss,ca))) for ca in contender_assignments],
                    key = lambda (ca,cf) : cf)[0]
    return winner

def breed(p,q):
    i = random.randint(0,len(p))
    r = p[:i] + q[i:]
    return r
    
def update_population(pop, tournament_factor, xss):
    popsize = len(pop)
    newpop = [breed(p,q) for (p,q) in [(select_from_population(pop,tournament_factor,xss),
                                        select_from_population(pop,tournament_factor,xss))
                                       for x in range(popsize)]]
#    print "size newpop: ", len(newpop)
    return newpop


def random_assignment(l,k):
    return [(lambda : random.randint(0,k))() for i in range(l)]

def ga(xss, k, popsize, generations, tournament_factor):
    L = len(xss)
    assignments = [random_assignment(L,k) for p in range(popsize)]
    for gen in range(generations):
        print gen
        assignments = update_population(assignments, k, xss)
    best_assignment = select_from_population(assignments, popsize, xss)
    return assign(xss,best_assignment)

def kmeans_iterate(xss,k,n):
    total = 0
    i = 0
    old_clusters = kmeans(xss,k)
    old_fit = fitness(old_clusters)
    while i < n:
        i += 1
        total += 1
        print i, total
        new_clusters = kmeans(xss,k)
        new_fit = fitness(new_clusters)
        if(new_fit < old_fit):
            old_clusters = new_clusters
            old_fit = new_fit
            print old_fit
            i = 0
    return old_clusters

def compile_dictionaries():
    "gather up all the x** dictionaries, check for redundancies"
    print "getting unique refseqs"
    with open("unique_refseqs.txt") as f:
        uniques = [line.strip() for line in f.readlines()]
    ur_dict = {}
    for filename in os.listdir('pickles'):
        if "pickle" in filename:
            print filename
            with open("pickles/" + filename) as handle:
                print "beginning unpickling"
                d = pickle.load(handle)
                print "finished unpickling"
            for key in d:
                if key in uniques:
                    uniques.remove(key)
                    ur_dict[key] = d[key]
            print len(uniques)
    return ur_dict

def uniquify_dictionaries():
    "iterate through x** pickle objects, create unique versions"
    print "getting unique refseqs"
    with open("unique_refseqs.txt") as f:
        uniques = [line.strip() for line in f.readlines()]
    for filename in os.listdir('pickles'):
        if "pickle" in filename:
            print filename
            unique_dict = {}
            with open("pickles/" + filename) as handle:
                print "beginning unpickling"
                d = pickle.load(handle)
                print "finished unpickling"
            for key in d:
                if key in uniques:
                    uniques.remove(key)
                    unique_dict[key] = d[key]
            with open("pickles/" + filename + "unique.pickle",'w') as handle:
                pickle.dump(unique_dict,handle)
            print len(uniques)

def gc_box(seq,window=50):
    return [len([b for b in box if b == 'a' or b == 't'])/window
            for box in [seq[i:i+window] for i in range(len(seq)-window)]]

def add_vectors(vs):
    return map(lambda x: x[0], zip(map(sum, zip(*vs))))

def indicate_at(seq):
    return [(b == 'a' or b == 't') * 1 for b in seq]
