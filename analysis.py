from __future__ import division
from scipy.stats import gaussian_kde
import numpy as np
from matplotlib import pyplot as plt
import random

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

def index_distribution(promoter, score_dictionary, cutoff = 0, entropy=True):
    return [position for tf in score_dictionary[promoter]
            for (position, score) in score_dictionary[promoter][tf]
            if score > (tf_from_id(tf).entropy if entropy else cutoff)]

def tf_from_id(tf):
    return [mr for mr in more_refined if mr.ID == tf][0]

def plot_density(xs):
    density = gaussian_kde(map(float,xs))
    xs = np.arange(min(xs),max(xs),1)
    plt.plot(xs,density(xs))

def find_peaks(xs):
    return [i for i in range(len(xs)) if xs[i-1] < xs[i] > xs [i+1]]

def bin_indices(indices,n):
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
#    print "initial assignment: ", initial_assignment
    clusters = [[xs for i, xs in enumerate(xss) if initial_assignment[i] == j] for j in range(k)]
#    print "clusters: ", clusters
    centroids = map(centroid,clusters)
#    print "centroids: ", centroids
    updated_assignment = [index_of_closest_centroid(xs, centroids) for xs in xss]
#    print "updated assignment: ", updated_assignment
    while initial_assignment != updated_assignment:
        r += 1
        initial_assignment = updated_assignment
        clusters = [[xs for i, xs in enumerate(xss) if initial_assignment[i] == j] for j in range(k)]
        centroids = map(centroid,clusters)
        updated_assignment = [index_of_closest_centroid(xs, centroids) for xs in xss]
#        print updated_assignment
#        print "round: ",r
        for i,cluster in enumerate(clusters):
#            print "cluster: ", i
#            print "size: ", len(cluster)
            cent = centroids[i]
            distances = map(lambda x: distance(x,cent), cluster)
#            print "avg: ", sum(distances)/len(distances)
#            print "max: ", max(distances)
#        print
#    print "finished on round: ",r
    for i,cluster in enumerate(clusters):
#        print "cluster: ", i
#        print "size: ", len(cluster)
        cent = centroids[i]
        distances = map(lambda x: distance(x,cent), cluster)
#        print "avg: ", sum(distances)/len(distances)
#        print "max: ", max(distances)
    return assign(xss,updated_assignment)

def cluster_fitness(cluster):
    cent = centroid(cluster)
    return sum(map(lambda x: distance(x,cent)**2, cluster))
        
def fitness(clusters):
    return sum(map(cluster_fitness,clusters))
    
def assign(data,assignment):
    """group data by assignment"""
#    print min(assignment), max(assignment)
    return [[datum for i, datum in enumerate(data)
             if assignment[i] == j]
            for j in range(min(assignment),max(assignment) + 1)]

def select_from_population(pop,tournament_factor,xss):
    contender_assignments = [random.choice(pop) for i in range(tournament_factor)]
#    print len(contender_assignments), len(contender_assignments[0])
    contender_clusters = [assign(xss,ca) for ca in contender_assignments]
    (winner, fitness_winner) = sorted([(ca,fitness(assign(xss,ca))) for ca in contender_assignments],
                    key = lambda (ca,cf) : cf)[0]
#    print "select_from_population returns size: ", len(winner)
    return winner

def breed(p,q):
#    print p, q
#    print "breed gets: ", len(p), len(q)
    i = random.randint(0,len(p))
    r = p[:i] + q[i:]
#    print "breed returns size: ", len(r)
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

def levenshtein(s1, s2):#from wikibooks
    if len(s1) < len(s2):
        return levenshtein(s2, s1)
    if not s1:
        return len(s2)
 
    previous_row = xrange(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1 # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1       # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
 
    return previous_row[-1]
