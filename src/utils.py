from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import gaussian_kde


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

def plot_density(xs):
    density = gaussian_kde(map(float,xs))
    xs = np.arange(min(xs),max(xs),1)
    plt.plot(xs,density(xs))

def add_vectors(vs):
    return map(lambda x: x[0], zip(map(sum, zip(*vs))))

def safe_log2(x):
    """Implements log2, but defines log2(0) = 0"""
    return math.log(x,2) if x > 0 else 0

def split_on(xs, pred):
    """Split xs into a list of lists each beginning with the next x
    satisfying pred, except possibly the first"""
    indices = [i for (i,v) in enumerate(xs) if pred(v)]
    return [xs[i:j] for (i,j) in zip([0]+indices,indices+[len(xs)]) if i != j]

def separate(pred, lst):
    """separates lst into a list of elements satisfying pred and a list of 
    elements not satisfying it.
    """
    sheep = []
    goats = []
    for elem in lst:
        if pred(elem):
            sheep.append(elem)
        else:
            goats.append(elem)
    return (sheep, goats)

def normalize(xs):
    return map(lambda(x): x/float(sum(xs)),xs)
