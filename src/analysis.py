from __future__ import division
import random, os, pickle
import utils

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

def indicate_at(seq):
    return [(b == 'a' or b == 't') * 1 for b in seq]
