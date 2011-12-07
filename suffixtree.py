def lcp(xs):
    if any([x=="" for x in xs]):
        return ""
    else:
        first = xs[0][0]
        if all([x.startswith(first) for x in xs]):
            return first + lcp([x[1:] for x in xs])
        else:
            return ""

def non_trivially_startswith(word,prefix):
    if prefix == "":
        return word == ""
    else:
        return word.startswith(prefix)
    
def chomp_prefix(pre,word):
    """Return suffix suf such that pre + suf = word"""
    assert(word.startswith(pre))
    return word[word.index(pre)+len(pre):]

def adjacent_pairs(xs):
    return zip([None] + xs, xs + [None])[1:len(xs)]

def suffixes(word):
    "return suffixes of word, including empty suffix"
    return [word[i:] for i in range(len(word)+1)]
    
class SuffixTree(object):
    def __init__(self, text_or_suffixes,node_text=None):
        print text_or_suffixes
        self.node_text = node_text
        if isinstance(text_or_suffixes,str):
            self.text = text_or_suffixes
            self.alphabet = list(set(self.text))
            self.suffixes = suffixes(self.text)
        else:
            self.suffixes = text_or_suffixes
            self.alphabet = list(set("".join(self.suffixes)))
        self.partition = [[s for s in self.suffixes if s.startswith(c)]
                          for c in self.alphabet]
        self.children = {}
        for p in self.partition:
            if p:
                prefix = lcp(p)
                self.children[prefix] = SuffixTree([chomp_prefix(prefix,word) for word in p],prefix)

    def traverse(self,word):
        if word == "":
            return self
        else:
            next_link_list = [child for child in self.children.keys()
                              if word.startswith(child)]
            if next_link_list == []:
                return False #couldn't traverse suffix tree on input
            else:
                [child] = next_link_list
                return self.children[child].traverse(chomp_prefix(child,word))
    
    def pprint(self):
        print self.suffixes
        for suffix in self.children.keys():
            print suffix
            self.children[suffix].pprint()

def lexicographical_cmp(xs,ys):
    print(xs,ys)
    if not xs or not ys:
        return cmp(len(ys),len(xs))
    elif cmp(xs[0],ys[0]) != 0:
        return cmp(xs[0],ys[0])
    else:
        return lexicographical_cmp(xs[1:],ys[1:])
    
class ESA(object):
    def __init__(self,word):
        raw_suffixes = suffixes(word)
        self.suffixes = sorted(raw_suffixes,cmp=lexicographical_cmp)
        self.suf = [raw_suffixes.index(s) for s in self.suffixes]
        self.lcp = [0] + [len(lcp([x,y])) for (x,y) in adjacent_pairs(self.suffixes)]
        self.skp = [([j for j in range(i+1,len(word)+1)
                      if esa.lcp[j] < esa.lcp[i]] +
                     [len(word)+1])[0]
                    for i in range(len(word)+1)]
        self.suffix_tree = Suffix_Tree(word)

class ST(object):
    def __init__(self,suffixes,prefix=""):
        self.prefix = prefix
        if len(suffixes) == 1:
            self.leaf = self.prefix + suffixes[0]
        else:
            self.alphabet = list(set(("".join(suffixes)))) + [""]
            self.partition = [[s for s in suffixes if non_trivially_startswith(s,c)]
                              for c in self.alphabet]
            self.children = {}
            for p in self.partition:
                if p:
                    prefix = lcp(p)
                    self.children[prefix] = ST([chomp_prefix(prefix,word) for word in p],self.prefix+prefix)

    def set_index(self,master_suffixes):
        """Recursively set the index of the suffix table that the leaf identifies"""
        print self.prefix
        if self.has_leaf():
            self.suffix_index = master_suffixes.index(self.leaf)
        else:
            for child in self.children.values():
                child.set_index(master_suffixes)

    def set_skp(self,skp_table):
        """Recursively set the index of the suffix table identified by skp_table"""
        if self.has_leaf():
            self.skp = skp_table[self.suffix_index]
        else:
            for child in self.children.values():
                child.set_skp(skp_table)

    def has_leaf(self):
        return hasattr(self,leaf)
