import re, string, math

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

class TransfacTable(object):
    """Represents matrix.dat"""
    def __init__(self,filename):
        """Accept a filename (matrix.dat) and return a
        dictionary representation thereof"""
        lines = open(filename).readlines()
        accession_chunks = split_on(lines,matches_accession_number)[1:]
        self.entries = map(self.parse_accession_chunk, accession_chunks)

    def parse_accession_chunk(self, chunk):
        """Accept an accession chunk represented as a list of entries
        corresponding to lines of the input file and return a PSSM and
        a dictionary associating the two letter prefixes of the file
        with the line contents"""
        return TranscriptionFactor(chunk)
    
class TranscriptionFactor(object):
    """Represents the entries in matrix.dat"""
    def __init__(self, chunk):
        """Accepts a chunk of lines from matrix.dat and returns an
        object representation"""
        column_lines, other_lines = separate(matches_column, chunk)
        self.matrix = self.matrix_from_lines(column_lines)
        self.set_entropy(self.matrix)
        usable_lines = filter(self.usable_line,other_lines)
        self.attr_dict = {}
        for line in usable_lines:
            tlid = line[:2]
            attr = line[2:].strip()
            if not tlid in self.attr_dict:
                self.attr_dict[tlid] = [attr]
            else:
                self.attr_dict[tlid].append(attr)
        self.attr_dict["pssm"] = self.pssm_from_matrix(self.matrix)
        for key in self.attr_dict:
            val = self.attr_dict[key]
            exec("""storable_val = val if len(val) > 1 else val[0];self.{0} = storable_val""".format(key,val))
        
    def pssm_from_matrix(self, matrix,background_probs = (0.25,)*4):
        """Accept count matrix (as nested list) and return pssm (as nested list)"""
        def convert_column(col):
            return [safe_log2(c/p) for (c,p) in zip(normalize(col),background_probs)]
        return [convert_column(col) for col in matrix]

    def set_entropy(self,matrix):
        """Accept count matrix as nested list and set entropy"""
        def col_entropy(col):
            return -sum(c*safe_log2(c) for c in normalize(col))
        self.entropy = sum(col_entropy(col) for col in matrix)
        print self.entropy
        
    def usable_line(self,line):
        bad_tlids = ['XX', #separates entries in chunk
                     '//', #separates chunks in file
                     'P0'] #signals beginning of PSWM
        return not any(line.startswith(tlid) for tlid in bad_tlids)

    def matches_text(self,text):
        """Search object attributes for text"""
        def contains_text(attr):
            if type(attr) is list:
                return any(text.upper() in line.upper() for line in attr)
            else:
                return text.upper() in attr.upper()
        attrs = [attr for attr in dir(self) if attr[0] in string.uppercase]
        return any(contains_text(getattr(self,attr)) for attr in attrs)

    def matrix_from_lines(self,lines):
        """Convert raw column lines into matrix.  Assumes all lines are
        column lines"""
        return [map(float,matches_column(line).groups()) for line in lines]
