class scar(object):

    def __init__(self, cigar, seq, pos):
        self.cigar = cigar
        self.seq = seq
        self.pos = pos
        return

    def __hash__(self):
        return hash((self.cigar, self.pos, self.seq))

    def __eq__(self, other):
        return (self.cigar, self.pos, self.seq) == (other.cigar, other.pos, other.seq)

    def __cmp__(self, other):
        return cmp(self.pos, other.pos)

    def __str__(self):
        return "%s\t%s\t%s" % (self.cigar, self.pos, self.seq)

