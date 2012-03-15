#! /usr/bin/env python

class Affine:
    def __init__(self, a, X, b):
        self.a = a
        self.X = X
        self.b = b
    def __repr__(self):
        return str(self)
    def __str__(self):
        if self.a == 1 and self.b == 0:
            return 'Aff(%s)' % (self.X)
        elif self.a == 1:
            return 'Aff(%s+%d)' % (self.X, self.b)
        elif self.b == 0:
            return 'Aff(%d%s)' % (self.a, self.X)
        else:
            return 'Aff(%d%s+%d)' % (self.a, self.X, self.b)
    def __eq__(self, other):
        return self.a == other.a \
                and self.X == other.X \
                and self.b == other.b
class Subscript:
    def __init__(self, L, R):
        self.L = L
        self.R = R
    def is_ziv(self):
        return self.L.X is None or self.R.X is None
    def is_siv(self):
        return not(self.is_ziv()) and self.L.X == self.R.X
    def is_miv(self):
        return not(self.is_ziv()) and not(self.is_siv())
    def __repr__(self):
        return str(self)
    def __str__(self):
        return "Sub<%s %s>" % (self.L, self.R)
    def __eq__(self, other):
        return self.L == other.L and self.R == other.R
def partition( subs, indexes ):
    """
    :param subs: the list of subscripts
    :param indexes: the loop index variables
    :returns: a list of lists containing the subscript partitionings
    """
    results = []
    def has_index(s, idx):
        for it in s:
            if it.L.X == idx or it.R.X == idx:
                return True
        return False
    for s in subs:
        results.append( [ s ] )
    for idx in indexes:
        k = None
        killed = []
        for j in xrange(0, len(results)):
            r = results[j]
            if has_index(r, idx):
                if k is None:
                    k = j
                else:
                    results[k] = results[k] + r
                    killed.append( j )
        for x in killed:
            del results[x]
        del killed
    return results
def siv_test( sub ):
    if sub.L.X != sub.R.X:
        raise Exception("%s is not SIV" % sub)
    return sub.L.b - sub.R.b

# A[I+1][I][K] = A[I][J][J]
example = [ \
    Subscript( Affine(1,'I',1), Affine(1,'I',0) ),
    Subscript( Affine(1,'I',0), Affine(1,'J',0) ),
    Subscript( Affine(1,'K',0), Affine(1,'J',0) ),
    ]
for p in partition( example, ['I','J','K'] ):
    print "Partition=", p
    for it in p:
        if it.is_siv():
            print "SIV.d=", siv_test( it )

