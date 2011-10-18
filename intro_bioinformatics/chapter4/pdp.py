# problem 4.2
# Partial Digest Problem
# Steven Skiena's algorithm in 1990
import copy
debug = False

def PartialDigest(L):
    width = max(L)
    for i, x in enumerate(L):
        if x == width:
            L = L[:i]+L[i+1:]
            break
    X = [0, width]
    Place(L, X)

# generate a set whose elements
# are absolute difference from
# each element in L to y
def diffGen(y, X):
    return [abs(y-x) for x in X]

# subtract multiset S from multiset L
def diffSet(S, L):
    for s in S:
        for i,l in enumerate(L):
            if s == l:
                L = L[:i]+L[i+1:]
                break
    return L

# used deepcopy to cope with lists
def Place(L, X):
    if L==[]:
        X1 = copy.deepcopy(X)
        X1.sort()
        print X1
        return True

    width = max(X)
    if debug:
        print L
        print X
        print '------------'

    L1 = copy.deepcopy(L)
    X1 = copy.deepcopy(X)

    y = max(L)
    diff = diffGen(y, X)
    if set(diff) <= set(L):
        X += [y]
        L = diffSet(diff, L)
        flag = Place(L,X)

    L = L1
    X = X1
    y = width - max(L)
    diff = diffGen(y, X)
    if set(diff) <= set(L):
        X += [y]
        L = diffSet(diff,L)
        flag = Place(L,X)

    return


if __name__ == "__main__":

    L = [1,1,1,2,2,3,3,3,4,4,5,5,6,6,6,9,9,10,11,12,15]
    print 'query', L, '\n'
    PartialDigest(L)


