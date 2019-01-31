#!/usr/bin/env python3
from math import sqrt

def d_squared_euclidean(a,b):
    assert len(a) == len(b), "Can't compute Euclidean distance on different-dimensional points"
    return sum((a[i]-b[i])**2 for i in range(len(a)))

def d_euclidean(a,b):
    return sqrt(d_squared_euclidean(a,b))

def align(pos1, pos2, d=d_euclidean):
    M = [[(float('inf'),None) for j in range(len(pos2)+1)] for i in range(len(pos1)+1)] # (score,backtrack) tuples, where 0 = left, 1 = up, 2 = diag
    M[0][0] = (0,None) # top-left corner
    for i in range(1,len(M)): # left column (iterate over rows)
        M[i][0] = (0,1)
    for j in range(1,len(M[0])): # top row (iterate over columns)
        M[0][j] = (0,0)
    for i in range(1,len(M)):
        for j in range(1,len(M[i])):
            assert len(pos1[i-1]) == len(pos2[j-1]), "Must have same number of points in each file"
            opts = [
                float('inf'), # index 0: left (gap in pos1) TODO
                float('inf'), # index 1: up (gap in pos2) TODO
                M[i-1][j-1][0] + sum(d(pos1[i-1][c],pos2[j-1][c]) for c in range(len(pos1[i-1]))) # index 2: diagonal (match)
            ]
            # incorporate free gaps at ends of sequences TODO
            for o,v in enumerate(opts):
                if v < M[i][j][0]:
                    M[i][j] = (v,o)
    # TODO need to think of gap penalty. one idea: by introducing a gap, the "gap" is actually some average between the two things it's between. How to handle multiple gaps, though?
    align1 = list(); align2 = list(); i = len(pos1); j = len(pos2)
    while i > 0 and j > 0:
        choice = M[i][j][1]
        assert choice in {0,1,2}, "Invalid choice: %s" % str(choice)
        if choice == 0: # left (use pos2, gap in pos1)
            align1.append('---')
            align2.append(pos2[j-1]); j -= 1
        elif choice == 1: # up (use pos1, gap in pos2)
            align1.append(pos1[i-1]); i -= 1
            align2.append('---')
        else: # diag (use pos1 and pos2)
            align1.append(pos1[i-1]); i -= 1
            align2.append(pos2[j-1]); j -= 1
    while i > 0: # move all the way up
        align1.append(pos1[i-1]); i -= 1
        align2.append('---')
    while j > 0: # move all the way left
        align1.append('---')
        align2.append(pos2[j-1]); j -= 1
    return align1[::-1],align2[::-1]

if __name__ == "__main__":
    from sys import argv
    if len(argv) != 3:
        print("USAGE: %s <positions_1> <positions_2>" % argv[0]); exit(1)
    pos = [list(),list()]
    for fi,f in enumerate(argv[1:]):
        for l in open(f):
            pos[fi].append(list())
            cols = [float(e) for e in l.strip().split(',')]
            for i in range(0,len(cols),3):
                pos[fi][-1].append(tuple(cols[i:i+3]))
    a1,a2 = align(pos[0],pos[1])
    for l in a1:
        if isinstance(l,str):
            print(l.strip())
        else:
            print(','.join('%f,%f,%f'%xyz for xyz in l))
    print('===')
    for l in a2:
        if isinstance(l,str):
            print(l.strip())
        else:
            print(','.join('%f,%f,%f'%xyz for xyz in l))
