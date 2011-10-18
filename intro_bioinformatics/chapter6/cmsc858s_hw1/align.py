#!/usr/local/bin/python
# implementation of alignment problems described in lecture notes
# enable debug mode will show the computation record table and trace table
debug = False
import copy
from sys import *

## Global variables needed by several alignment algorithms
global g_match_score
global g_mismatch_score
global g_gap_open
global g_gap_extend

def score(ch1, ch2):
    global g_match_score
    global g_mismatch_score

    if ch1 == ch2:
        return g_match_score
    else:
        return g_mismatch_score

def score_gap(ch1, ch2):
    global g_match_score
    global g_mismatch_score
    
    if ch1 == ch2:
        return g_match_score
    else:
        return g_mismatch_score

# give all possible alignment result
# x, y are given as the stop coordinate on the trace table
def recover(trace, x, y, S, T, buff):
    if trace[x][y] == ['e']:
        print ''.join(buff[0])
        print ''.join(buff[1])
        print 
        #print '---------------------------------------\n'
        return

    for elem in trace[x][y]:
        if elem == 'l':
            buff[0] = ['-'] + buff[0]
            buff[1] = [T[y-1]] + buff[1]
            recover(trace, x, y-1, S, T, buff)
            buff[0] = buff[0][1:]
            buff[1] = buff[1][1:]
        elif elem == 'u':
            buff[0] = [S[x-1]] + buff[0]
            buff[1] = ['-'] + buff[1]
            recover(trace, x-1, y, S, T, buff)
            buff[0] = buff[0][1:]
            buff[1] = buff[1][1:]
        else:
            buff[0] = [S[x-1]] + buff[0]
            buff[1] = [T[y-1]] + buff[1]
            recover(trace, x-1, y-1, S, T, buff)
            buff[0] = buff[0][1:]
            buff[1] = buff[1][1:]

def result(trace, S, T):
    x = len(trace)-1
    y = len(trace[0])-1
    #print '\n-----------Possible alignment(s)-----------\n'
    recover(trace, x, y, S, T, [[], []])

# global alignment
def g_align(S, T):
    # table records values
    # trace records results
    table = []
    trace = []
    for i in range(0, len(S)+1):
        table += [[]]
        trace += [[]]
        for j in range(0, len(T)+1):
            table[i] += [0]
            trace[i] += [[]]

    for i in range(0, len(S)+1):
        for j in range(0, len(T)+1):
            #trace[i][j] = []
            if i==0 and j==0:
                table[i][j] = 0
                trace[i][j] += ['e']
            elif i==0:
                table[i][j] = table[i][j-1]-1
                trace[i][j] += ['l']
            elif j==0:
                table[i][j] = table[i-1][j]-1
                trace[i][j] += ['u']
            else:
                # take the max of
                vu = table[i-1][j] + score(S[i-1], '-')
                vl = table[i][j-1] + score('-', T[j-1])
                vd = table[i-1][j-1] + score(S[i-1], T[j-1])
                table[i][j] = max(max(vl, vu), max(vu, vd))
                if table[i][j] == vl:
                    trace[i][j] += ['l']
                if table[i][j] == vu:
                    trace[i][j] += ['u']
                if table[i][j] == vd:
                    trace[i][j] += ['d']

    if debug:
        for i in range(0, len(S)+1):
            print table[i]
        for i in range(0, len(S)+1):
            print trace[i]

    # display the result
    result(trace, S, T)


# Smith-Waterman local alignment
def l_align(S, T):
    table = []
    trace = []
    for i in range(0, len(S)+1):
        table += [[]]
        trace += [[]]
        for j in range(0, len(T)+1):
            table[i] += [0]
            trace[i] += [[]]

    best_score = -1
    best_i = -1
    best_j = -1
    for i in range(0, len(S)+1):
        for j in range(0, len(T)+1):
            if i == 0 or j == 0:
                table[i][j] = 0
                trace[i][j] += ['e']
            else:
                vu = table[i-1][j] + score(S[i-1], '-')
                vl = table[i][j-1] + score('-', T[j-1])
                vd = table[i-1][j-1] + score(S[i-1], T[j-1])
                v0 = 0

                table[i][j] = max( max(vu, vl), max(v0, vd) )

                if table[i][j] > best_score:
                    best_score = table[i][j]
                    best_i = [i]
                    best_j = [j]
                elif table[i][j] == best_score:
                    best_i += [i]
                    best_j += [j]

                if table[i][j] == vl:
                    trace[i][j] += ['l']
                if table[i][j] == vu:
                    trace[i][j] += ['u']
                if table[i][j] == vd:
                    trace[i][j] += ['d']
                if table[i][j] == v0:
                    trace[i][j] += ['e']

    if debug:
        for entry in table:
            print entry
        print '\n'
        for entry in trace:
            print entry
        print '\n'

    for [ier, jer] in zip(best_i, best_j):
        recover(trace, ier, jer, S, T, [[],[]])


# global alignment with affine gap penalty
infty = 10e9
def g_align_gap(S, T):    
    global g_gap_open
    global g_gap_extend
    Wg = -g_gap_open
    Ws = -g_gap_extend

    # initialization
    V = []
    trace = []
    F = []
    E = []
    for i in range(0, len(S)+1):
        V += [[]]
        trace += [[]]
        F += [[]]
        E += [[]]
        for j in range(0, len(T)+1):
            V[i] += [0]
            trace[i] += [[]]
            F[i] += [0]
            E[i] += [0]

    # iteration
    V[0][0] = 0
    trace[0][0] = ['e']

    for i in range(0, len(S)+1):
        for j in range(0, len(T)+1):
            if i == 0:
                if j!=0:
                    V[i][j] = -infty
                    F[i][j] = -infty
                    E[i][j] = -Wg - (j+1)*Ws
                    trace[i][j] = ['l']

            elif j == 0:
                V[i][j] = -infty
                F[i][j] = -Wg - (i+1)*Ws
                E[i][j] = -infty
                trace[i][j] = ['u']

            else:
                # compute the general recurrence
                mv = score_gap(S[i-1], T[j-1])
                G =  V[i-1][j-1] + mv
                vl = E[i-1][j-1] + mv
                vu = F[i-1][j-1] + mv
                V[i][j] = max( max(G, vu), vl)                
                F[i][j] = max( max(F[i-1][j]-Ws, E[i-1][j]-Wg-Ws),
                               V[i-1][j]-Wg-Ws)
                E[i][j] = max( max(E[i][j-1]-Ws, F[i][j-1]-Wg-Ws),
                               V[i][j-1]-Wg-Ws)
                
                # track the result
                curr_max = max(V[i][j], max(F[i][j], E[i][j]))
                if curr_max == E[i][j]:
                    trace[i][j] += ['l']
                if curr_max == F[i][j]:
                    trace[i][j] += ['u']
                if curr_max == V[i][j]:
                    trace[i][j] += ['d']

    debug = 0
    if debug:
        for entry in V:
            print entry
        print '\n'
        for entry in F:
            print entry
        print '\n'
        for entry in E:
            print entry
        print '\n'
        for entry in trace:
            print entry
        print '\n'
        
    ns = len(S)
    nt = len(T)
    print 'Score =', max(V[ns][nt], max(E[ns][nt], F[ns][nt]))
    result(trace, S, T)


# optimal fitting alignment as presented in blue book p215

# the fitting score for two nuclitides
def fit_score(ch1, ch2):
    if ch1 == '-' or ch2 == '-' or ch1 != ch2:
        return -1
    else:
        return 1

# the alignment tool for computing the result
# all valid alignment will be displayed
def fit_align(S, T):
    table = []
    trace = []

    for i in range(0, len(S)+1):
        table += [[]]
        trace += [[]]
        for j in range(0, len(T)+1):
            table[i] += [0]
            trace[i] += [[]]

    for i in range(0, len(S)+1):
        for j in range(0, len(T)+1):
            if i == 0 or j == 0:
                if i == 0 and j == 0:
                    table[i][j] = 0
                    trace[i][j] += ['e']
                else:
                    if i == 0: # matching nothing with prefix of T
                        table[i][j] = -7e21
                        trace[i][j] += ['e']
                    else: # matching prefix of i to nothing
                        table[i][j] = 0
                        trace[i][j] += ['e']

            else:
                vu = table[i-1][j] + fit_score(S[i-1], '-')
                vl = table[i][j-1] + fit_score('-', T[j-1])
                vd = table[i-1][j-1] + fit_score(S[i-1], T[j-1])

                table[i][j] = max( max(vu, vl), vd )

                if table[i][j] == vl:
                    trace[i][j] += ['l']
                if table[i][j] == vu:
                    trace[i][j] += ['u']
                if table[i][j] == vd:
                    trace[i][j] += ['d']

    if debug:
        for entry in table:
            print entry
        print '\n'
        for entry in trace:
            print entry
        print '\n'

    best_score = -1
    best_i = -1
    for i in range(0, len(S)+1):
        this_score = table[i][len(T)]
        if this_score > best_score:
            best_score = this_score
            best_i = i
    print 'Fitting Alignment for the two sequences:'
    print 'S:', S
    print 'T:', T
    print 'best alignment score:', best_score

    recover(trace, best_i, len(T), S, T, [[],[]])


# overlapping alignment as presented in blue book p215

# the overlapping score for two nuclitides
def overlap_score(ch1, ch2):
    if ch1 == '-' or ch2 == '-' or ch1 != ch2:
        return -1
    else:
        return 2

# the alignment tool for computing the result
# all valid alignment will be displayed
def overlap_align(S, T):
    table = []
    trace = []

    for i in range(0, len(S)+1):
        table += [[]]
        trace += [[]]
        for j in range(0, len(T)+1):
            table[i] += [0]
            trace[i] += [[]]

    for i in range(0, len(S)+1):
        for j in range(0, len(T)+1):
            if i == 0 or j == 0:
                if i == 0 and j == 0:
                    table[i][j] = 0
                    trace[i][j] += ['e']
                else:
                    if i == 0:
                        table[i][j] = -7e21
                        trace[i][j] += ['e']
                    else:
                        table[i][j] = 0
                        trace[i][j] += ['e']

            else:
                vu = table[i-1][j] + overlap_score(S[i-1], '-')
                vl = table[i][j-1] + overlap_score('-', T[j-1])
                vd = table[i-1][j-1] + overlap_score(S[i-1], T[j-1])

                table[i][j] = max( max(vu, vl), vd )

                if table[i][j] == vl:
                    trace[i][j] += ['l']
                if table[i][j] == vu:
                    trace[i][j] += ['u']
                if table[i][j] == vd:
                    trace[i][j] += ['d']

    if debug:
        for entry in table:
            print entry
        print '\n'
        for entry in trace:
            print entry
        print '\n'

    best_score = -1
    best_j = -1
    for j in range(0, len(T)+1):
        this_score = table[len(S)][j]
        if this_score > best_score:
            best_score = this_score
            best_j = j
    print 'Overlapping Alignment for the two sequences:'
    print 'S:', S
    print 'T:', T
    print 'best alignment score:', best_score

    recover(trace, len(S), best_j, S, T, [[],[]])


# semiglobal alignment as presented in blue book p216

# the semiglobal score for two nuclitides
def semiglobal_score(ch1, ch2):
    if ch1 == '-' or ch2 == '-' or ch1 != ch2:
        return -1
    else:
        return 1

# the alignment tool for computing the result
# all valid alignment will be displayed
def semiglobal_align(S, T):
    table = []
    trace = []

    for i in range(0, len(S)+1):
        table += [[]]
        trace += [[]]
        for j in range(0, len(T)+1):
            table[i] += [0]
            trace[i] += [[]]

    for i in range(0, len(S)+1):
        for j in range(0, len(T)+1):
            if i == 0 or j == 0:
                if i == 0 and j == 0:
                    table[i][j] = 0
                    trace[i][j] += ['e']
                else:
                    if i == 0:
                        table[i][j] = 0
                        trace[i][j] += ['l']
                    else:
                        table[i][j] = 0
                        trace[i][j] += ['u']

            else:
                vu = table[i-1][j] + semiglobal_score(S[i-1], '-')
                vl = table[i][j-1] + semiglobal_score('-', T[j-1])
                vd = table[i-1][j-1] + semiglobal_score(S[i-1], T[j-1])

                table[i][j] = max( max(vu, vl), vd )

                if table[i][j] == vl:
                    trace[i][j] += ['l']
                if table[i][j] == vu:
                    trace[i][j] += ['u']
                if table[i][j] == vd:
                    trace[i][j] += ['d']

    if debug:
        for entry in table:
            print entry
        print '\n'
        for entry in trace:
            print entry
        print '\n'

    best_jscore = -1
    best_j = -1
    best_iscore = -1
    best_i = -1
    for i in range(0, len(S)+1):
        this_score = table[i][len(T)]
        if this_score > best_iscore:
            best_iscore = this_score
            best_i = i
    for j in range(0, len(T)+1):
        this_score = table[len(S)][j]
        if this_score > best_jscore:
            best_jscore = this_score
            best_j = j
    best_score = max(best_jscore, best_iscore)

    print 'Semiglobal Alignment for the two sequences:'
    print 'S:', S
    print 'T:', T
    print 'best alignment score:', best_score

    if best_score == best_iscore:
        for i in range(best_i+1, len(S)):
            trace[i][len(T)] = ['u']
        recover(trace, len(S), len(T), S, T, [[],[]])
    if best_score == best_jscore:
        for j in range(best_j+1, len(T)):
            trace[len(S)][j] = ['l']
        recover(trace, len(S), len(T), S, T, [[],[]])




# inexact repeat alignment
def is_overlap(start_list, end):
    tmp = []
    for start in start_list:
        if (start[0] > end[1] or start[1] > end[0])\
            and start not in tmp:
            tmp += [start]
    return tmp

# the alignment tool for computing the result
# all valid alignment will be displayed
def inexact_repeat_align(S):
    table = []
    trace = []
    start_pos = [[[]]*(len(S)+1)]*(len(S)+1)

    best_score = -7e14
    best_j = -1
    best_i = -1

    # initialize the tables
    for i in range(0, len(S)+1):
        table += [[]]
        trace += [[]]
        for j in range(0, len(S)+1):
            table[i] += [0]
            trace[i] += [[]]

    for i in range(0, len(S)+1):
        for j in range(0, len(S)+1):
            if i == 0 or j == 0:
                table[i][j] = 0
                trace[i][j] += ['e']
                start_pos[i][j] = [[i,j]]

            else:
                sp_l = is_overlap(start_pos[i][j-1], [i,j])
                sp_u = is_overlap(start_pos[i-1][j], [i,j])
                sp_d = is_overlap(start_pos[i-1][j-1], [i,j])
                vu = (table[i-1][j] + score(S[i-1], '-'))\
                     if sp_u != [] else -7e14
                vl = (table[i][j-1] + score('-', S[j-1]))\
                     if sp_l != [] else -7e14
                vd = (table[i-1][j-1] + score(S[i-1], S[j-1]))\
                     if sp_d !=[] else -7e14

                table[i][j] = max( max(vu, vl), max(vd, 0) )

                if table[i][j] == vl:
                    trace[i][j] += ['l']
                    start_pos[i][j] += sp_l

                if table[i][j] == vu:
                    trace[i][j] += ['u']
                    start_pos[i][j] += sp_u

                if table[i][j] == vd:
                    trace[i][j] += ['d']
                    start_pos[i][j] += sp_d

                if table[i][j] == 0:
                    trace[i][j] += ['e']
                    start_pos[i][j] += [[i,j]]

            this_score = table[i][j]
            if this_score > best_score:
                best_score = this_score
                best_i = i
                best_j = j

    if debug:
        for entry in table:
            print entry
        print
        for entry in trace:
            print entry
        print

    print 'Inexact Repeat problem for the two sequences:'
    print 'S:', S
    print 'best alignment score:', best_score

    recover(trace, best_i, best_j, S, S, [[],[]])


# shortest supersequence problem presented in blue book p217
# the alignment tool for computing the result
# all valid alignment will be displayed
def super_sequence(S, T):
    table = []
    trace = []

    for i in range(0, len(S)+1):
        table += [[]]
        trace += [[]]
        for j in range(0, len(T)+1):
            table[i] += [0]
            trace[i] += [[]]

    for i in range(0, len(S)+1):
        for j in range(0, len(T)+1):
            if i == 0 or j == 0:
                if i == 0 and j == 0:
                    table[i][j] = 0
                    trace[i][j] += ['e']
                else:
                    if i == 0:
                        table[i][j] = j
                        trace[i][j] += ['l']
                    else:
                        table[i][j] = i
                        trace[i][j] += ['u']

            else:
                vu = table[i-1][j] + 1
                vl = table[i][j-1] + 1
                vd = 7e21
                if S[i-1] == T[j-1]:
                    vd = table[i-1][j-1] + 1

                table[i][j] = min( min(vu, vl), vd )

                if table[i][j] == vl:
                    trace[i][j] += ['l']
                if table[i][j] == vu:
                    trace[i][j] += ['u']
                if table[i][j] == vd:
                    trace[i][j] += ['d']

    if debug:
        for entry in table:
            print entry
        print '\n'
        for entry in trace:
            print entry
        print '\n'

    print 'Shortest Supersequence for the two sequences:'
    print 'S:', S
    print 'T:', T
    print 'shortest length:', table[-1][-1]

    result(trace, S, T)


# virus infection problem presented in blue book p219

def copy_check(num, char):
    constraint = {'A':[1,5],
                  'C':[1,10],
                  'G':[1, infty],
                  'T':[1, infty]}
    thresh = constraint[char]
    if num >= thresh[0] and num <= thresh[1]:
        return True
    return False

# try to see whether the string T is an infected version of S
def is_infected(S, T):
    infty = 7e21
    num_row = len(S) + 1
    num_col = len(T) + 1
    table = []
    copy = []

    # initialize the tables
    for i in range(0, num_row):
        table += [[]]
        copy += [[]]
        for j in range(0, num_col):
            table[i] += [False]
            copy[i] += [0]
    table[0][0] = True

    for i in range(1, num_row):
        for j in range(1, num_col):
            vl = table[i][j-1] if T[j-1] == S[i-1] else False
            vd = table[i-1][j-1] if T[j-1] == S[i-1] else False

            table[i][j] = False
            if vl and copy_check(copy[i][j-1]+1, S[i-1]):
                table[i][j] = True
                copy[i][j] = copy[i][j-1]+1
            if vd:
                table[i][j] = True
                copy[i][j] = 0

    if debug:
        for entry in table:
            print entry
        print '\n'
        for entry in copy:
            print entry
        print '\n'

    print 'Infection check for the two sequences:'
    print 'S:', S
    print 'T:', T
    print 'result:', table[-1][-1]
    print


# solve the homo-edit distance problem
def homo_dist(S, T):
    num_row = len(S) + 1
    num_col = len(T) + 1
    infty = 7e21
    table = []
    state = []

    for i in range(0, num_row):
        table += [[]]
        state += [[]]
        for j in range(0, num_col):
            table[i] += [infty]
            state[i] += ['m']

    table[0][0] = 0

    for i in range(1, num_row):
        for j in range(1, num_col):

            if j == 1:
                vl = i
            else:
                vl = table[i][j-1] +\
                     (0 if state[i][j-1]=='i' and T[i-2] == T[i-1] else 1)
            if i == 1:
                vu = j
            else:
                vu = table[i-1][j] +\
                     (0 if state[i-1][j]== 'd' and S[i-2] == S[i-1] else 1)
            vd = table[i-1][j-1] if S[i-1] == T[j-1] else infty

            table[i][j] = min( min(vl, vd), vu )
            if table[i][j] == vl:
                state[i][j] = 'i'
            if table[i][j] == vd:
                state[i][j] = 'd'
            if table[i][j] == vu:
                state[i][j] = 'm'

    for entry in table:
        print entry
    for entry in state:
        print entry


    print 'Homo-Edit Distance of the two strings:'
    print 'S:', S
    print 'T:', T
    print 'distance:', table[-1][-1]


# find the longest increasing or decreasing subsequence
def get_lids(S, increasing = True):
    n = len(S)
    infty = 7e21
    tmp = []
    tmpS = []
    for i, elem in enumerate(S):
        tmpS += [int(elem)]
        tmp += [int(elem)]
    S = tmpS
    T = tmp
    T.sort()
    if not increasing:
        T.reverse()
    table = []
    trace = []

    # initialization
    for i in range(0, n+1):
        table += [[]]
        trace += [[]]
        for j in range(0, n+1):
            if i == 0 or j == 0:
                trace[i] += [['e']]
            else:
                trace[i] += [[]]
            table[i] += [0]
    table[0][0] = 0

    for i in range(1, n+1):
        for j in range(1, n+1):
            vl = table[i][j-1]
            vu = table[i-1][j]
            vd = table[i-1][j-1] + (1 if S[i-1] == T[j-1] else 0)

            table[i][j] = max(vl, max(vu, vd))

            if table[i][j] == vl:
                trace[i][j] += ['l']
            if table[i][j] == vu:
                trace[i][j] += ['u']
            if table[i][j] == vd and S[i-1] == T[j-1]: # for a potential error
                trace[i][j] += ['d']

    if debug:
        for entry in table:
            print entry
        print
        for entry in trace:
            print entry

    for i, elem in enumerate(S):
        S[i] = str(elem)
        T[i] = str(T[i])
    if increasing:
        print 'Longest Increasing Subsequence of sequence S'
    else:
        print 'Longest Decreasing Subsequence of sequence S'
    print 'S:', S
    print 'length:', table[-1][-1]
    result(trace, S, T)
    print


# the longest 2-increasing problem
def get_l2inc(S):
    n = len(S) + 1
    m = n

    # initialization for 1st iteration
    T = []
    for elem in S:
        T += [int(elem)]
    T.sort()
    print T
    for i,elem in enumerate(T):
        T[i] = int(elem)

    tmp = []
    for elem in S:
        tmp += [int(elem)]
    S = tmp

    record = []
    for i in range(0, m):
        record += [[]]
        for j in range(0, n):
            record[i] += [[]]
            for k in range(0, n-1):
                record[i][j] += [0]

    # first iteration
    for i in range(1, m):
        for j in range(1, n):
            vl = sum(record[i][j-1])
            vu = sum(record[i-1][j])
            if S[i-1] == T[j-1]:
                for k in range(0, n-1):
                    record[i][j][k] = record[i-1][j-1][k]
                record[i][j][j-1] = 1
            vd = sum(record[i][j])

            maxer = max( vl, max(vu, vd) )

            if maxer == vl:
                for k in range(0, n-1):
                    record[i][j][k] = record[i][j-1][k]
            elif maxer == vu:
                for k in range(0, n-1):
                    record[i][j][k] = record[i-1][j][k]
            elif maxer == vd:
                if S[i-1] != T[j-1]:
                    for k in range(0, n-1):
                        record[i][j][k] = record[i-1][j-1][k]

    # initiatlization for 2nd iteration
    to_del = record[-1][-1]
    tmp = []
    for i, elem in enumerate(S):
        if to_del[elem-1] == 0:
            tmp += [elem]
    S = tmp
    tmp = []
    for i, elem in enumerate(T):
        if to_del[elem-1] == 0:
            tmp += [elem]
    T = tmp
    m = len(S) + 1
    n = len(T) + 1
    print to_del
    print S, T

    # second iteration
    table = []
    trace = []
    for i in range(0, m):
        table += [[]]
        trace += [[]]
        for j in range(0, n):
            table[i] += [0]
            if i == 0 or j == 0:
                trace[i] += [['e']]
            else:
                trace[i] += [[]]

    for i in range(1, m):
        for j in range(1, n):
            vl = table[i][j-1]
            vu = table[i-1][j]
            vd = table[i-1][j-1] + 1 if S[i-1] == T[j-1] else 0

            table[i][j] = max( vl, max(vu, vd) )

            if table[i][j] == vl:
                trace[i][j] += ['l']
            if table[i][j] == vu:
                trace[i][j] += ['u']
            if table[i][j] == vd:
                trace[i][j] += ['d']

    if False:
        for entry in table:
            print entry
        print
        for entry in trace:
            print entry

    for i, elem in enumerate(S):
        S[i] = str(elem)
    for i, elem in enumerate(T):
        T[i] = str(elem)

    print 'longest 2 increasing problem for sequence S'
    print 'S:', ''.join(S[0:len(S)/2])
    print 'length:', table[-1][-1]
    print 'base case:'
    base_case = []
    for i, elem in enumerate(to_del):
        if elem:
            base_case += [str(i+1)]
    print base_case
    result(trace, S, T)

# largest non-interleaving bond problem

# show whether the two nucleotide match
def rna_match(ch1, ch2):
    pair = set([ch1, ch2])
    if pair == set(['A', 'U']) or pair == set(['C', 'G']):
        return 1
    else:
        return 0

def rna_nibond(S):
    # initialization
    T = []
    for elem in S:
        T += [elem]
    m = len(S)

    trace = []
    table = []
    for i in range(0, m+1):
        table += [[]]
        trace += [[]]
        for j in range(0, m+1):
            table[i] += [0]
            trace[i] += [[]]

    # the loops
    for p in range(1, m):
        for i in range(1, m-p+1):
            j = i+p
            vm = table[i+1][j-1] + rna_match(S[i-1], T[j-1])
            best_match = 0
            for k in range(i, j):
                this_match = table[i][k] + table[k+1][j]
                if this_match > best_match:
                    best_match = this_match
            table[i][j] = max( best_match, vm )

    if debug:
        for entry in table:
            print entry
        print

    print 'Largest Set of Non-interleaving Bonds of RNA S'
    print 'S:', S
    print 'size:', table[1][m]




if __name__ == '__main__':
    """
    print example of global alignment from a trasparency shown in Tuesday\'s Lecture'
    S = 'agaagggcttccgcgacggcgacgttcgggggcctttttcttttgcggttt'
    T = 'agaagggcttccgcgacggcgacgttgagggggctcttttcttttgcggttt'
    print 'S:', S
    print 'T:', T
    g_align(S, T)

    print 'example of local alignment from lecture notes'
    S1 = 'abcxdex'
    T1 = 'xxxcde'
    print 'S:', S1
    print 'T:', T1
    l_align(S1, T1)
    """
    import sys
    args = sys.argv[1:]
    if len(args) != 6:
        print 'Usage: align <match> <mismatch> <gap open> <gap extend> <sequence1.txt> <sequence2.txt>'
    else:
        global g_match_score
        global g_mismatch_score
        global g_gap_open
        global g_gap_extend
        g_match_score = int(args[0])
        g_mismatch_score = int(args[1])
        g_gap_open = int(args[2])
        g_gap_extend = int(args[3])     

        #print g_match_score, g_mismatch_score, g_gap_open, g_gap_extend
        
        S = ''   
        S_file = open(args[4], 'rt')
        for line in S_file:
            S = S + line
        if S[-1] == '\n':
            S = S[:-1]
        
        T = ''
        T_file = open(args[5], 'rt')
        for line in T_file:
            T = T + line        
        if T[-1] == '\n':
            T = T[:-1]

        #print 'global alignment with affine gap penalty'
        #print 'S:', S
        #print 'T:', T
        g_align_gap(S, T)
    

    # S = 'GTAAATTGGCCAAGGTTAC'
    # T = 'TAGATA'
    #fit_align(S, T)
    #overlap_align('AAATTT', 'TATATA')
    #semiglobal_align('ACGTCAT','TCATGCA')

    #super_sequence('BLUE', 'ABLE')

    #inexact_repeat_align(T)
    #inexact_repeat_align(S)

    #is_infected('ATAGCTC','AAATAAAGGGGCCCCCTTTTTTTCC')
    #is_infected('ATTC', 'AAATTTTCC')

    #homo_dist('ACGCTA', 'AGGCCCGCTTTA')

    #get_lids('1273598', True)
    #get_lids('54123', False)

    #get_l2inc('821657439')

    #rna_nibond('GUGUACACG')
