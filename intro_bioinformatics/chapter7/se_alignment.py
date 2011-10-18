# space efficient alignment

# matching score for two nucleotides
def score(ch1, ch2):
    if ch1 == '-' or ch2 == '-':
        return -1
    if ch1 == ch2:
        return 1
    else:
        return -1

# linear-space global alignment algorithm
# each column has a return value
# local path constructed with linear space requirement
def se_g_align(S, T):
    infty = 7e21
    m = len(S)
    n = len(T)

    # base case
    if m == 0 or n == 0:
        return []

    # initialize the table
    table = []
    for i in range(0, m+2):
        table += [[]]
        for j in range(0, n+2):
            if i == 0 or j == 0 or i == m+1 or j == n+1:
                table[i] += [-infty]
            else:
                table[i] += [0]
    table[0][0] = 0
    table[-1][-1] = 0

    # the middle column
    pivot = n/2 if n > 1 else 1

    # calculate the left half
    for i in range(1, m+1):
        for j in range(1, pivot+1):
            vl = table[i][j-1]
            vd = table[i-1][j]
            vu = table[i-1][j-1] + score(S[i-1], T[j-1])
            table[i][j] = max( vl, max(vd, vu) )

    left_tmp = []
    for i in range(1, m+1):
        left_tmp += [table[i][pivot]]

    # calculate the right half
    rev_i_range = range(1, m+1)
    rev_i_range.reverse()
    rev_j_range = range(pivot, n+1)
    rev_j_range.reverse()

    for i in rev_i_range:
        for j in rev_j_range:
            vl = table[i][j+1]
            vu = table[i+1][j]
            vd = table[i+1][j+1] + score(S[i-1], T[j-1])
            table[i][j] = max( vl, max(vu, vd) )

    # update the pivot column
    for i in range(1, m+1):
        table[i][pivot] += left_tmp[i-1]

    # find the best value in the pivot column
    best = -infty
    best_i = 0
    best_j = pivot
    for k in range(1, m+1):
        curr = table[k][pivot]
        if curr > best:
            best = curr
            best_i = k

    if False:
        for entry in table:
            print entry
        print

    # divide and conquer
    if n > 1:
        left = se_g_align(S[:best_i-1], T[:best_j-1])
        right = se_g_align(S[best_i:], T[best_j:])
        for i, elem in enumerate(right):
            right[i] = [elem[0]+best_i, elem[1]+best_j]

        return left + [[best_i, best_j]] + right

    # another base case
    else:
        return [[best_i, best_j]]



if __name__ == "__main__":

    print se_g_align('TACG', 'GTACG')
