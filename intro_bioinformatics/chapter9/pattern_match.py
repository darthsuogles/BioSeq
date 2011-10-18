# match a pattern given a text

# naive matching algorithm
def simple_match(pattern, text):
    l = len(pattern)
    n = len(text)

    for i in range(0, n-l+1):
        flag = True
        for k in range(0, l):
            if pattern[k] != text[i+k]:
                flag = False
                break
        if flag:
            return True
    return False


# slightly advanced matching
# this will tolerate k mismatches
def better_match(pattern, text, k):
    l = len(pattern)
    n = len(text)

    for i in range(0, n-l+1):
        mis_count = 0
        for j in range(0, l):
            if pattern[j] != text[i+j]:
                mis_count += 1
        if mis_count <= k:
            return True
    return False


# a better algorithm that tolerate
# insertions and deletions
def tolerant_match(pattern, text, k):
    l = len(pattern)
    n = len(text)
    table = []
    for i in range(0, n+1):
        table += [[k]]
        for j in range(1, l+1):
            if i == 0:
                table[i] += [k-j]
            table[i] += [0]

    for i in range(1, n+1):
        for j in range(1, l+1):
            nl = table[i][j-1] - 1
            nu = table[i-1][j] - 1
            nd = table[i-1][j-1] + ( 0 if pattern[j-1] == text[i-1] else -1 )

            table[i][j] = max(nl, max(nu, nd))

    for i in range(1, n+1):
        if table[i][l] >= 0:
            return True
    return False




if __name__ == "__main__":
    print "simple match"
    print simple_match("AA", "ATTAA")
    print simple_match("AA", "ATTTA")
    print

    print "better than simple match"
    print better_match("AA", "ATTAA", 0)
    print better_match("AA", "ATTTA", 1)
    print

    print "tolerant match"
    print tolerant_match("AA", "ATTT", 0)
    print tolerant_match("AA", "ATTT", 1)
    print tolerant_match("AA", "A", 1)
    print tolerant_match("TT", "AAAA", 2)
