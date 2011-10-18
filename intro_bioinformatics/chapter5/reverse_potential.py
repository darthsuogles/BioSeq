# genome reverse problem
# a greedy potential function solution

# calculate the potential of each permutation
def potential(seq):

    sum_pot = 0
    for i, elem in enumerate(seq):
        tmp =  abs(1 + i - elem)
        sum_pot += tmp+1 if tmp else 0
    return sum_pot


# reverse a sequence wrt the position given
def reverser(seq, i, j):

    rev = seq[i-1:j]
    rev.reverse()
    seq = seq[:i-1] + rev + seq[j:]

    return seq


# a greedy search method
def potential_rev(seq, depth=0):

    curr_pot = potential(seq)
    next_pot = curr_pot
    next_seq = ""

    if curr_pot == 0:
        print seq
        print depth+1, "steps in total\n"
        return depth + 2

    length = len(seq)
    for i in range(1, length):
        for j in range(i+1, length+1):
            this_seq = reverser(seq, i, j)
            this_pot = potential(this_seq)
            if this_pot < next_pot:
               next_pot = this_pot
               next_seq = this_seq

    print seq
    potential_rev(next_seq, depth+1)


# an a-star search method

# heuristic function
def astar_heuristic(seq):
    return num_bp(seq)[0]/2

def astar_terminate(seq):
    for i, elem in enumerate(seq):
        if elem != i+1:
            return False
    return True

def astar_str(seq):
    str_code = ''
    for elem in seq:
        str_code += str(elem)
    return str_code

def astar_rev(seq):
    heap = [[seq, astar_heuristic(seq), 0]]
    length = len(seq)
    visited = {astar_str(seq):True}

    while len(heap)>0:
        curr_seq = heap[0][0]
        last_dist = heap[0][1]
        last_step = heap[0][2]
        heap = heap[1:]

        for i in range(1, length):
            for j in range(i+1, length+1):
                this_seq = reverser(curr_seq, i, j)
                this_step = last_step + 1
                this_dist = astar_heuristic(this_seq) + this_step

                stred = astar_str(this_seq)
                if stred in visited:
                    continue
                visited[stred] = True

                if astar_terminate(this_seq):
                    print 'optimal solution takes', last_step+1, 'steps'
                    return

                flag = True
                for k, elem in enumerate(heap):
                    if elem[1] > this_dist:
                        flag = False
                        heap = heap[:k]+\
                               [[this_seq, this_dist, this_step]]\
                               +heap[k:]
                        break
                if flag:
                    heap += [[this_seq, this_dist, this_step]]


# breakpoint method in the textbook

# calculate the number of breakpoints in the sequence
circular = True # for sorting circular genome
def num_bp(seq):

    global circular

    num = 0
    bp_list = []
    curr_sign = ''

    for i in range(0, len(seq)-1):
        diff = seq[i]-seq[i+1]
        if circular:
            if seq[i] == 1 and seq[i+1] == len(seq):
                diff = 1
            elif seq[i] == len(seq) and seq[i+1] == 1:
                diff = -1

        if abs(diff) != 1:
            bp_list += [i+1]
            curr_sign = ''
            num += 1
            continue

        this_sign = '-' if diff>0 else '+'

        if curr_sign == '':
            curr_sign = this_sign

        if curr_sign != this_sign:
            bp_list += [i+1]
            curr_sign = ''
            num += 1

    if circular:
        if num == 0:
            for i in range(0, len(seq)-1):
                diff = seq[i] - seq[i+1]
                if abs(diff) != 1:
                    continue
                if diff > 0:
                    num = 1
                break

    else:
        if seq[0] != 1:
            num += 1
        if seq[-1] != len(seq):
            num += 1

    return [num, bp_list]


# search the result
def bp_rev(seq, depth = 0):

    print seq
    curr_bp = num_bp(seq)

    if curr_bp[0] == 0:
        print depth, 'iterations in total\n'
        return

    pos_list = [0]+curr_bp[1]+[len(seq)]

    # validate that we can revert this sequence
    flag = False
    for i in range(0, len(pos_list)-1):
        diff = seq[pos_list[i+1]-1] - seq[pos_list[i]]
        if diff >= 0:
            continue
        else:
            flag = True
            break

    if flag:
        # find the reversing position
        best_bp = curr_bp[0]
        best_seq = seq
        for i in range(0, len(seq)-1):
            for j in range(i+1, len(seq)):
                curr_seq = reverser(seq, i+1, j+1)
                curr_bp = num_bp(curr_seq)[0]
                if curr_bp < best_bp:
                    best_bp = curr_bp
                    best_seq = curr_seq

        #print best_seq
        bp_rev(best_seq, depth+1)

    # reverse an increasing strip
    else:
        for i in range(0, len(pos_list)-1):
            if pos_list[i+1] - pos_list[i] > 1:
                curr_seq = reverser(seq, pos_list[i]+1, pos_list[i+1])
                #print curr_seq
                bp_rev(curr_seq, depth+1)
                break



# breadth first search method to find the best reversal
def bfs_rev(seq):

    curr_bp = num_bp(seq)

    visited = [seq]
    parent = ['root']
    curr_pos = 0
    curr_seq = seq
    flag = True

    while flag:

        for i in range(1, len(curr_seq)):

            for j in range(i+1, len(curr_seq)+1):
                this_seq = reverser(curr_seq, i, j)
                #if this_seq not in visited: # the checking takes too long
                visited += [this_seq]
                parent += [curr_pos]
                if num_bp(this_seq)[0] == 0:
                    flag = False
                    break

            if not flag:
                break

        curr_pos += 1
        curr_seq = visited[curr_pos]


    cursor = len(visited)-1
    stack = []

    while cursor != 'root':
        stack = [visited[cursor]] + stack
        cursor = parent[cursor]

    depth = len(stack)-1
    while stack != []:
        print stack[0]
        stack = stack[1:]
    print 'optimal solution takes', depth, 'steps in total'


# swap sorting
def bfs_swap_sort(seq, target):
    curr_pot = potential(seq)
    visited = [seq]
    parent = [-1]
    curr_pos = 0

    flag = True
    while flag:
        curr_seq = visited[curr_pos]
        #print curr_seq, curr_pos, len(visited)

        for i in range(0, len(seq)-1):
            this_seq = curr_seq[:i]+[curr_seq[i+1]]\
                       +[curr_seq[i]]+curr_seq[i+2:]
            if this_seq in visited:
                continue
            visited += [this_seq]
            parent += [curr_pos]

            if this_seq == target:
                flag = False
                break
        curr_pos += 1

    stack = []
    curr_idx = len(parent) - 1
    while curr_idx != -1:
        stack = [visited[curr_idx]] + stack
        curr_idx = parent[curr_idx]

    depth = len(stack)-1
    while stack != []:
        print stack[0]
        stack = stack[1:]
    print 'solution takes', depth, 'steps'


# brute force multiple break point distance

# check the number of break points over all
# the sequences towards the base sequence
def mbp_check(seq_list, base):
    mapper = [0]*(len(base)+1)
    for i in range(0,len(base)):
        mapper[base[i]] = i+1

    total_bp = 0
    for seq in seq_list:
        last_elem = mapper[seq[0]]
        for elem in seq[1:]:
            mapped_elem = mapper[elem]
            if abs(mapped_elem - last_elem) != 1:
                total_bp += 1
            last_elem = mapped_elem
        if mapper[seq[0]] != 1:
                total_bp += 1
        if mapper[seq[-1]] != len(base):
                total_bp += 1

    return total_bp


# enumerate all permutation sequences and use each as the base
best_bp = -1
best_seq = ''
def bf_mbp_test(seq_list, seq, start=0):
    if start == len(seq):
        global best_bp, best_seq
        this_bp = mbp_check(seq_list, seq)
        if this_bp < best_bp:
            best_bp = this_bp
            best_seq = seq
        return

    for i in range(start, start+len(seq[start:])):
        bf_mbp_test(seq_list, seq[:start]+[seq[i]]+seq[start:i]+seq[i+1:],\
                    start+1)

# the integrated solution for brute force mbp
def bf_mbp(seq_list):
    global best_bp
    global best_seq
    best_bp = 999999999
    best_seq = ''

    bf_mbp_test(seq_list, range(1, len(seq_list[0])+1))

    print best_seq
    print 'optimal solution has', best_bp, 'break points\n'
    best_bp = 999999999
    best_seq = ''

# a greedy method for mbp problem
def greedy1_mbp(seq_list):
    n = len(seq_list[0])
    best_seq = range(1,n+1)
    best_dist = mbp_check(seq_list, best_seq)

    for seq in seq_list:
        this_dist = mbp_check(seq_list, seq)
        if this_dist < best_dist:
            best_dist = this_dist
            best_seq = seq

    print best_seq
    print "[greedy1] solution has a distance of", best_dist



# single source smallest reversal distance
def revdist_check(seq_list, seq):
    n = len(seq_list[0])
    num_seq = len(seq_list)
    #seq = range(1, n+1)
    is_found = [False]*num_seq
    num_found = 0

    visited = [seq]
    parent = [-1]
    sl_parent = [-1]*num_seq
    curr_pos = 0

    while True:
        curr_seq = visited[curr_pos]
        for i in range(1, n):
            for j in range(i, n+1):
                this_seq = reverser(curr_seq, i, j)
                if this_seq not in visited:
                    visited += [this_seq]
                    parent += [curr_pos]
                for k, seqor in enumerate(seq_list):
                    if seqor == this_seq and not is_found[k]:
                        is_found[k] = True
                        num_found += 1
                        sl_parent[k] = curr_pos

                    if num_found == num_seq:
                        ancestor = []
                        for elem in sl_parent:
                            tmp = elem
                            tmp_list = []
                            while tmp != -1:
                                tmp_list = [tmp] + tmp_list
                                tmp = parent[tmp]
                            ancestor += [[-1]+tmp_list]
                        idx = 0
                        last_common = 0
                        while True:
                            flag = False
                            if idx < len(ancestor[0]):
                                common_anc = ancestor[0][idx]
                                for anc_list in ancestor[1:]:
                                    if idx == len(anc_list) or \
                                       anc_list[idx] != common_anc:
                                        flag = True
                                        break
                            else:
                                flag = True
                            if flag:
                                total_depth = 0
                                for x in ancestor:
                                    total_depth += len(x)-1
                                #print visited[last_common]
                                #print 'solution has', total_depth,\
                                #      'reversions in total\n'
                                #break
                                return total_depth
                            else:
                                last_common = common_anc
                                idx += 1
                        return
        curr_pos += 1

# the brute force solution to minimun reverse distance problem
best_revdist = -1
best_revseq = ''
def bf_revdist_test(seq_list, seq, start=0):
    if start == len(seq):
        global best_revdist, best_revseq
        this_revdist = revdist_check(seq_list, seq)
        if this_revdist < best_revdist:
            best_revdist = this_revdist
            best_revseq = seq
        return

    for i in range(start, start+len(seq[start:])):
        bf_revdist_test(seq_list, seq[:start]+[seq[i]]+seq[start:i]+seq[i+1:],\
                    start+1)

# the integrated solution for brute force sol to reversion distance problem
def bf_revdist(seq_list):
    global best_revdist, best_revseq
    best_revdist = 999999999999999
    n = len(seq_list[0])
    bf_revdist_test(seq_list, range(1, n+1))
    print best_revseq
    print 'optimal solution has', best_revdist, 'in total'
    best_revdist = -1
    best_revseq = ''


if __name__ == "__main__":

    circular = False

    #potential_rev([3,4,6,5,8,1,7,2])
    #bp_rev([3,4,6,5,8,1,7,2])
    #astar_rev([3,4,6,5,8,1,7,2])

    #swap_sort([3,1,2,4])
    #bfs_swap_sort([3,4,6,5,8,1,7,2], [1,2,3,4,5,6,7,8])

    seq_list = [ [1,2,4,3,5,6],
                 [1,4,3,2,5,6],
                 [1,2,3,4,6,5] ]
    base = [1,2,3,4,5,6]
    print mbp_check(seq_list, [1,2,3,4,6,5])
    print mbp_check(seq_list, [1,4,3,2,5,6])
    print mbp_check(seq_list, base)

    '''seq_list = [ [1,2,3,5,7,6,4,8],
                 [2,4,3,5,8,7,1,6],
                 [4,5,6,7,3,1,8,2],
                 [2,6,8,5,4,3,7,1],
                 [3,4,5,6,8,7,1,2],
                 [8,7,5,4,3,1,2,6]
                 ]'''

    #bf_mbp(seq_list)
    #greedy1_mbp(seq_list)

    #bf_revdist(seq_list)
