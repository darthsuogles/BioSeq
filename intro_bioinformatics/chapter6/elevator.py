# a simple program to solve the elevator problem 

def push(initial, target):
    val = [7e14]*51
    val[initial] = 0
    modify = True
    while modify:
        modify = False
        for n in range(1, 51):
            val1 = 7e14
            val2 = 7e14
            if n > 11 and val[n-11]+11 < val[n]:
                modify = True
                val[n] = val[n-11]+11
            if n < 45 and val[n+6]+6 < val[n]:
                modify = True
                val[n] = val[n+6]+6
    print val[target]

if __name__ == "__main__":

    push(32,33)
    push(11,22)
    push(11,5)
    push(50,49)
