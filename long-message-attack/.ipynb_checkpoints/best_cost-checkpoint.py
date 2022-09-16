def cost_n_l(n, l):
    """
    Given
      n: number of bits of the compression functtion
      l: 2^l is the number blocks
    """
    import os
    
    os.system(f"perf stat -o stats.txt -a -r 2 -e power/energy-pkg/ ./long_message_attack {n} {l}") #"perf stat -a -e power/energy-pkg/ ./long_message_attack 30 15"

    # open file
    f = open("stats.txt", "r")
    lines = f.readlines()
    
    idx = -1
    for i in range(len(lines)):
        if "Joules" in lines[i]:
            idx = i
            break
    
    # clean line  
    line = lines[idx]
    print(line)
    line = line.strip()
    idx2 = line.index('g/')
    line = line[:idx2+1]
    idx_joule = line.index("J")
    joules = float(line[:idx_joule-1])
    
    return joules

def min_cost(n):
    cost = cost_n_l(2, 1)
    l = 1
    
    for i in range(1, n):
        tmp_cost = cost_n_l(n, i)
        if tmp_cost < cost:
            cost = tmp_cost
            l = i
    
    return l, cost



if __name__ == "__main__":
    print("For each n, we try to find the best l the minimizes the energy consumption")
    
    for n in range(16, 40):
        l, cost = min_cost(n)
        print(f"n={n}, l={l}, cost={cost} Joules")
