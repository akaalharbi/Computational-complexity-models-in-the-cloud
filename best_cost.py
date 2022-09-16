def cost_n_l(n, l):
    """
    Given
      n: number of bits of the compression functtion
      l: 2^l is the number blocks
    """
    import os
    command = f"perf stat -o stats.txt -a -r 2 -e power/energy-pkg/ ./long_message_attack {n} {l} > tmp"
    os.system(command)

    # open file
    with open("stats.txt", "r") as f:
        lines = f.readlines()

        idx = -1
        for i in range(len(lines)):
            if "Joules" in lines[i]:
                idx = i
                break

        # clean line  
        line = lines[idx]
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

    with open(f"satistics/all_stats.txt", "w") as f:
        for n in range(16, 54):
            
            f.write(f"n={n}: cost(n, l) is the lth entry")
            f.write("[inf")
            for l in range(1, min(n, 24)):
                joules = cost_n_l(n, l)
                f.write(f", {joules}")
            f.write("]\n")

        




