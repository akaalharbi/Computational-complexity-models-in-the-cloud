# Configure the arguments of the attack in config.h

import argparse


def get_free_memory():
    """ Return available memory in kB"""
    import re
    with open("/proc/meminfo") as f:
        meminfo = f.read()
        matched = re.search(r"MemFree:\s+(\d+)", meminfo)

    return int(matched.groups()[0])


def get_ram_size():
    """Return ram size in kB."""
    import re
    with open("/proc/meminfo") as f:
        meminfo = f.read()
        matched = re.search(r"MemTotal:\s+(\d+)", meminfo)

    return 1000*int(matched.groups()[0])


def compute_dict_size():
    """Compute how many slots in the dictionary this device should get."""
    from math import log2
    ram_size = get_ram_size()
    leave_free_memory = (10**9) * 7  # 6 GB

    nslots = (ram_size - leave_free_memory)/4
    print(f"This server will hold {nslots}=2^{log2(nslots)} elements")

    return int(nslots)


def parse_config(N=None, L=None, NSERVERS=None, DIFFICULTY=None):
    """Edit N, L, NSERVERS, and DIFFICULTY."""
    print(f"N={N}, L={L}, NSERVERS={NSERVERS}, DIFFICUTLY={DIFFICULTY}")

    with open("include/config.h", "r") as f:
        lines = f.readlines()

    print(len(lines))
    for i in range(len(lines)):
        if lines[i][:10] == "#define L " and L is not None:
            lines[i] = f"#define L {L}\n"
            print("L done")

        if lines[i][:10] == "#define N " and N:
            lines[i] = f"#define N {N}\n"
            print("N done")

        if lines[i][:17] == "#define NSERVERS " and NSERVERS is not None:
            lines[i] = f"#define NSERVERS {NSERVERS}\n"
            print("nservers done")

        if lines[i][:19] == "#define DIFFICULTY " and DIFFICULTY is not None:
            lines[i] = f"#define DIFFICULTY {DIFFICULTY}\n"
            print("difficulty done")

        if lines[i][:17] == "#define TOTAL_RAM":
            memory_left = 7*10**9
            lines[i] = f"#define TOTAL_RAM {get_ram_size() - memory_left}LL\n"
            print(f"Updated RAM will used = {get_ram_size() - memory_left}kb")

    with open("include/config.h", "w") as f:
        f.writelines(lines)




if __name__ == "__main__":
    # parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("-N",
                        type=int,
                        help="Number of bytes to be attacked")

    parser.add_argument("-L",
                        type=int,
                        help="Number of distingusihed hashes to search among them")

    parser.add_argument("-s", "--NSERVERS",
                        type=int,
                        help="How many server? (i.e. receivers)")

    parser.add_argument("-d", "--DIFFICULTY",
                        type=int,
                        help="How many bits are zero")

    args = parser.parse_args()

    parse_config(N=args.N,
                 L=args.L,
                 NSERVERS=args.NSERVERS,
                 DIFFICULTY=args.DIFFICULTY)
