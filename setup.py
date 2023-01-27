import argparse


def get_free_memory():
    """ Return available memory in kB"""
    import re
    with open("/proc/meminfo") as f:
        meminfo = f.read()
        matched = re.search(r"MemFree:\s+(\d+)", meminfo)

    return int(matched.groups()[0])


def get_ram_size():
    """ Return ram size in kB """
    import re
    with open("/proc/meminfo") as f:
        meminfo = f.read()
        matched = re.search(r"MemTotal:\s+(\d+)", meminfo)

    return int(matched.groups()[0])


available_mem = get_free_memory()


def parse_config(N=None, L=None, NSERVERS=None, DIFFICULTY=None):
    """Edit N, L, NSERVERS, and DIFFICULTY."""
    print(N, L, NSERVERS, DIFFICULTY)

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

    with open("include/config.h", "w") as f:
        f.writelines(lines)


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





print(args)

parse_config(N=args.N,
             L=args.L,
             NSERVERS=args.NSERVERS,
             DIFFICULTY=args.DIFFICULTY)
