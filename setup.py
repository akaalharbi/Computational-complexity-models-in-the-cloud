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
    import re



    with open("include/config.h", "r") as f:
        lines = f.readlines()


    for i in range(len(lines)):
        if lines[i][2:] == "L " and L:
            lines[i] = f"L {L}"

        if lines[i][2:] == "N " and N:
            lines[i] = f"N {N}"

        if lines[i][9:] == "NSERVERS " and NSERVERS:
            lines[i] = f"NSERVERS {L}"

        if lines[i][11:] == "DIFFICULTY " and DIFFICULTY:
            lines[i] = f"DIFFICULTY {L}"

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






parse_config(args.N, args.L, args.NSERVERS, args.DIFFICULTY)
