
# Configure the arguments of the attack in config.h

import argparse


def get_free_memory():
    """ Return available memory in kB"""
    import re
    with open("/proc/meminfo") as f:
        meminfo = f.read()
        matched = re.search(r"MemFree:\s+(\d+)", meminfo)
    print(int(matched.groups()[0]))

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


def parse_config(N=None,
                 NSERVERS=None,
                 NRECEIVERS=None,
                 RAM=None,
                 DIFFICULTY=None,
                 INTERVAL=None):

    """Edit N, L, NSERVERS, NRECEIVERS, RAM, and DIFFICULTY."""

    with open("include/config.h", "r") as f:
        lines = f.readlines()

    print(len(lines))
    for i in range(len(lines)):
        if lines[i][:10] == "#define N " and N:
            lines[i] = f"#define N {N}\n"
            print(f"Configured to attack {N} bytes!")

        if lines[i][:17] == "#define NSERVERS " and NSERVERS is not None:
            lines[i] = f"#define NSERVERS {NRECEIVERS}\n"
            print(f"Assumed there are {NSERVERS} servers!")

        if lines[i][:17] == "#define NSERVERS " and NRECEIVERS is not None:
            # we are using a stupid convention that the number of receivers
            # as the number of servers!
            lines[i] = f"#define NSERVERS {NRECEIVERS}\n"
            print(f"Assumed there are {NRECEIVERS} receivers!")

        if lines[i][:28] == "#define NRECEIVERS_PER_NODE ":
            if NRECEIVERS is None:
                nrecv_per_node = 1 # nservers = nreceivers

            if NSERVERS is not None and NRECEIVERS is not None:
                nrecv_per_node = int(NRECEIVERS//NSERVERS)

            if NSERVERS is None and NRECEIVERS is not None:
                raise ValueError("You need to defien --nservers if you are going to use --receivers")

            lines[i] = f"#define NRECEIVERS_PER_NODE {nrecv_per_node}\n"
            print(f"Assumed there are {nrecv_per_node} receivers/node!")

        if lines[i][:19] == "#define DIFFICULTY " and DIFFICULTY is not None:
            lines[i] = f"#define DIFFICULTY {DIFFICULTY}\n"
            print("difficulty done")

        if lines[i][:17] == "#define INTERVAL " and INTERVAL is not None:
            lines[i] = f"#define INTERVAL (1<<{INTERVAL}LL)\n"
            print(f"Interval between hashes has been updated to 2^{INTERVAL}")

        if lines[i][:17] == "#define TOTAL_RAM" and RAM is not None:
            lines[i] = f"#define TOTAL_RAM {RAM}LL\n"
            print(f"Updated RAM will used for dicts in a node  {RAM} bytes")

    with open("include/config.h", "w") as f:
        f.writelines(lines)


def run_phase_ii():
    import os

    os.system("ulimit -c unlimited")
    os.system("echo 'core' | sudo-g5k tee /proc/sys/kernel/core_pattern")

    os.system("cd lib/sha256_intel_avx/ && sudo-g5k apt install nasm && make clean && make all && cd ../../")

    os.system("make clean && make all")
    os.system("mpirun -machinefile $OAR_NODEFILE  -map-by node:PE=1 ./phase_ii")


if __name__ == "__main__":
    # parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("-N",
                        type=int,
                        help="Number of bytes to be attacked")

    parser.add_argument("--nservers",
                        type=int,
                        help="How many servers?")

    parser.add_argument("--receivers",
                        type=int,
                        help="How many receivers in total? (there can be more than one receiver per node)")

    parser.add_argument("--interval",
                        type=int,
                        help="for testing purpose only! deliberately change the interval even if it is not correct")

    parser.add_argument("-r", "--ram",
                        type=int,
                        help="what is the available memory for the dictionary per server?")

    parser.add_argument("-d", "--difficulty",
                        type=int,
                        help="How many bits are zero")

    args = parser.parse_args()


# def parse_config(N=None,
#                  NSERVERS=None,
#                  NRECEIVERS=None,
#                  RAM=None,
#                  DIFFICULTY=None,
#                  INTERVAL=None):
    parse_config(N=args.N,
                 NSERVERS=args.nservers,
                 NRECEIVERS=args.receivers,
                 RAM=args.ram,
                 DIFFICULTY=args.difficulty,
                 INTERVAL=args.interval)

    # run phase_iii
    run_phase_ii()
