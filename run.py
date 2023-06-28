"""Configure the arguments of the attack in config.h
python run.py  --nservers 32 --receivers 32 -N 12 --ram 64000000000 --interval 30
"""
import argparse


def get_free_memory():
    """Return available memory in kB."""
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
    """Allow coredumps, compile all files, then run them using mpirun.

    Note: this will be rewritten using python subprocess to allow using
    perf.
    """
    import os

    os.system("ulimit -c unlimited")
    # @todo: remove the following two lines 
    os.system("echo 'core' | sudo-g5k tee /proc/sys/kernel/core_pattern")
    os.system("sudo-g5k apt install nasm")
    os.system("cd lib/sha256_intel_avx/ && make clean && make all && cd ../../")
    os.system("make clean && make all")
    # here is
    os.system("mpirun -machinefile $OAR_NODEFILE  -map-by node:PE=1 ./phase_ii")


def clean_hostfile():
    """Remove all repeated hosts names in $OAR_NODEFILE

    The goal is to run one process per node.
    """
    import os
    file_path = os.environ["OAR_NODEFILE"]
    lines_set = set([])
    with open(file_path, 'r') as f:
        for line in f:
            lines_set.add(line)
    lines = list(lines_set)
    with open("data/tmp_hosts", 'w') as f:
        for line in lines:
            f.write(line)


def run_perf():
    """Record the energy consumption."""
    import subprocess

    # make the python wrapper for perf an executable
    subprocess.run("chmod +x src/perf_bench.py", shell=True)
    # we know it's long commnad but what can we do!
    # cmd = 'mpirun -machinefile data/tmp -mca mtl psm2 -mca pml ^ucx,ofi -mca btl ^ofi,openib -map-by node:PE=1 "python src/perf_bench.py"'
    cmd = ['mpirun',
           '-machinefile',
           'data/tmp',
           '-mca',
           'mtl',
           'psm2',
           '-mca',
           'pml',
           '^ucx,ofi',
           '-mca',
           'btl',
           '^ofi,openib',
           '-map-by',
           'node:PE=1',
           'python src/perf_bench.py']

    subprocess.run(cmd, shell=True)


if __name__ == "__main__":
    from multiprocessing import Process
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
    p1 = Process(target=run_phase_ii)
    p2 = Process(target=run_perf)

    p1.start()
    p2.start()

    # todo use timeout to exit the process when collecting enough candidates
    p1.join()
    p2.join()
