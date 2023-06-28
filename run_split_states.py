"""
Allow coredumps, install nasm, compile the program, and run split_states.
"""


def run_phase_ii():
    """`mpirun` command is adapted for Grid5000.fr, grenbole, dahu."""
    import os

    os.system("ulimit -c unlimited")
    os.system("echo 'core' | sudo-g5k tee /proc/sys/kernel/core_pattern")
    os.system("sudo-g5k apt install nasm")
    os.system("cd lib/sha256_intel_avx/ && make clean && make all && cd ../../")

    os.system("make clean && make all")
    os.system("mpirun -machinefile $OAR_NODEFILE  -mca mtl psm2 -mca pml ^ucx,ofi -mca btl ^ofi,openib -map-by node:PE=1 ./split_states")


if __name__ == "__main__":
    run_phase_ii()
