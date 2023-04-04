# extrae depends on spack openmpi installation
module load openmpi/4.1.3_gcc-10.2.0
# install binutils and libxml2 for extrae profiling
sudo-g5k apt install libxml2-dev
sudo-g5k apt install binutils-dev
sudo-g5k apt install libunwind-dev
sudo-g5k apt-get install libtbb-dev
# dyninst prerequisite

echo "deb http://deb.debian.org/debian bullseye-backports main" | sudo-g5k tee /etc/apt/sources.list.d/backports.list
sudo-g5k update
sudo-g5k apt-get install -t bullseye-backports libelf-dev -y
sudo-g5k apt-get install -t bullseye-backports libdwarf-dev -y
sudo-g5k apt-get install build-essential cmake libboost-all-dev libelf-dev libiberty-dev -y



# core dump in the local directory if an error occur
ulimit -c unlimited
echo 'core' | sudo-g5k tee /proc/sys/kernel/core_pattern

# compile
cd lib/sha256_intel_avx/
sudo-g5k apt install nasm
make clean && make all
cd ../../
make clean && make all

# mpirun -machinefile $OAR_NODEFILE ./verify_data
mpirun -machinefile $OAR_NODEFILE  -map-by node:PE=1 ./phase_ii

