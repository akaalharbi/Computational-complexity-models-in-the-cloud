cd lib/sha256_intel_avx/
sudo-g5k apt install nasm
make clean && make all
cd ../../
make clean && make all
mpirun -machinefile $OAR_NODEFILE ./verify_data
