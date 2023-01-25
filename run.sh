make purge
make all

./phase_i
mpirun --oversubscribe -np 20 ./phase_ii
cat data/messages/* > data/messages/archive
./phase_iii
