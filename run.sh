make BUILD_LIBS
make purge
make all

./phase_i
cp data/send/digests/* data/receive/digests/
mpirun --oversubscribe -np 20 ./phase_ii
cat data/send/messages/* > data/receive/messages/archive
./phase_iii
