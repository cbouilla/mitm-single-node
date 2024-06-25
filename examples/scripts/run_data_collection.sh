sudo-g5k apt install binutils-dev
sudo-g5k apt install nasm
sudo-g5k apt install libtbb-dev
sudo-g5k apt-get install libboost-stacktrace-dev
ml cmake/3.23.3_gcc-10.4.0

echo "Running parallel_run_experiments_collision$1.py"
python "parallel_run_experiments_collision$1.py"
