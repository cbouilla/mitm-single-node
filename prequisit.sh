#!/bin/bash
sudo-g5k apt install -y libtbb-dev
ulimit -S -c 0
sudo-g5k sysctl -w kernel.core_pattern=/coredumps/core-%e-%s-%u-%g-%p-%t
