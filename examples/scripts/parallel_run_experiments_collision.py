"""
Edits `demos/sha2_claw_demo.cpp` macros specifically, NBYTES_A, NBYTES_B,
NBYTES_C, for values in the range (1, 4). For each triple, run the program
100 times (since it relies on random inputs).
"""

import os
# import itertools as itr
# from tqdm import tqdm
import subprocess
import time
import math
import random
import argparse


# Read the problem size from the command line
parser = argparse.ArgumentParser()

parser.add_argument("-n",
                    type=int,
                    help="f:{0, 1}^n -> {0, 1}^n, run many collision finding experiments for this function")

args = parser.parse_args()




#bits_range = list(range(18, 19))
bits_range = list(range(args.n, args.n+1))

# let's focus when they are equal
all_triples = [(i, i, i) for i in bits_range]  # itr.product(bytes_range, repeat=3)
print(f"all triples= {all_triples}")
nruns = 20000  # How many times we run the code for the same triple value
difficulty_range = 48  # i.e. difficulty between 0 and difficulty_range included



def is_memory_enough(command):
    import psutil
    import re

    # Get the available memory at the system
    memory = psutil.virtual_memory()
    available_mem = memory.available
    
    # see how many bytes the problem requires
    pattern = r'./sha2_collision_demo_(\d+)_(\d+) (\d+) (\d+)'
    results = re.search(pattern, command)

    if results:
        nbits = int(results.group(1))
        log2_ram = int(results.group(3))

        print(f"command: {command}")


        # memory estimated to be used by the program
        space = (2**log2_ram) * (8 + 8 + (nbits_C+7)//8 )
        print(f"nbits={nbits}, log2_ram={log2_ram}, space={space} bytes, available={available_mem} bytes")
        if 2**(log2_ram+2 - 3) < 0.85 * available_mem:
            return True

    return False
    

# source: chatgpt
def run_commands_in_parallel(commands, max_processes, timeout=8*60):
    """Run shell commands in parallel without exceeding the max number of processes, with a timeout for each."""
    processes = []
    start_times = {}

    # Define timeout in seconds (10 minutes * 60 seconds/minute)


    # Helper function to clean up finished processes and handle timeouts
    def clean_up_processes():
        nonlocal processes
        current_time = time.time()
        finished_processes = []
        for process in processes:
            if process.poll() is not None:  # None means the process hasn't finished yet
                finished_processes.append(process)
            elif (current_time - start_times[process]) > timeout:
                process.terminate()  # Terminate processes that exceed the timeout
                finished_processes.append(process)
        processes = [p for p in processes if p not in finished_processes]

    for command in commands:
        # Clean up finished processes
        clean_up_processes()

        # Wait if we have reached the maximum number of parallel processes
        not_enough_memory = not (is_memory_enough(command))
        print(f"number of processes now = {len(processes)}")
        while (len(processes) >= max_processes) or not_enough_memory:
            time.sleep(1)  # Sleep briefly to avoid constant CPU usage
            clean_up_processes()

        # Start a new subprocess
        print(f"running {command}")
        process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        processes.append(process)
        start_times[process] = time.time()  # Track the start time of this process

    # Wait for all processes to complete
    while processes:
        time.sleep(1)  # Sleep to reduce CPU usage
        clean_up_processes()

    print("All processes have completed.")



commands = []
# create a list commands to be runned
for triple in all_triples:
    nbits_C = triple[-1]
    #difficulty_range = (nbits_C // 2) + 3
    #for difficulty in range(min(difficulty_range + 1, nbits_C//2)):
    #for difficulty in range(difficulty_range):
    for log2_ram in range(4,  (nbits_C//2) +  4):
        # theta = alpha * sqrt(w/n), difficulty = 1/theta, alpha = 2.25,
        # w = 2^log2_ram, n = 2^nbits_C
        difficulty = int((nbits_C - log2_ram)//2 + math.log2(2.225))
        run_cmd = f"./sha2_collision_demo_{nbits_C}_{difficulty} {log2_ram} {difficulty}"
        for _ in range(nruns):
            commands.append(run_cmd)

commands[::-1]
random.shuffle(commands)

# Running in Parallel
timeout = 3*60*60
max_processes =  os.cpu_count() 
print(f"Going to use {max_processes} parallel processes...")
run_commands_in_parallel(commands, max_processes, timeout)


print("DONE!")
