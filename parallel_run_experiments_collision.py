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

bits_range = list(range(8, 10))
# let's focus when they are equal
all_triples = [(i, i, i) for i in bits_range]  # itr.product(bytes_range, repeat=3)
nruns = 30  # How many times we run the code for the same triple value
difficulty_range = 48  # i.e. difficulty between 0 and difficulty_range included


def edit_claw_demo(triple):
    """Replace the nbytes_c, nbytes_b, nbytes_a from."""
    keywords = ["NBITS_A", "NBITS_B", "NBITS_C"]

    with open("demos/sha2_collision_demo.cpp", "r") as claw,\
         open("demos/tmp_col.cpp", "w") as tmp:
        for line in claw:
            for i in range(len(keywords)):
                if line.startswith("#define " + keywords[i]):
                    line = "#define " + keywords[i] + " " + str(triple[i])\
                         + " \n"
            tmp.write(line)

    os.replace("demos/tmp_col.cpp", "demos/sha2_collision_demo.cpp")


def print_errors_if_any(result):
    """Print errors if any from subprocess run."""
    if result.stderr:
        print("Error:", result.stderr.decode())


def compile_project():
    """Compile the project."""
    clean_cmake = "rm -rf CMakeCache.txt  CMakeFiles/"
    compile_cmd = "cmake -S . -B . -D AES_IMPL=aesni && make"

    result = subprocess.run(clean_cmake, shell=True,
                            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    print_errors_if_any(result)

    result = subprocess.run(compile_cmd, shell=True,
                            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    print_errors_if_any(result)



def rename_executable(nbits_A_C, difficulty):
    move_cmd = f"cp sha2_collision_demo sha2_collision_demo_{nbits_A_C}_{difficulty}"
    result = subprocess.run(move_cmd, shell=True,
                            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    print_errors_if_any(result)



def run_project(log2_nbytes, difficulty, timeout=300):
    """Run the code."""
    run_cmd = f"./sha2_collision_demo {log2_nbytes} {difficulty}"
    try:
        result = subprocess.run(run_cmd,
                                shell=True,
                                timeout=timeout,
                                stdout=subprocess.DEVNULL,
                                stderr=subprocess.PIPE)

    except subprocess.TimeoutExpired:
        print("Command timed out after 1 hour.")
    else:
        print_errors_if_any(result)


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
        while len(processes) >= max_processes:
            time.sleep(1)  # Sleep briefly to avoid constant CPU usage
            clean_up_processes()

        # Start a new subprocess
        process = subprocess.Popen(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        processes.append(process)
        start_times[process] = time.time()  # Track the start time of this process

    # Wait for all processes to complete
    while processes:
        time.sleep(1)  # Sleep to reduce CPU usage
        clean_up_processes()

    print("All processes have completed.")


# First compile all the triples sequentially
for triple in all_triples:
    print(f"Triple: {triple}")
    edit_claw_demo(triple)
    compile_project()

    nbits_C = triple[-1]
    for difficulty in range(min(difficulty_range + 1, nbits_C)):
        rename_executable(nbits_C, difficulty)

    # for difficulty in range(min(difficulty_range + 1, nbits_C)):
    #     print(f"difficulty={difficulty}, (|C|, |A|, |B|)={triple}")
    #     try:
    #         for _ in tqdm(range(nruns)):
    #             run_project(int(nbits_C/2), difficulty)
    #     except:
    #         print("Think about installing tqdm by pip install tqdm")
    #         for _ in range(nruns):
    #             run_project(difficulty)

commands = []
repaeat = 40

for triple in all_triples:
    nbits_C = triple[-1]
    for difficulty in range(min(difficulty_range + 1, nbits_C//2)):
        for log2_ram in range(nbits_C//2, nbits_C):
            run_cmd = f"./sha2_collision_demo_{nbits_C}_{difficulty} {log2_ram} {difficulty}"
            commands.append(run_cmd)


# Running in Parallel
timeout = 8*60
max_processes = os.cpu_count()
print(f"Going to use {max_processes} parallel processes...")
run_commands_in_parallel(commands, max_processes, timeout)


# Finally collect all data in a single file
os.system("cat data/collision_* > data/collision_summary.csv")
print("DONE!")
