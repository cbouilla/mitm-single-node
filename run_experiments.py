"""
Edits `demos/sha2_claw_demo.cpp` macros specifically, NBYTES_A, NBYTES_B,
NBYTES_C, for values in the range (1, 4). For each triple, run the program
100 times (since it relies on random inputs).
"""

import os
import itertools as itr
from tqdm import tqdm
import subprocess

bytes_range = list(range(1, 5))
all_triples = itr.product(bytes_range, repeat=3)
nruns = 100  # How many times we run the code for the same triple value


def edit_claw_demo(triple):
    """Replace the nbytes_c, nbytes_b, nbytes_a from."""
    keywords = ["NBYTES_A", "NBYTES_A", "NBYTES_C"]

    with open("demos/sha2_claw_demo.cpp", "r") as claw,\
         open("demos/tmp.cpp", "w") as tmp:
        for line in claw:
            for i in range(len(keywords)):
                if line.startswith("#define " + keywords[i]):
                    line = "#define " + keywords[i] + " " + str(triple[i])\
                         + " \n"
            tmp.write(line)

    os.replace("demos/tmp.cpp", "demos/sha2_claw_demo.cpp")


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


def run_project():
    """Run the code."""
    run_cmd = "./sha2_claw_demo"
    result = subprocess.run(run_cmd, shell=True,
                            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    print_errors_if_any(result)


try:  # if tqdm was installed.
    for triple in tqdm(all_triples):
        edit_claw_demo(triple)
        compile_project()
        for _ in tqdm(range(nruns)):
            run_project()

except ImportError:
    print("Consider installing `tqdm` for a nicer visuals "
          "via  pip install tqdm")
    for triple in tqdm(all_triples):
        edit_claw_demo(triple)
        compile_project()
        for _ in tqdm(range(nruns)):
            run_project()
