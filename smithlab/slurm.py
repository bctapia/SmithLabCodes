"""
Copyright 2025 Brandon C. Tapia

MIT License
"""


def write_batch(
    batch_out, N, n, partition, cpus_per_task=None, mem_per_cpu=None, out_file="output_%j.txt", command=None):

    with open(batch_out, "w", encoding="utf-8") as f:
        f.write(
            f"#!/bin/bash\n#SBATCH -N {N}\n#SBATCH -n {n}\n#SBATCH --partition={partition}\n"
        )
        if cpus_per_task is not None:
            f.write(f"#SBATCH --cpus_per_task={cpus_per_task}\n")
        if out_file is not None:
            f.write(f"#SBATCH -o {out_file}\n")
        if command:
            f.write(command)
