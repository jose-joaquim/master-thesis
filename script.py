import polars as pl
import subprocess
import pdb

path_program = "./build/vns_pure"
prefix_instance = "./instances/du/U_{}/MD-VRBSP_U_{}_{}.txt"
seed = 1

shuffled_instances = (2 ** pl.int_range(3, 11, eager=True).shuffle(seed=seed)).to_list()

for L in shuffled_instances:
    shuffled_order = pl.int_range(1, 31, eager=True).shuffle(seed=seed).to_list()
    for i in shuffled_order:
        run = subprocess.run(
            [path_program, prefix_instance.format(L, L, i)],
            capture_output=True,
            text=True,
            cwd="./",
        )

        output = [x.split(",") for x in run.stdout.split("\n")]
        obj = [float(x[0]) for x in output[:-1]]
        at = [float(x[1]) for x in output[:-1]]

        df = pl.DataFrame({"objective": obj, "at": at})
        df.write_csv(f"progress{L}_{i}.csv")
