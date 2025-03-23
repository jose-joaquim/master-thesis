import polars as pl
import subprocess
import pdb

path_program = "./build/vns_pure"
prefix_instance = "./instances/dr/U_{}/MD-VRBSP_U_{}_{}.txt"
prefix_result = "results/dr/U_{}"
seed = 1

shuffled_instances = (2 ** pl.int_range(3, 11, eager=True).shuffle(seed=seed)).to_list()

for L in shuffled_instances:
    output_path = prefix_result.format(L)
    instance_path = prefix_instance.format(L, L, i)
    shuffled_order = pl.int_range(1, 3, eager=True).shuffle(seed=seed).to_list()
    initial_objs = []
    best_objs = []
    for i in shuffled_order:
        run = subprocess.run(
            [path_program, instance_path],
            capture_output=True,
            text=True,
            cwd="./",
        )

        output = [x.split(",") for x in run.stdout.split("\n")]
        obj = [float(x[0]) for x in output[:-1]]
        at = [float(x[1]) for x in output[:-1]]

        assert len(obj) > 0

        initial_objs.append(obj[0])
        best_objs.append(obj[-1])
        df = pl.DataFrame({"objective": obj, "at": at})
        df.write_csv(f"{prefix_result}/progress{L}_{i}.csv")

    master_df = pl.DataFrame({"ch": initial_objs, "vns": best_objs})
    master_df.write_csv(f"{prefix_result}/consolidated{L}.csv")
