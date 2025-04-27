import polars as pl
import subprocess
import sys, os

path_program = "./build/vns_pure"
prefix_instance = "./instances/dr/U_{}/MD-VRBSP_U_{}_{}.txt"
prefix_result = "results/dr/U_{}"
seed = 1

shuffled_instances = (2 ** pl.int_range(3, 11, eager=True).shuffle(seed=seed)).to_list()

for L in shuffled_instances:
    output_path = prefix_result.format(L)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    shuffled_order = pl.int_range(16, 31, eager=True).shuffle(seed=seed).to_list()
    initial_objs = []
    best_objs = []
    for i in shuffled_order:
        print(f"Running for {L}_{i}")
        instance_path = prefix_instance.format(L, L, i)
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
        df.write_csv(f"{output_path}/progress{L}_{i}.csv")

    existing_df = None
    if os.path.isfile(f"{output_path}/consolidated{L}.csv"):
        existing_df = pl.read_csv(f"{output_path}/consolidated{L}.csv")

    master_df = pl.DataFrame({"ch": initial_objs, "vns": best_objs})
    master_df = master_df.with_columns(
        improvement=(pl.col("vns") - pl.col("ch")) / pl.col("ch")
    )

    if not existing_df.is_empty():
        master_df = pl.concat([existing_df, master_df])

    master_df.write_csv(f"{output_path}/consolidated{L}.csv")
