import polars as pl
import subprocess
import os

vns_program_path = "./build/vns_pure"
brkga_program_path = "./codes-wo-cmake/decoder/bin/decoder_brkga"
prefix_instance = "./instances/dr/U_{}/MD-VRBSP_U_{}_{}.txt"
prefix_result = "results/dr/U_{}"
seed = 1

# Run first

shuffled_instances = (2 ** pl.int_range(3, 12, eager=True).shuffle(seed=seed)).to_list()
shuffled_order = pl.int_range(1, 31, eager=True).shuffle(seed=seed).to_list()

instances_backup = [1024, 2048]

order_backup = [
    25,
    7,
    28,
    19,
    1,
    # 5,
    # 29,
    # 15,
    # 21,
    # 30,
    # 14,
    # 9,
    # 16,
    # 23,
    # 24,
    # 11,
    # 18,
    # 4,
    # 12,
    # 27,
    # 26,
    # 20,
    # 17,
    # 8,
    # 3,
    # 13,
    # 10,
    # 6,
    # 2,
    # 22,
]

# for L in shuffled_instances:
for L in instances_backup:
    output_path = prefix_result.format(L)
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    ch_objs = []
    vns_objs = []
    brkga_objs = []
    instance = []
    for i in order_backup:
        print(f"Running for {L}_{i}")
        instance_path = prefix_instance.format(L, L, i)
        ch_json = f"{output_path}/solutionCH{L}_{i}.json"
        vns_json = f"{output_path}/solutionVNS{L}_{i}.json"
        brkga_json = f"{output_path}/solutionBRKGA{L}_{i}.json"

        run = subprocess.run(
            [vns_program_path, instance_path, ch_json, vns_json],
            capture_output=True,
            text=True,
            cwd="./",
        )

        output = [x.split(",") for x in run.stdout.split("\n")]
        obj = [float(x[0]) for x in output[:-1]]
        at = [float(x[1]) for x in output[:-1]]

        assert len(obj) > 0

        instance.append(i)
        ch_objs.append(obj[0])
        vns_objs.append(obj[-1])
        df = pl.DataFrame({"objective": obj, "at": at})
        df.write_csv(f"{output_path}/progress{L}_{i}.csv")

        run = subprocess.run(
            [brkga_program_path, instance_path, brkga_json],
            capture_output=True,
            text=True,
            cwd="./",
        )

        brkga_objs.append(-1 * float(run.stdout))

    master_df = pl.DataFrame(
        {"instance": instance, "ch": ch_objs, "vns": vns_objs, "brkga": brkga_objs}
    )
    master_df = master_df.with_columns(
        vns_over_ch=(pl.col("vns") - pl.col("ch")) / pl.col("ch"),
        vns_over_brkga=(pl.col("vns") - pl.col("brkga")) / pl.col("brkga"),
    )

    existing_df = None
    if os.path.isfile(f"{output_path}/consolidated{L}.csv"):
        existing_df = pl.read_csv(f"{output_path}/consolidated{L}.csv")
        master_df = pl.concat([existing_df, master_df])

    master_df.write_csv(f"{output_path}/consolidated{L}.csv")
