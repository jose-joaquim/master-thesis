import polars as pl
import json

seed = 1
shuffled_instances = (2 ** pl.int_range(3, 12, eager=True).shuffle(seed=seed)).to_list()

frequency_dict = {}
dfs = []
for algo in ["CH", "BRKGA", "VNS"]:
    for inst in shuffled_instances:
        frequency = {20: 0, 40: 0, 80: 0, 160: 0}
        for i in range(1, 35):
            with open(f"results/dr_backup/U_{inst}/solution{algo}{inst}_{i}.json") as f:
                sol_json = json.load(f)
                for channel in sol_json["channels"]:
                    frequency[channel[0]] += 1
                    dfs.append(
                        pl.DataFrame(
                            {
                                "algo": algo,
                                "group": inst,
                                "bandwidth": channel[0],
                            }
                        )
                    )

        frequency_dict[algo, inst] = frequency

master_dfs = pl.concat(dfs)
# master_dfs_tmp.group_by("algo", "group", "bandwidth").agg(pl.len() / 34).with_columns(
#     len=pl.when(pl.col.algo == "CH", pl.col.group == 16)
#     .then(pl.lit(0.001))
#     .otherwise("len")
# ).write_csv("seila2.csv")
breakpoint()
