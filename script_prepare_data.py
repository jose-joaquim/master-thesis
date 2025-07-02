import polars as pl

seed = 1
shuffled_instances = (2 ** pl.int_range(3, 12, eager=True).shuffle(seed=seed)).to_list()

dfs = {}
for inst in shuffled_instances:
    dfs[inst] = pl.read_csv(f"results/dr/U_{inst}/consolidated{inst}.csv").with_columns(
        Group=pl.lit(inst)
    )

concat = pl.concat(dfs.values())
master = (
    concat.drop(["vns_over_ch", "vns_over_brkga"])
    .unpivot(index=["instance", "Group"], on=["ch", "vns", "brkga"])
    .rename({"value": "Result", "instance": "Instance", "variable": "Algorithm"})
)
