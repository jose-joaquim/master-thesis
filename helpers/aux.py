import polars as pl

for i in [8, 16, 32, 64, 128, 256, 512, 1024, 2048]:
    df = pl.read_csv(f"results/dr_backup/U_{i}/consolidated{i}.csv")

    df = df.with_columns(
        ((pl.col.ch - pl.col.brkga) / pl.col.brkga).alias("ch_over_ga"),
        ((pl.col.vns - pl.col.brkga) / pl.col.brkga).alias("vns_over_ga"),
    )

    print(
        f'{i}\t\t\t{df["ch_over_ga"].mean() * 100}\t\t\t{df["vns_over_ga"].mean() * 100}'
    )
