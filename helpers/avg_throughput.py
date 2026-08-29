import polars as pl

for i in [8, 16, 32, 64, 128, 256, 512]:
    df = pl.read_csv(f"gmm{i}.out", has_header=False)
    print(df.mean().item())
