import json
import polars as pl

L = 128
ans = []
for algo in ["BRKGA", "VNS"]:
    inner_ans = []
    for i in range(1, 35):
        with open(f"results/dr_backup/U_{L}/solution{algo}{L}_{i}.json", "r") as _file:
            sol = json.load(_file)

            sum = 0
            for ch in range(len(sol["channels"])):
                sum += len(sol["channels"][ch][1])

            inner_ans.append(sum)

    print(inner_ans)
    ans.append(inner_ans)

df = (
    pl.DataFrame({"vns": ans[1], "brkga": ans[0]})
    .with_columns(dev=(pl.col.brkga - pl.col.vns) / pl.col.vns)
    .with_columns(vns_found_more=pl.col.dev < 0.0)
)
print(df["vns_found_more"].value_counts())
breakpoint()
