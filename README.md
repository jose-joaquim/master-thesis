# master-thesis

Solver inputs are JSON files (for example,
`instances/du/U_8/MD-VRBSP_U_8_1.json`). Pass the file path as the first
argument to `vns_pure`, `vrbsp_ip`, `mdvrbsp_col_gen`, or `brkga_decoder`.
The MD-VRBSP executable constructs the same `.json` path from its size,
instance index, and instance-directory arguments.

The loader uses nlohmann/json's `from_json` hook to parse `N`, `gamma`,
`noise`, `power_sender`, `receivers`, `senders`, `data_rates`, `SINR`, and
`overlap`. It derives distance, interference, bandwidth requirements and
channel sets. Missing `alfa` defaults to 3 and missing `T` defaults to `N`.
`to_json` serializes the input fields for round-tripping.
