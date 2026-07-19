from pathlib import Path

from lambforce_ec.io import save_artery_npz, save_result
from lambforce_ec.synthetic import make_synthetic_artery
from lambforce_ec.workflow import run_case

runtime = Path("runtime_synthetic")
runtime.mkdir(exist_ok=True)
record = make_synthetic_artery()
save_artery_npz(record, runtime / "synthetic_input.npz")
result = run_case(record)
save_result(result, runtime / "result")
print(result.metadata["configuration_checksum"])
