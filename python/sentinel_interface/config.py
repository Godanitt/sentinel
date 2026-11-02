"""The small user-facing configuration; paths are relative to the repository."""
from pathlib import Path
import math
import yaml

ROOT = Path(__file__).resolve().parents[2]

def load(path):
    path = Path(path).resolve()
    data = yaml.safe_load(path.read_text())
    schema = {
        "grid": {"size_mm", "electrons_per_cell"},
        "bounds": {"unit", "r", "z"},
        "electric_field": {"file"},
        "output": {"directory", "save_trajectories", "save_endpoints", "save_summary"},
    }
    if not isinstance(data, dict) or set(data) != set(schema):
        raise ValueError("El YAML debe contener grid, bounds, electric_field y output.")
    for key, fields in schema.items():
        if not isinstance(data[key], dict) or set(data[key]) != fields:
            raise ValueError(f"Campos incorrectos en {key}: se esperan {sorted(fields)}")
    def number(value):
        return type(value) in (int, float) and math.isfinite(value)
    size = data["grid"]["size_mm"]
    count = data["grid"]["electrons_per_cell"]
    if not number(size) or size <= 0:
        raise ValueError("grid.size_mm debe ser positivo y finito.")
    if type(count) is not int or count <= 0:
        raise ValueError("electrons_per_cell debe ser un entero positivo.")
    bounds = data["bounds"]
    if bounds["unit"] not in ("mm", "cm"):
        raise ValueError("bounds.unit debe ser mm o cm.")
    factor = 10 if bounds["unit"] == "cm" else 1
    for axis in ("r", "z"):
        values = bounds[axis]
        if (not isinstance(values, list) or len(values) != 2
                or not all(number(v) for v in values) or values[0] >= values[1]):
            raise ValueError(f"bounds.{axis} debe ser [mínimo, máximo], finito y creciente.")
        bounds[axis] = [v * factor for v in values]
    if bounds["r"][0] < 0:
        raise ValueError("El radio inicial no puede ser negativo.")
    for key in ("save_endpoints", "save_trajectories", "save_summary"):
        if type(data["output"][key]) is not bool:
            raise ValueError(f"output.{key} debe ser true o false.")
    for section, key in (("electric_field", "file"), ("output", "directory")):
        value = data[section][key]
        if not isinstance(value, str) or not value.strip():
            raise ValueError(f"{section}.{key} debe ser una ruta.")
        data[section][key] = (ROOT / value).resolve()
    if not data["electric_field"]["file"].is_file():
        raise ValueError(f"No existe el campo: {data['electric_field']['file']}")
    return data
