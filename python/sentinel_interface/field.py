"""Recognize COMSOL vector exports without changing coordinates or signs."""
import re
import numpy as np

def prepare(source, destination):
    with source.open() as stream:
        header = []
        for line in stream:
            if not line.startswith("%"):
                break
            header.append(line.strip())
    columns = next((line for line in header if re.match(r"%\s+r\s+z\s+", line)), "")
    if "Streamline" in columns:
        # Original SENTINEL examples: mm, kV/cm, inverted z + 1 mm.
        # The original Rust reader and conf/test.toml retain these conventions.
        if source.suffix != ".dat":
            raise ValueError("Para líneas de campo externas, configura explícitamente las unidades en conf/test.toml.")
        return "legacy", source
    names = re.split(r"\s{2,}", columns.lstrip("% "))
    required = ["r", "z", "es.Er (V/m)", "es.Ez (V/m)"]
    if not all(name in names for name in required):
        raise ValueError("Exportación incompatible: se requieren r, z, es.Er (V/m) y es.Ez (V/m). Ephi no sustituye a Ez.")
    unit_line = next((line for line in header if "Length unit:" in line), "")
    unit = unit_line.split(":")[-1].strip()
    if unit not in ("m", "mm", "cm"):
        raise ValueError("La cabecera COMSOL debe indicar Length unit: m, mm o cm.")
    data = np.loadtxt(source, comments="%", ndmin=2)
    if data.shape[1] != len(names):
        raise ValueError("El número de columnas no coincide con la cabecera COMSOL.")
    coordinates = data[:, [names.index("r"), names.index("z")]]
    if not np.isfinite(coordinates).all() or (coordinates[:, 0] < 0).any():
        raise ValueError("Coordenadas no finitas o radios negativos en el mapa.")
    if "es.Ephi (V/m)" in names:
        ephi = data[:, names.index("es.Ephi (V/m)")]
        if np.any(np.isfinite(ephi) & (ephi != 0)):
            raise ValueError("Este lector inicial requiere Ephi = 0; no se descartará una componente azimutal no nula.")
    vectors = data[:, [names.index("es.Er (V/m)"), names.index("es.Ez (V/m)")]]
    valid = np.isfinite(vectors).all(axis=1)
    if not valid.any():
        raise ValueError("El mapa no contiene vectores de campo válidos.")
    # Retain NaNs: the Rust transport stops at the nearest invalid sample.
    vectors[~valid] = np.nan
    coordinates *= {"m": 1000, "cm": 10, "mm": 1}[unit]
    np.savetxt(destination, np.column_stack((coordinates, vectors * 0.01)), fmt="%.17g")
    print(f"Campo COMSOL: {valid.sum()} vectores válidos; {(~valid).sum()} muestras sin datos. z sin invertir.", flush=True)
    return "vector", destination
