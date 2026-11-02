#!/usr/bin/env python3
"""Lanzador de campañas YAML (terminal, botón Run o IPython).

Edita CAMPAIGNS más abajo para definir las campañas que se ejecutan por defecto.
Sin argumentos: ejecuta todas esas campañas, una tras otra.
Campaña concreta: python scripts/run_simulation.py --config configs/mi_campana.yaml
Solo listar: python scripts/run_simulation.py --list-configs
En IPython: %run scripts/run_simulation.py

El campo, las condiciones y la carpeta de resultados se definen en cada YAML.
Las rutas relativas del YAML parten de la raíz del repositorio.
"""
from pathlib import Path
import sys

# Campañas que se ejecutan por defecto. Añade aquí rutas relativas a la raíz.
CAMPAIGNS = [
    "configs/testing_100um.yaml",
    "configs/testing_250um.yaml",
]

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))

if __name__ == "__main__":
    import importlib
    from sentinel_interface import runner
    importlib.reload(runner)
    runner.main(campaigns=CAMPAIGNS)
