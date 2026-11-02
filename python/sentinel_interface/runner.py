"""Prepare one run, invoke Rust once, and retain portable HDF5 outputs."""
import argparse
from datetime import datetime
from html import escape
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

import h5py
import numpy as np

from .config import ROOT, load
from .field import prepare
from .grid import cells, initial_batches

MAP_DTYPE = np.dtype([("event_id", "<u4"), ("cell_id", "<u8")])

def _notebook_api():
    try:
        from IPython import get_ipython
        from IPython.display import HTML, display
        shell = get_ipython()
        if shell is not None and (hasattr(shell, "kernel") or "zmq" in shell.__class__.__module__.lower()):
            return HTML, display
    except ImportError:
        pass
    return None


class NotebookProgress:
    """Update one self-contained HTML progress card in IPython output."""
    def __init__(self, HTML, display):
        self.HTML = HTML
        self.display = display
        self.handle = None
        self.tasks = {}

    @staticmethod
    def duration(seconds):
        seconds = max(0, int(seconds))
        return f"{seconds // 3600:02d}:{(seconds // 60) % 60:02d}:{seconds % 60:02d}"

    def update(self, phase, done, total, details):
        task = self.tasks.setdefault(phase, {"started": time.monotonic(), "done": 0, "total": total, "details": ""})
        task.update(done=done, total=total, details=details)
        cards = []
        for key in ("validation", "transport"):
            if key not in self.tasks:
                continue
            item = self.tasks[key]
            elapsed = time.monotonic() - item["started"]
            fraction = min(1.0, item["done"] / item["total"]) if item["total"] else 1.0
            remaining = elapsed * (1.0 - fraction) / fraction if fraction else 0
            title = "Validando el campo" if key == "validation" else "Transportando electrones"
            cards.append(
                '<section style="margin:14px 0 0">'
                '<div style="display:flex;justify-content:space-between;gap:16px;align-items:center;margin-bottom:7px">'
                f'<span style="font-weight:650;color:#e5e7eb">{title}</span>'
                f'<span style="font-variant-numeric:tabular-nums;color:#86efac">{fraction * 100:.1f}%</span></div>'
                '<div style="height:12px;border-radius:99px;overflow:hidden;background:#7f1d1d">'
                f'<div style="height:100%;width:{fraction * 100:.3f}%;border-radius:99px;'
                'background:linear-gradient(90deg,#22c55e,#86efac);transition:width .35s ease"></div></div>'
                '<div style="display:flex;justify-content:space-between;gap:10px;flex-wrap:wrap;'
                'margin-top:7px;color:#a1a1aa;font-size:12px;font-variant-numeric:tabular-nums">'
                f'<span>Transcurrido {self.duration(elapsed)} · Restante ~{self.duration(remaining)}</span>'
                f'<span>{escape(item["details"])}</span></div></section>'
            )
        body = (
            '<div style="font-family:ui-sans-serif,system-ui,-apple-system,Segoe UI,sans-serif;'
            'max-width:900px;padding:16px 18px;border:1px solid #27272a;border-radius:14px;'
            'background:linear-gradient(145deg,#18181b,#111827);color:#f4f4f5;'
            'box-shadow:0 8px 24px #0003">'
            '<div style="display:flex;align-items:center;gap:9px">'
            '<span style="width:9px;height:9px;border-radius:50%;background:#4ade80;'
            'box-shadow:0 0 12px #4ade8080"></span>'
            '<span style="font-size:13px;letter-spacing:.08em;font-weight:750;color:#d4d4d8">SENTINEL · SIMULACIÓN</span>'
            '</div>' + ''.join(cards) + '</div>'
        )
        rendered = self.HTML(body)
        if self.handle is None:
            self.handle = self.display(rendered, display_id=True)
        else:
            self.handle.update(rendered)


def run_transport(command):
    """Stream Rust progress into Rich, with a plain-text fallback."""
    process = subprocess.Popen(command, cwd=ROOT, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT, text=True, bufsize=1)
    notebook = _notebook_api()
    if notebook is not None:
        ui = NotebookProgress(*notebook)
        for line in process.stdout:
            if line.startswith("SENTINEL_PROGRESS\t"):
                _, phase, done_text, total_text, details = line.rstrip("\n").split("\t", 4)
                if phase == "validation":
                    details = "Celdas descartadas: " + details.rsplit("celdas descartadas ", 1)[-1]
                else:
                    parts = details.split(" | ")
                    cells = parts[0].removeprefix("celdas ")
                    discarded = parts[2].removeprefix("descartadas ") if len(parts) > 2 else "0"
                    electrons = parts[3].removeprefix("electrones ") if len(parts) > 3 else ""
                    details = f"Celdas {cells} · descartadas {discarded} · electrones {electrons}"
                ui.update(phase, int(done_text), int(total_text), details)
            elif line.strip():
                print(line, end="", flush=True)
        return_code = process.wait()
        if return_code:
            raise subprocess.CalledProcessError(return_code, command)
        return

    try:
        from rich.console import Console
        from rich.progress import (BarColumn, Progress, SpinnerColumn,
                                   TaskProgressColumn, TextColumn, TimeElapsedColumn,
                                   TimeRemainingColumn)
        from rich.theme import Theme
    except ImportError:
        print(f"Rich no está instalado en {sys.executable}; mostrando progreso simple. Instala con: {sys.executable} -m pip install -r python/requirements.txt", flush=True)
        _run_transport_fallback(process, command)
        return

    if Progress is not None:
        console = Console(theme=Theme({
            "progress.elapsed": "bright_green",
            "progress.remaining": "yellow",
            "progress.spinner": "bright_green",
        }))
        columns = [
            SpinnerColumn("dots", style="bold bright_green"),
            TextColumn("[bold]{task.description}[/bold]"),
            BarColumn(bar_width=30, style="bold red", complete_style="bold bright_green",
                      finished_style="bold green", pulse_style="yellow"),
            TaskProgressColumn(text_format="[bold bright_green]{task.percentage:>3.0f}%[/bold bright_green]"),
            TimeElapsedColumn(),
            TextColumn("[dim]{task.fields[info]}[/dim]"),
            TimeRemainingColumn(),
        ]
        with Progress(*columns, console=console, refresh_per_second=8, transient=False) as progress:
            task = None
            phase = None
            for line in process.stdout:
                if line.startswith("SENTINEL_PROGRESS\t"):
                    _, next_phase, done_text, total_text, details = line.rstrip("\n").split("\t", 4)
                    done, total = int(done_text), int(total_text)
                    if next_phase != phase:
                        description = ("[bright_green]Validando campo[/bright_green]" if next_phase == "validation"
                                       else "[bright_green]Transportando electrones[/bright_green]")
                        task = progress.add_task(description, total=total, info="")
                        phase = next_phase
                    if next_phase == "validation":
                        info = "[red]omitidas " + details.rsplit("celdas descartadas ", 1)[-1] + " celdas[/red]"
                    else:
                        parts = details.split(" | ")
                        cell_progress = parts[0].removeprefix("celdas ")
                        discarded = parts[2].removeprefix("descartadas ") if len(parts) > 2 else "0"
                        info = f"celdas {cell_progress} · [red]omitidas {discarded}[/red]"
                    progress.update(task, completed=done, info=info)
                elif line.strip():
                    console.print(line.rstrip())
    else:
        _run_transport_fallback(process, command)
        return
    return_code = process.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, command)


def _run_transport_fallback(process, command):
    for line in process.stdout:
        if line.startswith("SENTINEL_PROGRESS\t"):
            _, next_phase, done_text, total_text, details = line.rstrip("\n").split("\t", 4)
            description = "Validando campo" if next_phase == "validation" else "Transportando"
            print(f"{description}: {done_text}/{total_text} | {details}", flush=True)
        elif line.strip():
            print(line, end="", flush=True)
    return_code = process.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, command)

def executable():
    binary = ROOT / "target/release/simulate"
    sources = list((ROOT / "src").rglob("*.rs")) + [ROOT / "Cargo.toml", ROOT / "Cargo.lock"]
    if binary.exists() and binary.stat().st_mtime >= max(p.stat().st_mtime for p in sources):
        return binary
    cargo = shutil.which("cargo")
    if cargo:
        command = [cargo, "build", "--release", "--bin", "simulate"]
    elif shutil.which("nix"):
        command = ["nix", "develop", "--command", "cargo", "build", "--release", "--bin", "simulate"]
    else:
        raise RuntimeError("Compila una vez desde WSL: nix develop --command cargo build --release --bins")
    subprocess.run(command, cwd=ROOT, check=True)
    return binary

def write_initial(path, table, count):
    total = len(table) * count
    with h5py.File(path, "w") as file:
        file.attrs["position_unit"] = "mm"
        file.attrs["event_meaning"] = "unique electron identifier"
        file.create_dataset("cells", data=table)
        dataset = file.create_dataset("initial_positions", (total,), dtype=next(initial_batches(table[:1], 1)).dtype,
                                     chunks=True, compression="gzip")
        mapping = file.create_dataset("event_cells", (total,), dtype=MAP_DTYPE, chunks=True, compression="gzip")
        start = 0
        for batch in initial_batches(table, count):
            stop = start + len(batch)
            dataset[start:stop] = batch
            ids = np.empty(len(batch), dtype=MAP_DTYPE)
            ids["event_id"], ids["cell_id"] = batch["event"], batch["cell_id"]
            mapping[start:stop] = ids
            start = stop

def collect(raw, initial, output, config, table):
    counts = np.zeros(len(table), dtype=np.uint64)
    with h5py.File(raw) as source, h5py.File(initial) as inputs:
        processed = source["processed"]
        accepted_events = processed["event_id"][:]
        candidate_mapping = inputs["event_cells"]
        accepted_mapping = (candidate_mapping[accepted_events] if len(accepted_events)
                            else np.empty(0, dtype=candidate_mapping.dtype))
        for start in range(0, len(processed), 8192):
            block = processed[start:start+8192]
            event_ids = block["event_id"].astype(np.int64)
            expected = inputs["event_cells"][event_ids]
            if not np.array_equal(block, expected):
                raise RuntimeError("La correspondencia event_id/cell_id cambió durante la simulación.")
            counts += np.bincount(block["cell_id"].astype(np.int64), minlength=len(table)).astype(np.uint64)
        for name in ("endpoints", "trajectories"):
            if config["output"][f"save_{name}"]:
                with h5py.File(output / f"{name}.h5", "w") as destination:
                    source.copy(name, destination)
                    destination.create_dataset("event_cells", data=accepted_mapping)
                    destination.attrs["position_unit"] = "mm"
                    destination.attrs["time_unit"] = "us"
        if config["output"]["save_summary"]:
            dtype = np.dtype([
                ("cell_id", "<u8"), ("r_center_mm", "<f8"), ("z_center_mm", "<f8"),
                ("n_electrons_requested", "<u8"), ("n_electrons_generated", "<u8"),
                ("n_electrons_processed", "<u8"),
            ])
            summary = np.empty(len(table), dtype=dtype)
            for name in ("cell_id", "r_center_mm", "z_center_mm"):
                summary[name] = table[name]
            summary["n_electrons_requested"] = config["grid"]["electrons_per_cell"]
            summary["n_electrons_generated"] = counts
            summary["n_electrons_processed"] = counts
            with h5py.File(output / "summary.h5", "w") as destination:
                destination.create_dataset("summary", data=summary)
                destination.create_dataset("event_cells", data=accepted_mapping)

        # Keep only electrons from cells accepted by Rust in the portable
        # initial file. Their event IDs stay stable for joins with endpoints.
        filtered = initial.with_name("initial_positions.filtered.h5")
        with h5py.File(filtered, "w") as destination:
            for key, value in inputs.attrs.items():
                destination.attrs[key] = value
            inputs.copy("cells", destination)
            positions = inputs["initial_positions"]
            mapping = inputs["event_cells"]
            out_positions = destination.create_dataset(
                "initial_positions", shape=(len(accepted_events),), dtype=positions.dtype,
                chunks=True, compression="gzip")
            out_mapping = destination.create_dataset(
                "event_cells", shape=(len(accepted_events),), dtype=mapping.dtype,
                chunks=True, compression="gzip")
            for start in range(0, len(accepted_events), 8192):
                ids = accepted_events[start:start + 8192].astype(np.int64)
                stop = start + len(ids)
                out_positions[start:stop] = positions[ids]
                out_mapping[start:stop] = mapping[ids]
    filtered.replace(initial)
    return int(counts.sum())

def run(config_path):
    config_path = Path(config_path).resolve()
    config = load(config_path)
    print(f"Campaign YAML: {config_path}", flush=True)
    print(f"Electric field: {config['electric_field']['file']}", flush=True)
    campaign_name = config_path.stem
    base = ROOT / "outputs" / campaign_name / datetime.now().strftime("%Y_%m_%d")
    base.mkdir(parents=True, exist_ok=True)
    print(f"Output base directory: {base}", flush=True)
    table = cells(config)
    binary = executable()
    run_number = 1
    while True:
        output = base / f"run{run_number:02d}"
        try:
            output.mkdir()
            break
        except FileExistsError:
            run_number += 1
    (output / "configs").mkdir()
    shutil.copy2(config_path, output / "configs/initial_conditions.yaml")
    initial = output / "initial_positions.h5"
    total = len(table) * config["grid"]["electrons_per_cell"]
    print(f"Number of cells: {len(table)}\nElectrons per cell: {config['grid']['electrons_per_cell']}\nTotal electrons: {total}", flush=True)
    started = time.perf_counter()
    try:
        with tempfile.TemporaryDirectory(prefix="sentinel_") as work:
            work = Path(work)
            mode, field = prepare(config["electric_field"]["file"], work / "field.dat")
            write_initial(initial, table, config["grid"]["electrons_per_cell"])
            args = [str(binary), "--input", str(initial), "--field", str(field),
                    "--output", str(work / "transport.h5"), "--mode", mode]
            for name in ("endpoints", "trajectories"):
                if config["output"][f"save_{name}"]:
                    args.append(f"--{name}")
            run_transport(args)
            processed_count = collect(work / "transport.h5", initial, output, config, table)
        elapsed = time.perf_counter() - started
        field_path = config["electric_field"]["file"]
        with field_path.open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        (output / "run.txt").write_text(
            f"campaign_config: {config_path}\nfield: {field_path}\nfield_sha256: {digest}\nmode: {mode}\n"
            f"electrons_processed: {processed_count}\nsimulation_time_s: {elapsed:.6f}\n")
    except Exception as error:
        (output / "FAILED.txt").write_text(str(error) + "\n")
        raise RuntimeError(f"{error}\nEjecución incompleta conservada en: {output}") from error
    print(f"Simulation time: {elapsed:.3f} s\nElectrons processed: {processed_count}\nElectrons per second: {processed_count/elapsed:.1f}\nOutput directory: {output}")
    return output

def available_configs():
    """Every YAML under configs/ is an independently selectable campaign."""
    return sorted(path for path in (ROOT / "configs").rglob("*")
                  if path.is_file() and path.suffix.lower() in (".yaml", ".yml"))


def main(argv=None, campaigns=None):
    parser = argparse.ArgumentParser(description="Transporte 3D de electrones mediante SENTINEL")
    parser.add_argument("--config", type=Path, nargs="+", help="Uno o más YAML de campaña; por defecto se ejecutan las campañas configuradas en el lanzador.")
    parser.add_argument("--list-configs", action="store_true", help="Listar las campañas de configs/ sin simular.")
    if argv is None:
        # Jupyter añade --f=<kernel.json> al argv del proceso. No es un
        # argumento de la aplicación, así que lo retiramos antes de parsear.
        argv = sys.argv[1:]
        argv = [arg for arg in argv if not arg.startswith("--f=")]
        if "-f" in argv:
            index = argv.index("-f")
            del argv[index:index + 2]
    args = parser.parse_args(argv)
    try:
        configured = ([Path(path) for path in campaigns] if campaigns is not None
                      else available_configs())
        configured = [path if path.is_absolute() else ROOT / path for path in configured]
        if args.list_configs:
            for path in configured:
                print(path.relative_to(ROOT))
            return
        configs = args.config if args.config is not None else configured
        for config in configs:
            run(config)
    except KeyboardInterrupt:
        print("\nEjecución cancelada.", file=sys.stderr)
    except (ValueError, OSError, RuntimeError, subprocess.CalledProcessError) as error:
        sys.exit(f"Error: {error}")
