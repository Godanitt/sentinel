[![built with nix](https://builtwithnix.org/badge.svg)](https://builtwithnix.org)

# sentinel

**S**implified **E**lectro**N** **T**ransport **I**n **N**obl**E** **L**iquids is a Rust toolkit for simulating electron transport and detector readout in a noble-liquid detector. It models drift through a field map in a tapered cylindrical geometry, converts ionization energy deposits into SiPM hits using light tables, and digitizes those hits into DAQ-like waveforms.

The repository contains four command-line programs:

| Program | Purpose | Input | Output |
| --- | --- | --- | --- |
| `trajectories` | Propagate sampled electrons and write every drift point. | Field/configuration files | HDF5 `/trajectories` or whitespace-separated text |
| `endpoints` | Propagate sampled electrons and write only start/end points. | Field/configuration files | HDF5 `/endpoints` or whitespace-separated text |
| `detsim` | Turn ionization deposits into simulated S1/S2 SiPM hits. | MC HDF5 + light-table HDF5 | HDF5 `/MC/sensor_hits` |
| `waveforms` | Turn SiPM hits into 125-channel digitized waveforms. | HDF5 `/MC/sensor_hits` | HDF5 `/Waveforms/*` |

The normal detector-readout pipeline is:

```text
MC /MC/ionization_hits + light tables /S1/LT and /S2/LT
                              │
                              ▼
                         detsim
                              │  /MC/sensor_hits
                              ▼
                         waveforms
                              │  /Waveforms/adcs
                              ▼
                         DAQ-like data
```

`trajectories` and `endpoints` are independent transport tools and are useful for checking field maps, geometry, diffusion, and drift times.

## Installation

The supported development environment is Nix on `x86_64-linux`. The flake supplies the Rust toolchain pinned by [rust-toolchain.toml](../rust-toolchain.toml), HDF5 development libraries, `just`, `cargo-nextest`, and `bacon`.

### Recommended: Nix and direnv

Install Nix with flakes enabled. The [Determinate Nix Installer](https://determinate.systems/nix-installer) provides instructions for Linux, macOS, and WSL. Install `direnv` if it is not already available:

```sh
nix profile add nixpkgs#direnv
```

Enable its shell hook once. For zsh:

```sh
eval "$(direnv hook zsh)"
```

For bash:

```sh
eval "$(direnv hook bash)"
```

You probably want to enable the hook automatically whenever you start a new shell. Add the appropriate command to your shell configuration file:

```sh
# zsh
echo 'eval "$(direnv hook zsh)"' >> ~/.zshrc
# or: echo 'eval "$(direnv hook zsh)"' >> ~/.zsh_profile

# bash
echo 'eval "$(direnv hook bash)"' >> ~/.bashrc
# or: echo 'eval "$(direnv hook bash)"' >> ~/.bash_profile
```

From the repository root, allow the checked-in `.envrc` and enter the development environment:

```sh
direnv allow
```

The first activation may download flake inputs and build the environment, which might take a while. The `.envrc` runs `use flake`, so subsequent shells in this directory automatically receive the configured Rust and HDF5 environment.

If `direnv` is not available, activate the same environment manually:


```sh
nix develop
```

## Build

Build all four binaries in debug mode:

```sh
just --justfile dev/justfile build
```
The executables are written to `target/debug/` or `target/release/`. Use the `just --justfile dev/justfile release` wrapper to (compile and) run a release binary and see its command-line options:

```sh
just --justfile dev/justfile release trajectories --help
just --justfile dev/justfile release endpoints --help
just --justfile dev/justfile release detsim --help
just --justfile dev/justfile release waveforms --help
```

`just` provides short wrappers around the common Cargo commands:

```sh
just --justfile dev/justfile build
just --justfile dev/justfile debug trajectories --conf conf/test.toml --n-events 2 --batch-size 2 --output /tmp/trajectories.h5 --batch
just --justfile dev/justfile release endpoints --conf conf/test.toml --n-events 2 --batch-size 2 --output /tmp/endpoints.h5 --batch
just --justfile dev/justfile clean
```

All commands below assume they are run from the repository root. The programs resolve relative paths such as `data/...` from the current working directory, not from the directory containing the configuration file.

## Quick start: self-contained transport example

[conf/test.toml](../conf/test.toml) uses the bundled homogeneous field map, a small cone, and a fixed starting point. It is intended for a smoke test rather than a production detector configuration. The last batch is truncated to the requested event count:

```sh
mkdir -p outputs

just --justfile dev/justfile release endpoints \
  --conf conf/test.toml \
  --n-events 2 \
  --batch-size 2 \
  --output outputs/endpoints.h5 \
  --threads 2 \
  --batch
```

To write a human-readable text file instead:

```sh
just --justfile dev/justfile release trajectories \
  --conf conf/test.toml \
  --n-events 2 \
  --batch-size 2 \
  --output outputs/trajectories.txt \
  --threads 2 \
  --csv \
  --batch
```

The `--csv` output is whitespace-separated text, not comma-separated CSV. Its headers are `event x y z t` for trajectories and `event x0 y0 z0 x1 y1 z1 t` for endpoints.

HDF5 command-line tools are available in the Nix environment. For example:

```sh
h5dump -n outputs/endpoints.h5
```

Outputs are not overwritten by default. Add `--overwrite` when intentionally replacing an existing file; parent directories must already exist.

## Transport programs

### Geometry and field propagation

The transport code uses a tapered cylindrical (cone-like) active region:

```text
z = 0 .. neck_length + drift_length
r(z) = rmin + form_factor * max(z - neck_length, 0)
```

At every step, the tracker finds the nearest point in the supplied axisymmetric field map, interpolates the liquid-xenon drift velocity from [data/drift_velocity_gushchin.dat](../data/drift_velocity_gushchin.dat), and adds transverse and longitudinal diffusion. Positions are in millimetres and drift times are in microseconds after the configured field conversion.

`trajectories` writes one record for every propagated point. `endpoints` writes one record per sampled electron. Both use the same configuration and support:

| Option | Meaning | Default |
| --- | --- | --- |
| `-c, --conf PATH` | TOML configuration. | required |
| `-n, --n-events N` | Number of requested events. | required |
| `-b, --batch-size N` | Number of events processed per write batch. | `1000` |
| `-o, --output PATH` | Output HDF5 or text path. | required |
| `-t, --threads N` | Rayon worker threads. | `10` |
| `--csv` | Write whitespace-separated text instead of HDF5. | off |
| `--batch` | Suppress the interactive progress bar. | off |
| `--overwrite` | Allow an existing output file to be replaced. | off |

The HDF5 datasets have these fields:

```text
/trajectories: event, x, y, z, t
/endpoints:    event, x0, y0, z0, x1, y1, z1, t
```

The current batching loop processes a complete `--batch-size` block even for the final partial block. For an exact event count, choose a batch size that divides `--n-events`.

### Transport configuration

The transport configuration is loaded by both `trajectories` and `endpoints`:

| Key | Meaning |
| --- | --- |
| `rmin` | Radius of the neck, in mm. |
| `form_factor` | Radial increase per mm above the neck. |
| `neck_length` | Cylindrical neck length, in mm. |
| `drift_length` | Conical drift length, in mm. |
| `field_file` | COMSOL-style field export. The loader skips eight header lines and expects `r z streamline field` columns. |
| `field_to_mm` | Coordinate conversion applied to the field file. |
| `field_to_Vpercm` | Field-magnitude conversion applied to the field file. |
| `field_invert_z` | Stored in the configuration; the current binaries always enable z inversion when loading the field, so changing this value currently has no runtime effect. |
| `t_step` | Propagation step, in microseconds. |
| `generator` | Initial-position sampler; see below. |

Available generator values are:

```toml
generator = "Volume"
generator = "Surface"
generator = "CathodeCenter"
generator = "Center"
generator = { FixedPos = [0.0, 0.0, 9.0] }
generator = { Plane = 1.57079632679 } # angle in radians
generator = { Edge = 0.1 }             # distance from wall/cathode, in mm
generator = { FromFile = "input.h5" } # reads /MC/ionization_hits positions
```

`Volume`, `Surface`, `Plane`, and `Edge` sample the configured cone. `FromFile` consumes ionization hits from the HDF5 dataset `/MC/ionization_hits` in sequence.

## Detector simulation

`detsim` converts ionization deposits into S1 and S2 sensor hits. It requires:

1. An input HDF5 file containing `/MC/ionization_hits`.
2. A light-table HDF5 file containing `/S1/LT` and `/S2/LT`.
3. A detector configuration containing all transport fields plus `lifetime` in microseconds.

The repository includes [data/kr_hits.h5](../data/kr_hits.h5), a small example MC file. It does not include a light-table file. The required light-table records are:

```text
/S1/LT: sensor_id, z, x, y, nevt, prob
/S2/LT: sensor_id, x, y, nevt, prob
```

The input ionization records are:

```text
/MC/ionization_hits: event, track_id, x, y, z, e
```

Here `x`, `y`, and `z` are positions in mm and `e` is an energy deposit in MeV. Events must be grouped contiguously in the input file because the program groups consecutive records by event ID.

The production-like [conf/pcolina.toml](../conf/pcolina.toml) contains a machine-specific absolute `field_file` path. Make a local copy and update that path before running, for example:

```sh
cp conf/pcolina.toml conf/pcolina-local.toml
# edit field_file in conf/pcolina-local.toml
```

Then run:

```sh
mkdir -p outputs
just --justfile dev/justfile release detsim \
  --conf conf/pcolina-local.toml \
  --input data/kr_hits.h5 \
  --maps /absolute/path/to/light-tables.h5 \
  --output outputs/sensor_hits.h5 \
  --medium xenon \
  --threads 8 \
  --batch
```

The output dataset is:

```text
/MC/sensor_hits: event, sensor_id, t
```

`t` is in microseconds. The default and currently complete detector-simulation medium is Xenon. The CLI also lists Argon and Krypton, but their ionization/scintillation constants are not implemented in `detsim` and will panic if selected.

## Waveform digitization

`waveforms` consumes `/MC/sensor_hits` and applies the following chain independently for each event:

1. Bin sensor hits into time samples, including the configured pre-trigger.
2. Add Poisson dark noise.
3. Convolve with the bundled Xenon or Argon pulse shape.
4. Apply photoelectron gain fluctuations and ADC-per-photoelectron scaling.
5. Add white noise.
6. Quantize to `u16` and apply the configured pedestal convention.

Run it on the detector-simulation output with:

```sh
just --justfile dev/justfile release waveforms \
  --conf conf/electronics.toml \
  --input outputs/sensor_hits.h5 \
  --output outputs/waveforms.h5 \
  --medium xenon \
  --threads 8
```

The electronics configuration contains:

| Key | Meaning |
| --- | --- |
| `buffer_size_mus` | Waveform duration, in microseconds. |
| `sampling_time_mus` | Sampling interval, in microseconds. |
| `sipm_size_mm2` | SiPM active area used for dark-noise probability. |
| `noise_std` | White-noise standard deviation. |
| `adc_per_pe` | ADC scale per photoelectron. |
| `pe_sigma` | Relative gain spread. |
| `pretrigger_mus` | Time before the trigger, in microseconds. |
| `baseline` | `u16` pedestal value. |

The buffer duration must be a multiple of the sampling interval. The supplied configuration produces `6250` samples per channel (`100 µs / 0.016 µs`) for `125` fixed pCOLINA channels. The output datasets are:

```text
/Waveforms/event_id   # one u64 event ID per waveform
/Waveforms/channel_id # 125 u16 sensor IDs, in waveform channel order
/Waveforms/adcs       # u16 array with shape (events, 125, time samples)
```

Waveform input events must also be grouped contiguously. Events with no sensor hits are absent from the input and therefore absent from the waveform output. Xenon and Argon pulse-shape data are bundled; Krypton waveform processing is not implemented.

## Configuration files and bundled data

| Path | Role |
| --- | --- |
| [conf/test.toml](../conf/test.toml) | Small self-contained homogeneous-field transport test. |
| [conf/pcolina.toml](../conf/pcolina.toml) | Reduced pCOLINA geometry; includes detector lifetime but an absolute, machine-specific field path. |
| [conf/colina.toml](../conf/colina.toml) | Full-size geometry template; its field path is machine-specific and it does not include `lifetime`, so add that key before using it with `detsim`. |
| [conf/electronics.toml](../conf/electronics.toml) | Waveform buffer, sampling, noise, gain, and pedestal settings. |
| [data/homogeneous_field.dat](../data/homogeneous_field.dat) | Small COMSOL-style homogeneous field used by the test configuration. |
| [data/partial_efield.dat](../data/partial_efield.dat) | Larger partial field-map export. |
| [data/drift_velocity_gushchin.dat](../data/drift_velocity_gushchin.dat) | Drift-velocity lookup table used by the tracker. |
| `data/pulse_shape_vuv4_lxe.csv` | Xenon pulse shape. |
| `data/pulse_shape_vuv4_lar.csv` | Argon pulse shape. |
| [data/kr_hits.h5](../data/kr_hits.h5) | Example ionization-hit input. |

The field-map loader expects the numeric COMSOL export format used by the bundled files: eight header lines followed by four whitespace-separated columns. If a field map uses different units, adjust `field_to_mm` and `field_to_Vpercm` in a local configuration.

## Development and tests

Useful commands are:

```sh
just --justfile dev/justfile build
just --justfile dev/justfile test
just --justfile dev/justfile test-slow
just --justfile dev/justfile verbose '<nextest filter expression>'
```

`just --justfile dev/justfile test` uses `cargo nextest` and `just --justfile dev/justfile test-slow` uses a release build. Tests that exercise HDF5 require the development environment so that `HDF5_DIR` and the HDF5 libraries are available.

The waveform configuration tests match the current electronics configuration.
The simulation uses random sampling and has no command-line seed option, so repeated runs are not bit-for-bit reproducible.

## License

sentinel is distributed under the [GNU General Public License, version 3](../LICENSE).
