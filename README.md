# SENTINEL

Transporte de electrones en líquidos nobles con un motor Rust y una interfaz Python sencilla. Python prepara las posiciones iniciales; Rust mantiene la deriva, la difusión y la propagación en tres dimensiones.

## Empezar en Ubuntu / WSL

Desde la raíz del repositorio:

```bash
nix develop
cargo build --release --bins
/home/colina/.venv/bin/python -m pip install -r python/requirements.txt
/home/colina/.venv/bin/python scripts/run_simulation.py --config configs/initial_conditions.yaml
```

Se utiliza el entorno Python existente, sin crear otro. Si Cargo está disponible, el lanzador compila una vez cuando falta el binario release o cambia el código. Fuera del entorno Nix, el ejecutable debe poder encontrar sus bibliotecas HDF5; ejecuta el comando dentro de `nix develop`.

## Elegir una campaña

Cada YAML de `configs/` representa una campaña independiente. Duplica `initial_conditions.yaml` con el nombre que prefieras, por ejemplo `tpc_campo_a.yaml` y `tpc_campo_b.yaml`. En cada copia cambia `electric_field.file`, las condiciones iniciales y `output.directory` (por ejemplo `outputs/tpc_campo_a`). No hace falta editar Python.

Desde WSL, dentro de `nix develop`, ejecutar sin argumentos abre un selector:

~~~bash
/home/colina/.venv/bin/python scripts/run_simulation.py
~~~

También puedes elegir directamente o listar sin ejecutar:

~~~bash
/home/colina/.venv/bin/python scripts/run_simulation.py --config configs/initial_conditions.yaml
/home/colina/.venv/bin/python scripts/run_simulation.py --list-configs
~~~

Desde IPython, situado en la raíz del proyecto y usando el mismo entorno WSL/Nix. El lanzador descarta el `--f=...` que Jupyter añade al kernel:

~~~ipython
%run scripts/run_simulation.py
%run scripts/run_simulation.py --config configs/initial_conditions.yaml
~~~

El selector permite introducir varios números separados por comas o `*` para ejecutar todas las campañas. También se pueden indicar varios YAML en `--config`:

~~~bash
/home/colina/.venv/bin/python scripts/run_simulation.py --config configs/tpc_campo_a.yaml configs/tpc_campo_b.yaml
~~~

El botón **Run Python File** de VS Code también muestra el selector si utiliza el intérprete `/home/colina/.venv/bin/python` y la terminal del entorno Nix. Antes de simular se muestran el YAML elegido, el campo y la carpeta base de resultados.

## Organización

```text
configs/                    Solo YAML: condiciones iniciales y prueba pequeña
conf/                       TOML internos de los programas originales
data/fields/Full_output.txt  Exportación original de TPC_Avanzada
data/                       Tablas físicas y datos de ejemplo existentes
python/sentinel_interface/  Configuración, rejilla y ejecución
python/tests/               Pruebas de integración
scripts/run_simulation.py   Punto de entrada
src/                        Motor Rust y ejecutables
outputs/                    Resultados; cada ejecución en su propia carpeta
outputs/legacy/             Resultados anteriores conservados
dev/                        Recetas just y entrada Nix heredada
docs/original-cli.md         Referencia de los cuatro programas originales
```

Los manifiestos Cargo, el toolchain Rust y la flake Nix permanecen en la raíz para que las herramientas los encuentren automáticamente.

## Configuración y salidas

Edita [configs/initial_conditions.yaml](configs/initial_conditions.yaml). Las rutas relativas del YAML se resuelven desde la raíz del repositorio, incluso al lanzar Python desde otra carpeta. El tamaño de celda siempre se expresa en mm; `bounds.unit` admite mm o cm. Las celdas del borde se recortan si el intervalo no es múltiplo del tamaño.

Cada celda recibe exactamente `electrons_per_cell` posiciones aleatorias uniformes en volumen anular, con `r >= 0`, `phi = 0`, `x = r`, `y = 0`. La difusión posterior permite salir de ese plano. El identificador `event` identifica un electrón; `cell_id` identifica su celda. No dependen del orden de finalización de los hilos.

El ejemplo inicial usa una celda pequeña dentro del mapa, con 1000 electrones. Amplía los límites cuando hayas seleccionado la región de interés. Para una primera comprobación de solo tres electrones:

```bash
/home/colina/.venv/bin/python scripts/run_simulation.py --config configs/smoke_test.yaml
```

Los resultados quedan bajo `outputs/<campaña>/AAAA_MM_DD/runNN/`, usando el nombre del YAML como campaña (por ejemplo `outputs/initial_conditions/2026_05_12/run01/`). El contador continúa con `run02`, `run03`, etc. si ya existe una carpeta del día para esa campaña. Nunca se reutiliza una carpeta de ejecución:

- `initial_positions.h5`: datasets `initial_positions`, `cells` y `event_cells`, siempre conservados.
- `endpoints.h5`: dataset original `endpoints`, solo si se activa.
- `trajectories.h5`: dataset original `trajectories`, solo si se activa.
- `summary.h5`: recuentos administrativos por celda, solo si se activa.
- `configs/initial_conditions.yaml`: copia exacta de la configuración utilizada.
- `run.txt`: ruta y SHA256 del campo, modo y tiempo de ejecución.

Todas las salidas activadas incluyen `event_cells` con `event_id` y `cell_id`. Si alguna posición inicial generada en una celda cae fuera del soporte del mapa o sobre una muestra de campo inválida/nula, se descarta la celda completa y el resto continúa. `initial_positions.h5` conserva solo los electrones de celdas aceptadas; `summary.h5` refleja cuántos se generaron y procesaron por celda. Los nombres y tipos originales de endpoints y trayectorias se conservan. Las posiciones son mm y los tiempos µs. El tiempo del endpoint coincide con el último punto de su trayectoria. Si se activan ambas salidas, se simula una sola vez. Solo endpoints no almacena los pasos intermedios. Si falla el transporte de un electrón después de iniciar, se interrumpe la ejecución y se conserva `FAILED.txt`.

Lectura desde Python, con Pandas si está instalado:

```python
from pathlib import Path
import h5py
import pandas as pd

run = Path("outputs/initial_conditions/2026_05_12/run01")
with h5py.File(run / "endpoints.h5") as f:
    endpoints = pd.DataFrame.from_records(f["endpoints"][:])
    mapping = pd.DataFrame.from_records(f["event_cells"][:])
endpoints = endpoints.merge(mapping, left_on="event", right_on="event_id", validate="one_to_one")
```

## Campo y límites del modelo

El archivo COMSOL adjunto se conserva sin cambios. Su cabecera contiene `r, z, es.Er, es.Ephi, es.Ez, V, es.normE`: el lector convierte m a mm y V/m a V/cm, mantiene el eje z sin invertir y deriva en la dirección opuesta a `(Er, Ez)`. `Ephi` es cero en este mapa; una componente azimutal no nula se rechaza explícitamente. El módulo vectorial calcula la magnitud a partir de Er/Ez y conserva las transformaciones cilíndricas del motor. No incorpora modelos nuevos de transporte.

Los `bounds` del YAML **solo seleccionan posiciones iniciales**. No representan ánodo, cátodo ni paredes. Para COMSOL, el límite de transporte es el dominio numérico muestreado del mapa y las zonas cuyo punto más cercano tiene campo inválido o nulo. No se aplica el cono de ejemplo a TPC_Avanzada. La consulta del campo sigue siendo por vecino más próximo, como en SENTINEL: no es una reconstrucción de superficies ni una máscara de materiales. Los NaN no se rellenan.

El mapa disponible tiene 40000 muestras; 36916 contienen componentes finitas. Sus coordenadas cubren aproximadamente r = 0.1425–56.8575 mm y z = −1.74–101.74 mm. La existencia de campo no identifica LXe: puede incluir PTFE. Selecciona límites adecuados; no se garantiza excluir materiales. Los endpoints en límites numéricos no constituyen una clasificación física de colección. Un inicio inválido se rechaza; alcanzar 2000000 pasos produce un error explícito. El paso del modo vectorial es 0.01 µs, con las parametrizaciones de velocidad y difusión existentes.

Las líneas de campo de ejemplo siguen usando el lector original y las convenciones internas de [conf/test.toml](conf/test.toml), incluido el cono original, conversiones e inversión de z. Para otros archivos de líneas de campo hay que revisar explícitamente ese TOML. Los binarios originales `endpoints`, `trajectories`, `detsim` y `waveforms` siguen disponibles; consulta [su referencia](docs/original-cli.md).

El muestreo con phi fijo se limita al caso axisimétrico inicial. No incluye hilos, geometrías realistas, máscaras de materiales ni análisis físico posterior.

## Comprobaciones

```bash
cargo test --release
PYTHONPATH=python /home/colina/.venv/bin/python -m unittest discover -s python/tests -v
just --justfile dev/justfile test
```

Las pruebas de Python comprueban la rejilla, unidades, IDs, los ocho estados de los booleanos, coherencia endpoints/trayectorias, difusión en 3D, HDF5 y una simulación pequeña con el campo real.
