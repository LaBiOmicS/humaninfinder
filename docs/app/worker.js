// Web Worker: runs HumaninFinder inside Pyodide, off the page's main thread.
// Messages in:  {type: "run", id, fasta: ArrayBuffer, filename, args: [...]}
// Messages out: {type: "status", text} | {type: "log", text} | {type: "ready", versions}
//               {type: "done", id, files: {name: ArrayBuffer}, seconds} | {type: "error", id, text}
const PYODIDE_VERSION = "0.27.8";
const PYODIDE_URL = `https://cdn.jsdelivr.net/pyodide/v${PYODIDE_VERSION}/full/`;
const WHEELS = [
  "wheels/pyhmmer-0.12.3-cp312-abi3-pyodide_2024_0_wasm32.whl",
  "wheels/humanin_finder-1.0.0b3-py3-none-any.whl",
];

importScripts(PYODIDE_URL + "pyodide.js");

const post = (msg, transfer) => self.postMessage(msg, transfer || []);
let pyodideReady = null;

async function init() {
  post({ type: "status", text: "Loading Python (Pyodide " + PYODIDE_VERSION + ")…" });
  const py = await loadPyodide({
    indexURL: PYODIDE_URL,
    stdout: (s) => post({ type: "log", text: s }),
    stderr: (s) => post({ type: "log", text: s }),
  });
  post({ type: "status", text: "Loading Biopython and click…" });
  await py.loadPackage(["micropip", "biopython", "click"]);

  post({ type: "status", text: "Loading pyhmmer and HumaninFinder…" });
  const base = new URL(".", self.location.href).href;
  // psutil is not packaged for Pyodide; pyhmmer only calls cpu_count() when cpus=0
  py.FS.mkdirTree("/stubs");
  py.FS.writeFile("/stubs/psutil.py", "def cpu_count(logical=True):\n    return 1\n");
  py.globals.set("wheel_urls", py.toPy(WHEELS.map((w) => base + w)));
  await py.runPythonAsync(`
import sys
sys.path.insert(0, "/stubs")
import micropip
await micropip.install(list(wheel_urls), deps=False)
import pyhmmer, humanin_finder
versions = {"humanin_finder": humanin_finder.__version__, "pyhmmer": pyhmmer.__version__}
`);
  const versions = py.globals.get("versions").toJs({ dict_converter: Object.fromEntries });
  post({ type: "ready", versions });
  return py;
}

async function run(py, msg) {
  const t0 = performance.now();
  const work = "/work";
  // fresh working directory for every run
  py.runPython(`
import os, shutil
os.chdir("/")  # the previous run left us inside /work, which cannot be removed while current
shutil.rmtree("${work}", ignore_errors=True)
os.makedirs("${work}")
os.chdir("${work}")
`);
  const inputName = "/work/" + msg.filename.replace(/[^A-Za-z0-9._-]/g, "_");
  py.FS.writeFile(inputName, new Uint8Array(msg.fasta));
  const args = msg.args.map((a) => (a === "@INPUT" ? inputName : a));
  py.globals.set("cli_args", py.toPy(args));
  await py.runPythonAsync(`
from humanin_finder.cli import main
try:
    main.main(args=list(cli_args), prog_name="humanin-finder", standalone_mode=False)
except SystemExit as e:
    if e.code not in (0, None):
        raise RuntimeError(f"HumaninFinder exited with status {e.code}")
`);
  const files = {};
  const transfer = [];
  for (const f of py.FS.readdir(work)) {
    if (f === "." || f === ".." || "/work/" + f === inputName) continue;
    const data = py.FS.readFile(work + "/" + f);
    files[f] = data.buffer;
    transfer.push(data.buffer);
  }
  post({ type: "done", id: msg.id, files, seconds: (performance.now() - t0) / 1000 }, transfer);
}

self.onmessage = async (event) => {
  const msg = event.data;
  if (msg.type !== "run" && msg.type !== "init") return;
  try {
    if (!pyodideReady) pyodideReady = init();
    const py = await pyodideReady;
    if (msg.type === "run") await run(py, msg);
  } catch (err) {
    const text = String((err && err.message) || err);
    post({ type: "error", id: msg.id, text: text.split("\n").slice(-8).join("\n") });
  }
};
