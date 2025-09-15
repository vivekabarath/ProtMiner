import shutil, subprocess as sp
from pathlib import Path

def have(tool: str) -> bool:
    return shutil.which(tool) is not None

def which(tool: str):
    return shutil.which(tool)

def run_logged(cmd, log: Path) -> int:
    with log.open("a") as L:
        L.write("\n# CMD: " + " ".join(cmd) + "\n")
        proc = sp.run(cmd, stdout=sp.PIPE, stderr=sp.STDOUT, text=True)
        L.write(proc.stdout)
        return proc.returncode

def ensure_dir(p: Path):
    p.mkdir(parents=True, exist_ok=True)
