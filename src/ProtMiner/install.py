import sys, shutil, subprocess as sp

REQUIRED_PY = ["pandas", "Bio"]  # Bio from biopython

PKG_MAP = {
    "hmmsearch": {"conda":"hmmer", "apt":"hmmer", "brew":"hmmer"},
    "diamond":   {"conda":"diamond", "apt":"diamond-aligner", "brew":"diamond"},
    "interproscan.sh": {"conda":"interproscan"},
    "seqkit":    {"conda":"seqkit", "apt":"seqkit", "brew":"seqkit"},
    "tmhmm":     {},   # license-restricted
    "signalp":   {},   # license-restricted
    "deeptmhmm": {"conda":"deeptmhmm"}
}

def _pip_install(pkg: str) -> bool:
    try:
        sp.check_call([sys.executable, "-m", "pip", "install", "--upgrade", pkg])
        return True
    except Exception:
        return False

def ensure_python_modules(auto_install: bool):
    missing = []
    try:
        import pandas  # noqa
    except Exception:
        missing.append("pandas")
    try:
        import Bio.SeqIO  # noqa
    except Exception:
        missing.append("Bio")
    if missing and not auto_install:
        sys.stderr.write(
            f"[WARN] Missing Python modules: {', '.join(missing)}. "
            "Run with --auto-install-py to install via pip.\n"
        )
        return
    if missing and auto_install:
        for m in missing:
            pkg = "biopython" if m == "Bio" else m
            sys.stderr.write(f"[INFO] Installing Python package: {pkg}\n")
            if not _pip_install(pkg):
                sys.stderr.write(f"[ERROR] Failed to install {pkg}. Please install manually.\n")

def which_pkgmgr():
    for mgr in ["conda", "apt", "brew"]:
        if shutil.which(mgr):
            return mgr
    return None

def print_manual_hint(tool: str):
    if tool == "interproscan.sh":
        sys.stderr.write("InterProScan (large):\n  conda install -c bioconda -c conda-forge interproscan\n  https://interpro-docs.ebi.ac.uk/downloads/\n")
    elif tool == "tmhmm":
        sys.stderr.write("TMHMM requires registration; prefer DeepTMHMM via conda.\n  https://services.healthtech.dtu.dk/services/TMHMM-2.0/\n")
    elif tool == "signalp":
        sys.stderr.write("SignalP (license-restricted):\n  https://services.healthtech.dtu.dk/services/SignalP-6.0/\n")
    elif tool == "hmmsearch":
        sys.stderr.write("Install HMMER: conda install -c bioconda -c conda-forge hmmer | brew install hmmer | sudo apt install hmmer\n")
    elif tool == "diamond":
        sys.stderr.write("Install DIAMOND: conda install -c bioconda -c conda-forge diamond | brew install diamond | sudo apt install diamond-aligner\n")
    elif tool == "seqkit":
        sys.stderr.write("Install SeqKit: conda install -c bioconda -c conda-forge seqkit | brew install seqkit | sudo apt install seqkit\n")
    elif tool == "deeptmhmm":
        sys.stderr.write("Install DeepTMHMM: conda install -c bioconda -c conda-forge deeptmhmm\n")
    else:
        sys.stderr.write(f"Install '{tool}' via your package manager.\n")

def try_install(tool: str, auto_install_tools: bool):
    from .utils import have
    if have(tool):
        return
    sys.stderr.write(f"[INFO] '{tool}' not found.\n")
    if not auto_install_tools:
        print_manual_hint(tool); return
    mgr = which_pkgmgr()
    if mgr is None:
        sys.stderr.write("[WARN] No supported package manager found (conda/apt/brew).\n")
        print_manual_hint(tool); return
    pkg = PKG_MAP.get(tool, {}).get(mgr)
    if not pkg:
        print_manual_hint(tool); return
    try:
        if mgr == "conda":
            sp.check_call([mgr, "install", "-y", "-c", "bioconda", "-c", "conda-forge", pkg])
        elif mgr == "apt":
            sp.check_call("sudo apt update && sudo apt install -y " + pkg, shell=True)
        elif mgr == "brew":
            sp.check_call(["brew", "install", pkg])
    except Exception as e:
        sys.stderr.write(f"[WARN] Auto-install failed for {tool} ({e}).\n")
        print_manual_hint(tool)
