from __future__ import annotations
import time, shutil, subprocess, hashlib, re
from pathlib import Path
from typing import List, Optional, Dict

def now_ts() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")

def log(msg: str) -> None:
    print(f"[{now_ts()}] {msg}", flush=True)

class PipelineError(RuntimeError):
    pass

def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p

def which(cmd: str) -> Optional[str]:
    return shutil.which(cmd)

def require_tools(tools: List[str]) -> None:
    missing = [t for t in tools if which(t) is None]
    if missing:
        raise PipelineError("Missing required external tools in PATH: " + ", ".join(missing))

def md5_file(path: Path, chunk_size: int = 2**20) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        while True:
            b = f.read(chunk_size)
            if not b:
                break
            h.update(b)
    return h.hexdigest()

def run_cmd(cmd: List[str], cwd: Optional[Path] = None, log_path: Optional[Path] = None, check: bool = True) -> subprocess.CompletedProcess:
    if log_path:
        ensure_dir(log_path.parent)
        with log_path.open("a", encoding="utf-8") as fh:
            fh.write(f"\n[{now_ts()}] CMD: {' '.join(cmd)}\n")
    log("CMD: " + " ".join(cmd))
    cp = subprocess.run(cmd, cwd=str(cwd) if cwd else None, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if log_path:
        with log_path.open("a", encoding="utf-8") as fh:
            fh.write(cp.stdout)
            fh.write(cp.stderr)
    if check and cp.returncode != 0:
        raise PipelineError(f"Command failed ({cp.returncode}): {' '.join(cmd)}\n{cp.stderr[:2000]}")
    return cp

def safe_split_semi(s: str) -> List[str]:
    if s is None:
        return []
    s = str(s).strip()
    if not s or s.lower() in {"nan","none"}:
        return []
    parts = re.split(r"[;,\s]+", s)
    return [p for p in parts if p]
