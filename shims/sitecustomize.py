"""
Runtime patching for external `repex_tarean` without modifying it.

Python auto-imports `sitecustomize` (if present on PYTHONPATH) during startup.
We use that hook to monkeypatch `lib.r2py.create_connection()` so:
- Rserve runs on a caller-chosen RSERVE_PORT (no port collisions)
- Rserve uses a caller-chosen workdir (no shared /tmp/Rserv collisions)

This file must remain tiny and defensive: it should never break non-seqclust
Python runs.
"""

from __future__ import annotations

import os
import sys
import time
import subprocess


def _log(msg: str) -> None:
    path = os.environ.get("REPORTR_PATCH_LOG")
    if not path:
        return
    try:
        with open(path, "a") as f:
            f.write(msg.rstrip() + "\n")
    except Exception:
        return


def _is_seqclust_process() -> bool:
    argv0 = (sys.argv[0] or "").lower()
    return argv0.endswith("seqclust") or argv0.endswith("/seqclust") or "seqclust" in argv0


def _patch_repex_r2py() -> None:
    # Ensure the seqclust directory is on sys.path so `import lib` works even
    # when seqclust is invoked via an absolute path.
    try:
        seqclust_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
        if seqclust_dir and seqclust_dir not in sys.path:
            sys.path.insert(0, seqclust_dir)
    except Exception:
        pass

    # Import inside function so normal Python runs don't pay import cost.
    from lib import r2py  # type: ignore
    import config  # type: ignore
    import pyRserve  # type: ignore

    original_create = getattr(r2py, "create_connection", None)
    if not callable(original_create):
        _log("reportr sitecustomize: lib.r2py.create_connection not callable; skipping")
        return

    def create_connection_patched():
        env_port = os.environ.get("RSERVE_PORT")
        workdir = os.environ.get("RSERVE_WORKDIR") or os.environ.get("R_SESSION_TMPDIR")

        # If caller didn't request anything, fall back to upstream behavior.
        if not env_port and not workdir:
            return original_create()

        try:
            port = int(env_port) if env_port else None
        except ValueError:
            port = None

        # Upstream uses a random open port; if we want determinism, we must
        # explicitly pick the port.
        if port is None:
            return original_create()

        config.RSERVE_PORT = port
        print("Trying to start Rserve...",)

        cmd = ["R", "CMD", "Rserve", "--RS-port", str(port), "-q", "--no-save"]
        if workdir:
            cmd.extend(["--RS-set", f"workdir={workdir}"])

        # Rserve daemonizes; ignore output here.
        subprocess.run(cmd, check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        time.sleep(1)

        try:
            conn = pyRserve.connect(port=port)
            print("connection OK")
            conn.close()
            import atexit

            atexit.register(r2py.shutdown, port)
            return port
        except Exception:
            print("Connection with Rserve was not established!")
            raise

    r2py.create_connection = create_connection_patched
    _log("reportr sitecustomize: patched lib.r2py.create_connection")


def _main() -> None:
    # Opt-in flag (default on) + only patch actual seqclust processes.
    if os.environ.get("REPORTR_PATCH_RSERVE", "1") != "1":
        return
    if not _is_seqclust_process():
        return

    try:
        _log(f"reportr sitecustomize: argv0={sys.argv[0]!r} python={sys.executable!r}")
        _log(
            "reportr sitecustomize: env "
            f"RSERVE_PORT={os.environ.get('RSERVE_PORT')!r} "
            f"RSERVE_WORKDIR={os.environ.get('RSERVE_WORKDIR')!r}"
        )
        _patch_repex_r2py()
    except Exception as e:
        # Never brick the pipeline due to patching.
        _log(f"reportr sitecustomize: exception while patching; skipping: {e!r}")
        return


_main()

