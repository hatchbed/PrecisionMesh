#!/usr/bin/env python3
"""Phase C5 corpus harness: A/B-validate the UV-grid CDT tessellation path.

For every STEP part in test/step/ this runs precision_mesh twice with identical
parameters -- once with the CDT path enabled (default) and once with
--no-uv-tess (baseline) -- collecting the --validate-report YAML from each run,
then compares the pairs and writes a summary (markdown + CSV) with a per-part
verdict.

Usage:
    python3 test/validate_cdt.py [--only GLOB] [--exclude GLOB] [--timeout S]
                                 [--surface-error] [--binary PATH] [-- extra args...]

Exits nonzero if any part regresses (FAIL) or the CDT run errors.
"""

import argparse
import csv
import fnmatch
import re
import subprocess
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
STEP_DIR = REPO_ROOT / "test" / "step"
RESULTS_DIR = REPO_ROOT / "test" / "results"
DEFAULT_BINARY = REPO_ROOT / "builddir" / "precision_mesh"

STEP_EXTENSIONS = {".step", ".stp"}

# Warnings from the CDT path that indicate a hard failure even when the run
# completes (see uv_grid_retessellate in src/step_tessellation.cpp).  The
# border-halfedge audit is intentionally absent: its loop-size accounting
# mismatches benignly, and real gaps show up in open_boundary_edges anyway.
HARD_FAILURE_PATTERNS = [
    "add_face rejected",
    "CDT exception",
    "CDT unknown exception",
    "0 faces added",
    "CDT produced 0 interior triangles",
    "no mappable interior triangles",
]


def parse_yaml(path):
    """Minimal flat YAML reader for the --validate-report format.

    Handles nested maps by indentation (2 spaces per level) and scalar values:
    quoted strings, ints, floats, bools, null.  Not a general YAML parser.
    """
    root = {}
    stack = [(0, root)]  # (indent, dict)
    for raw in path.read_text().splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        indent = len(raw) - len(raw.lstrip())
        key, _, value = raw.strip().partition(":")
        value = value.strip()
        while stack and stack[-1][0] > indent:
            stack.pop()
        parent = stack[-1][1]
        if not value:
            child = {}
            parent[key] = child
            stack.append((indent + 2, child))
            continue
        if value.startswith('"') and value.endswith('"'):
            parent[key] = value[1:-1]
        elif value == "null":
            parent[key] = None
        elif value in ("true", "false"):
            parent[key] = value == "true"
        else:
            try:
                parent[key] = int(value)
            except ValueError:
                try:
                    parent[key] = float(value)
                except ValueError:
                    parent[key] = value
    return root


def run_once(binary, part, report_path, log_path, timeout, extra_args, no_uv_tess):
    cmd = [str(binary), "-i", str(part), "--validate-report", str(report_path)]
    if no_uv_tess:
        cmd.append("--no-uv-tess")
    cmd += extra_args
    start = time.monotonic()
    result = {"cmd": " ".join(cmd), "status": "ok", "elapsed": 0.0,
              "report": None, "hard_failures": []}
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              timeout=timeout, text=True)
        log_path.write_text(proc.stdout)
        if proc.returncode != 0:
            result["status"] = f"exit {proc.returncode}"
    except subprocess.TimeoutExpired as e:
        # TimeoutExpired.stdout is bytes even when text=True.
        out = e.stdout or b""
        log_path.write_text(out.decode(errors="replace") if isinstance(out, bytes) else out)
        result["status"] = "timeout"
    except OSError as e:
        result["status"] = f"error: {e}"
    result["elapsed"] = time.monotonic() - start

    if log_path.exists():
        log = log_path.read_text()
        result["hard_failures"] = [p for p in HARD_FAILURE_PATTERNS if p in log]
    if result["status"] == "ok" and report_path.exists():
        result["report"] = parse_yaml(report_path)
    elif result["status"] == "ok":
        result["status"] = "no report"
    return result


def judge(cdt, base):
    """Compare a CDT run against its --no-uv-tess baseline.  Returns (verdict, reasons)."""
    reasons = []

    if cdt["status"] != "ok":
        return "ERROR", [f"CDT run: {cdt['status']}"]
    if base["status"] != "ok":
        return "BASELINE-ERROR", [f"baseline run: {base['status']}"]

    cv = cdt["report"]["validation"]
    bv = base["report"]["validation"]
    ct = cdt["report"]["tessellation"]
    verdict = "PASS"

    def fail(msg):
        nonlocal verdict
        verdict = "FAIL"
        reasons.append(msg)

    def warn(msg):
        nonlocal verdict
        if verdict == "PASS":
            verdict = "WARN"
        reasons.append(msg)

    if cv["open_boundary_edges"] > bv["open_boundary_edges"]:
        fail(f"open edges {cv['open_boundary_edges']} > baseline {bv['open_boundary_edges']}")
    if cv["non_manifold_edges"] > bv["non_manifold_edges"]:
        fail(f"non-manifold {cv['non_manifold_edges']} > baseline {bv['non_manifold_edges']}")
    if cv["border_beyond_tol"] > bv["border_beyond_tol"]:
        fail(f"border verts beyond tol {cv['border_beyond_tol']} > {bv['border_beyond_tol']}")
    if cv["interior_beyond_tol"] > bv["interior_beyond_tol"]:
        fail(f"interior verts beyond tol {cv['interior_beyond_tol']} > {bv['interior_beyond_tol']}")
    if cv.get("flipped_triangles", 0) > bv.get("flipped_triangles", 0):
        fail(f"flipped triangles {cv.get('flipped_triangles', 0)} > "
             f"baseline {bv.get('flipped_triangles', 0)}")
    if cdt["hard_failures"]:
        fail("CDT warnings: " + "; ".join(cdt["hard_failures"]))

    for key, label in (("max_border_edge_dist", "border->edge"),
                       ("max_interior_face_dist", "interior->face")):
        c, b = cv[key], bv[key]
        if c > max(1e-6, 1.5 * b):
            warn(f"{label} max {c:.3g} vs baseline {b:.3g}")

    cse, bse = cv.get("surface_error"), bv.get("surface_error")
    if isinstance(cse, dict) and isinstance(bse, dict):
        if cse["max"] > 1.1 * bse["max"] + 1e-9:
            warn(f"surface error max {cse['max']:.3g} vs baseline {bse['max']:.3g}")

    if ct["cdt_faces_attempted"] > 0 and ct["cdt_faces_succeeded"] == 0:
        warn("CDT attempted but never succeeded (all fell back)")
    elif ct["cdt_faces_attempted"] == 0:
        reasons.append("no CDT-eligible faces")

    return verdict, reasons


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--only", default=None, help="only parts matching this glob")
    ap.add_argument("--exclude", default=None, help="skip parts matching this glob")
    ap.add_argument("--timeout", type=float, default=900.0, help="per-run timeout (s)")
    ap.add_argument("--surface-error", action="store_true",
                    help="also measure triangle-sample surface error (expensive)")
    ap.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    ap.add_argument("--results-dir", type=Path, default=RESULTS_DIR)
    ap.add_argument("extra", nargs="*",
                    help="extra precision_mesh args after '--' (e.g. --max-edge-length-percent 1)")
    args = ap.parse_args()

    if not args.binary.exists():
        sys.exit(f"binary not found: {args.binary} (build first)")
    parts = sorted(p for p in STEP_DIR.iterdir() if p.suffix.lower() in STEP_EXTENSIONS)
    if args.only:
        parts = [p for p in parts if fnmatch.fnmatch(p.name, args.only)]
    if args.exclude:
        parts = [p for p in parts if not fnmatch.fnmatch(p.name, args.exclude)]
    if not parts:
        sys.exit(f"no STEP parts matched in {STEP_DIR}")

    extra_args = list(args.extra)
    if args.surface_error:
        extra_args.append("--validate-surface-error")

    args.results_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for part in parts:
        slug = re.sub(r"[^A-Za-z0-9._-]+", "_", part.stem)
        print(f"=== {part.name}", flush=True)
        runs = {}
        for tag, no_uv in (("cdt", False), ("nouv", True)):
            report = args.results_dir / f"{slug}__{tag}.yaml"
            log = args.results_dir / f"{slug}__{tag}.log"
            report.unlink(missing_ok=True)
            runs[tag] = run_once(args.binary, part, report, log,
                                 args.timeout, extra_args, no_uv)
            print(f"    {tag:>4}: {runs[tag]['status']} ({runs[tag]['elapsed']:.1f}s)",
                  flush=True)
        verdict, reasons = judge(runs["cdt"], runs["nouv"])
        print(f"    => {verdict}" + (f": {'; '.join(reasons)}" if reasons else ""), flush=True)
        rows.append((part.name, runs["cdt"], runs["nouv"], verdict, reasons))

    write_summary(args.results_dir, rows, args.surface_error)
    print(f"\nsummary written to {args.results_dir / 'c5_summary.md'}")
    bad = [r for r in rows if r[3] in ("FAIL", "ERROR")]
    if bad:
        print(f"{len(bad)} part(s) FAILED: " + ", ".join(r[0] for r in bad))
        sys.exit(1)
    print(f"all {len(rows)} part(s) passed"
          + (" (with warnings)" if any(r[3] != "PASS" for r in rows) else ""))


def metric(run, section, key, fmt="{}"):
    if run["report"] is None:
        return "-"
    value = run["report"].get(section, {}).get(key)
    if value is None or isinstance(value, dict):
        return "-"
    return fmt.format(value)


def write_summary(results_dir, rows, surface_error):
    md = ["# C5 cylinder CDT corpus validation", "",
          "Per part: CDT path (default) vs `--no-uv-tess` baseline.", "",
          "| part | verdict | open edges (cdt/base) | non-manifold (cdt/base) "
          "| border max (cdt/base) | interior max (cdt/base) | cdt ok/try "
          + ("| surf err max (cdt/base) " if surface_error else "")
          + "| time s (cdt/base) | notes |",
          "|" + "---|" * (10 if surface_error else 9)]
    csv_rows = []
    for name, cdt, base, verdict, reasons in rows:
        pair = lambda sec, key, fmt="{}": (metric(cdt, sec, key, fmt) + " / "
                                           + metric(base, sec, key, fmt))
        surf = ""
        if surface_error:
            def se(run):
                v = (run["report"] or {}).get("validation", {}).get("surface_error")
                return f"{v['max']:.3g}" if isinstance(v, dict) else "-"
            surf = f"| {se(cdt)} / {se(base)} "
        md.append(
            f"| {name} | {verdict} "
            f"| {pair('validation', 'open_boundary_edges')} "
            f"| {pair('validation', 'non_manifold_edges')} "
            f"| {pair('validation', 'max_border_edge_dist', '{:.3g}')} "
            f"| {pair('validation', 'max_interior_face_dist', '{:.3g}')} "
            f"| {metric(cdt, 'tessellation', 'cdt_faces_succeeded')}/"
            f"{metric(cdt, 'tessellation', 'cdt_faces_attempted')} "
            f"{surf}"
            f"| {cdt['elapsed']:.0f} / {base['elapsed']:.0f} "
            f"| {'; '.join(reasons)} |")
        csv_rows.append({
            "part": name, "verdict": verdict,
            "cdt_status": cdt["status"], "base_status": base["status"],
            "cdt_open_edges": metric(cdt, "validation", "open_boundary_edges"),
            "base_open_edges": metric(base, "validation", "open_boundary_edges"),
            "cdt_non_manifold": metric(cdt, "validation", "non_manifold_edges"),
            "base_non_manifold": metric(base, "validation", "non_manifold_edges"),
            "cdt_border_max": metric(cdt, "validation", "max_border_edge_dist"),
            "base_border_max": metric(base, "validation", "max_border_edge_dist"),
            "cdt_interior_max": metric(cdt, "validation", "max_interior_face_dist"),
            "base_interior_max": metric(base, "validation", "max_interior_face_dist"),
            "cdt_border_beyond_tol": metric(cdt, "validation", "border_beyond_tol"),
            "cdt_interior_beyond_tol": metric(cdt, "validation", "interior_beyond_tol"),
            "cdt_flipped_triangles": metric(cdt, "validation", "flipped_triangles"),
            "base_flipped_triangles": metric(base, "validation", "flipped_triangles"),
            "cdt_faces_attempted": metric(cdt, "tessellation", "cdt_faces_attempted"),
            "cdt_faces_succeeded": metric(cdt, "tessellation", "cdt_faces_succeeded"),
            "cdt_runtime_s": f"{cdt['elapsed']:.1f}",
            "base_runtime_s": f"{base['elapsed']:.1f}",
            "notes": "; ".join(reasons),
        })
    (results_dir / "c5_summary.md").write_text("\n".join(md) + "\n")
    with open(results_dir / "c5_summary.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(csv_rows[0].keys()))
        writer.writeheader()
        writer.writerows(csv_rows)


if __name__ == "__main__":
    main()
