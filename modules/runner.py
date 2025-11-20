#!/usr/bin/env python3
import re
import shlex
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable, List, Tuple, Optional

# Define patterns to identify R1 and R2 FASTQ files
PAIR_PATTERNS = [
    (re.compile(r"(.+?)_R1(?:[_-]?\d+)?\\.fastq\\.gz$"), "_R2.fastq.gz"),
    (re.compile(r"(.+?)_1(?:[_-]?\d+)?\\.fastq\\.gz$"), "_2.fastq.gz"),
]

def walk_fastqs(root: Path):
    """Recursively finds all FASTQ files in the given root directory."""
    return [p for p in Path(root).rglob("*.fastq.gz")]

# Matches either SAMPLE_R1.fastq.gz / SAMPLE_R2.fastq.gz
# or SAMPLE_1.fastq.gz / SAMPLE_2.fastq.gz, with optional extra run chunk: _R1_001, _2-01, etc.
_R12 = re.compile(r"^(?P<base>.+?)_(?:R)?(?P<read>[12])(?:[_-]?\d+)?\.fastq\.gz$", re.IGNORECASE)

def pair_fastqs(root: Path | str) -> List[Tuple[str, Path, Optional[Path]]]:
    """
    Recursively find FASTQs under `root` and pair R1/R2 files.

    Returns a list of tuples: (sample_name, r1_path, r2_path_or_None)

    Pairing key is (directory, base) so identical basenames in different folders don't collide.
    """
    root = Path(root)

    r1_map: dict[tuple[Path, str], Path] = {}
    r2_map: dict[tuple[Path, str], Path] = {}

    # Walk *.fastq.gz
    for p in sorted(root.rglob("*.fastq.gz")):
        m = _R12.match(p.name)
        if not m:
            # Not an R1/R2 filename we recognize; skip silently
            continue
        base = m.group("base")
        read = m.group("read")
        key = (p.parent, base)

        if read == "1":
            # Only keep the first seen R1 for a key (avoid accidental overwrites)
            r1_map.setdefault(key, p)
        else:  # read == "2"
            r2_map.setdefault(key, p)

    # Build output as (sample_name, r1, r2_or_None)
    pairs: List[Tuple[str, Path, Optional[Path]]] = []
    for key, r1 in sorted(r1_map.items()):
        folder, base = key
        r2 = r2_map.get(key)
        # Sample name: include folder name to disambiguate identical basenames in different dirs
        sample_name = base
        pairs.append((sample_name, r1, r2))

    return pairs


def run_crispresso_parallel(
        *, 
        pairs, 
        amplicon: str, 
        guide: str, 
        coding_seq: str, 
        out_root: Path,
        threads_per_sample: int = 4, 
        parallel_samples: int = 4, 
        extra_args: str = "",
        default_min_aln_score: int = 60,
        plot_window_size: int = 20,
    ) -> Iterable[str]:
    """Run multiple CRISPResso2 analyses in parallel and yield log lines."""
    
    # Ensure output root exists
    out_root = Path(out_root)
    out_root.mkdir(parents=True, exist_ok=True)

    # Inner function that runs CRISPResso for one sample
    def one(sample: str, r1: Path, r2: Path):
        parent = r1.parent
        rep = parent.parent.name if parent.parent.exists() else ""
        if rep.startswith("Replica-"):
            sample_out = out_root / rep / parent.name / sample
        else:
            sample_out = out_root / sample
        sample_out.mkdir(parents=True, exist_ok=True)
        # Build the CRISPResso command
        cmd = [
            "CRISPResso",
            "-r1", str(r1),
            "-a", amplicon,
            "-g", guide,
            "-o", str(sample_out),
            "--n_processes", str(threads_per_sample),
        ]
        # If there’s a paired-end read, include r2
        if r2 is not None:
            cmd.extend(["-r2", str(r2)])
        # If a coding seq is provided, add it
        if coding_seq:
            cmd.extend(["--coding_seq", coding_seq])
        # Add alignment score and plot window size parameters
        cmd.extend([
            "--default_min_aln_score", str(default_min_aln_score),
            "--plot_window_size", str(plot_window_size),
        ])
        # Add any extra arguments
        if extra_args:
            cmd.extend(shlex.split(extra_args))

        # Run the command and capture output
        proc = subprocess.run(
            cmd, 
            stdout=subprocess.PIPE,   # capture output
            stderr=subprocess.STDOUT, # merge stderr into stdout
            text=True                 # get output as string
        )

        # Return sample name, exit code, output log, and output folder
        return sample, proc.returncode, proc.stdout, sample_out
    
    # Use ThreadPoolExecutor to run samples in parallel
    with ThreadPoolExecutor(max_workers=parallel_samples) as ex:
        # Submit one future per sample
        futs = [ex.submit(one, s, r1, r2) for (s, r1, r2) in pairs]
        # Process results as soon as each job finishes
        for fut in as_completed(futs):
            sample, code, _, _ = fut.result()
            if code != 0:
                yield f"[ERROR] {sample} failed."
        yield f"All CRISPResso2 jobs completed ({len(pairs)} samples)."
