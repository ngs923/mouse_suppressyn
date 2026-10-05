#!/usr/bin/env python3
"""
run_eve_rnaseq.py
=================

Pipeline for quantifying the expression of endogenous retrovirus (ERV)-derived
env genes in mouse placenta RNA-seq data.

  SRA download  ->  Cutadapt  ->  TopHat2  ->  Cufflinks (gEVE GTF)

This reproduces the following procedure:
  Mouse placenta RNA-seq datasets (SRR392616, SRR392617) were downloaded from
  the NCBI SRA, adapter sequences and low-quality bases were trimmed with
  Cutadapt v1.8.3, the processed reads were aligned to GRCm38.p1 with
  TopHat2 v2.1.1, and transcript abundance was quantified with Cufflinks v2.2.1
  using the gEVE annotation Mmus38.geve.m_v1.gtf.

Examples:
  python3 run_eve_rnaseq.py check
  python3 run_eve_rnaseq.py inspect-gtf --geve-gtf ref/Mmus38.geve.m_v1.gtf
  python3 run_eve_rnaseq.py all --config config.json --threads 16

Steps can also be run individually:
  python3 run_eve_rnaseq.py fetch
  python3 run_eve_rnaseq.py trim
  python3 run_eve_rnaseq.py index
  python3 run_eve_rnaseq.py align
  python3 run_eve_rnaseq.py quant
  python3 run_eve_rnaseq.py summarize
"""

from __future__ import annotations

import argparse
import collections
import gzip
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import textwrap
import time
from pathlib import Path

LOG = logging.getLogger("eve")

# --------------------------------------------------------------------------
# Default settings. Override them with a JSON file (see --config).
# --------------------------------------------------------------------------
DEFAULTS = {
    # SRA run accessions to analyse
    "srr_ids": ["SRR392616", "SRR392617"],

    # Working directory
    "outdir": "eve_rnaseq",

    # Reference genome: an uncompressed GRCm38.p1 FASTA file
    "genome_fasta": "ref/GRCm38.p1.fa",

    # Bowtie2 index prefix. If empty, it is derived from genome_fasta.
    # e.g. "ref/GRCm38.p1" -> ref/GRCm38.p1.1.bt2 and friends
    "bowtie2_index": "",

    # gEVE protein-coding EVE annotation.
    # Download it from http://geve.med.u-tokai.ac.jp
    "geve_gtf": "ref/Mmus38.geve.m_v1.gtf",

    # Cutadapt: Illumina TruSeq adapters (change these to match your data)
    "adapter_r1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    "adapter_r2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
    "quality_cutoff": 20,      # cutadapt -q
    "min_length": 30,          # cutadapt -m

    # TopHat2
    "tophat_max_multihits": 20,     # -g : matters a lot for repetitive ERVs
    "tophat_library_type": "fr-unstranded",
    "tophat_mate_inner_dist": 50,   # -r, paired-end only
    "tophat_mate_std_dev": 50,      # --mate-std-dev, paired-end only

    # Cufflinks
    "cufflinks_multi_read_correct": True,   # -u
    "cufflinks_frag_bias_correct": True,    # -b (uses genome_fasta)
    "cufflinks_library_type": "fr-unstranded",

    # summarize: regex used to pick out env genes (matched against the
    # gene_id and the raw attribute field)
    "env_pattern": r"env",

    # Resources
    "threads": 8,
    "sra_tmpdir": "",   # empty means the system default
}

TOOLS = {
    "prefetch":      "sra-tools  (conda install -c bioconda sra-tools)",
    "fasterq-dump":  "sra-tools",
    "pigz":          "pigz (optional; gzip is used when absent)",
    "cutadapt":      "cutadapt v1.8.3 (conda install -c bioconda cutadapt=1.8.3)",
    "bowtie2":       "bowtie2 (required by TopHat2)",
    "bowtie2-build": "bowtie2",
    "tophat2":       "tophat v2.1.1 (conda install -c bioconda tophat=2.1.1); needs python2",
    "cufflinks":     "cufflinks v2.2.1 (conda install -c bioconda cufflinks=2.2.1)",
    "samtools":      "samtools",
}

OPTIONAL_TOOLS = {"pigz", "prefetch", "samtools"}


# --------------------------------------------------------------------------
# Infrastructure
# --------------------------------------------------------------------------
class PipelineError(RuntimeError):
    pass


def setup_logging(logfile: Path | None, verbose: bool) -> None:
    LOG.setLevel(logging.DEBUG if verbose else logging.INFO)
    fmt = logging.Formatter("[%(asctime)s] %(levelname)-7s %(message)s", "%Y-%m-%d %H:%M:%S")
    sh = logging.StreamHandler(sys.stderr)
    sh.setFormatter(fmt)
    LOG.addHandler(sh)
    if logfile:
        logfile.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(logfile, encoding="utf-8")
        fh.setFormatter(fmt)
        LOG.addHandler(fh)


def run(cmd: list[str], *, dry_run: bool = False, cwd: Path | None = None,
        stdout_path: Path | None = None) -> None:
    """Run an external command, raising PipelineError if it fails."""
    pretty = " ".join(str(c) for c in cmd)
    if stdout_path:
        pretty += f"  > {stdout_path}"
    LOG.info("RUN  %s", pretty)
    if dry_run:
        return

    t0 = time.time()
    out = open(stdout_path, "wb") if stdout_path else None
    try:
        proc = subprocess.run([str(c) for c in cmd], cwd=str(cwd) if cwd else None,
                              stdout=out, check=False)
    finally:
        if out:
            out.close()

    if proc.returncode != 0:
        raise PipelineError(f"Command failed with exit code {proc.returncode}:\n  {pretty}")
    LOG.debug("     done (%.1f s)", time.time() - t0)


def tool_path(name: str) -> str | None:
    return shutil.which(name)


def require(name: str) -> str:
    p = tool_path(name)
    if not p:
        raise PipelineError(
            f"Required tool '{name}' is not on PATH. {TOOLS.get(name, '')}")
    return p


def done_marker(path: Path) -> Path:
    return path.parent / f".{path.name}.done"


def is_done(marker: Path, force: bool) -> bool:
    if force:
        return False
    return marker.exists()


def mark_done(marker: Path, dry_run: bool) -> None:
    if not dry_run:
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.write_text(time.strftime("%Y-%m-%d %H:%M:%S\n"))


# --------------------------------------------------------------------------
# Configuration
# --------------------------------------------------------------------------
class Config:
    def __init__(self, data: dict):
        self.__dict__.update(data)
        self.outdir = Path(self.outdir)
        self.genome_fasta = Path(self.genome_fasta)
        self.geve_gtf = Path(self.geve_gtf)

    # Directory layout
    @property
    def dir_sra(self) -> Path:      return self.outdir / "00_sra"
    @property
    def dir_fastq(self) -> Path:    return self.outdir / "01_fastq"
    @property
    def dir_trim(self) -> Path:     return self.outdir / "02_trimmed"
    @property
    def dir_align(self) -> Path:    return self.outdir / "03_tophat2"
    @property
    def dir_quant(self) -> Path:    return self.outdir / "04_cufflinks"
    @property
    def dir_result(self) -> Path:   return self.outdir / "05_results"
    @property
    def dir_log(self) -> Path:      return self.outdir / "logs"

    def index_prefix(self) -> Path:
        if self.bowtie2_index:
            return Path(self.bowtie2_index)
        # genome.fa -> genome
        return self.genome_fasta.with_suffix("")

    def dump(self) -> str:
        d = {k: (str(v) if isinstance(v, Path) else v)
             for k, v in self.__dict__.items()}
        return json.dumps(d, indent=2, ensure_ascii=False)


def load_config(path: str | None, overrides: dict) -> Config:
    data = dict(DEFAULTS)
    if path:
        with open(path, encoding="utf-8") as fh:
            user = json.load(fh)
        unknown = set(user) - set(DEFAULTS)
        if unknown:
            LOG.warning("Unknown keys in config (ignored): %s", ", ".join(sorted(unknown)))
        data.update({k: v for k, v in user.items() if k in DEFAULTS})
    data.update({k: v for k, v in overrides.items() if v is not None})
    return Config(data)


# --------------------------------------------------------------------------
# step: check
# --------------------------------------------------------------------------
def step_check(cfg: Config, args) -> None:
    """Verify that the required tools and reference files are available."""
    LOG.info("=== Checking dependencies ===")
    missing_required, missing_optional = [], []
    for name, hint in TOOLS.items():
        p = tool_path(name)
        if p:
            LOG.info("  [OK]       %-14s %s", name, p)
        elif name in OPTIONAL_TOOLS:
            missing_optional.append((name, hint))
            LOG.warning("  [optional] %-14s not found  -- %s", name, hint)
        else:
            missing_required.append((name, hint))
            LOG.error("  [MISSING]  %-14s not found  -- %s", name, hint)

    # TopHat2 depends on python2, so warn about it separately
    if tool_path("tophat2") and not tool_path("python2"):
        LOG.warning("  tophat2 is present but python2 is not. "
                    "TopHat2 2.1.1 requires python2.")

    LOG.info("=== Checking reference files ===")
    for label, path in [("genome FASTA", cfg.genome_fasta),
                        ("gEVE GTF", cfg.geve_gtf)]:
        if path.exists():
            LOG.info("  [OK]       %-14s %s (%.1f MB)", label, path,
                     path.stat().st_size / 1e6)
        else:
            LOG.error("  [MISSING]  %-14s %s", label, path)
            missing_required.append((label, "see the README"))

    idx = cfg.index_prefix()
    if Path(f"{idx}.1.bt2").exists() or Path(f"{idx}.1.bt2l").exists():
        LOG.info("  [OK]       bowtie2 index  %s.*.bt2", idx)
    else:
        LOG.warning("  [TO BUILD] bowtie2 index %s.*.bt2 -- the 'index' step will create it", idx)

    if missing_required:
        raise PipelineError(
            "Some requirements are missing:\n" +
            "\n".join(f"  - {n}: {h}" for n, h in missing_required))
    LOG.info("All checks passed.")


# --------------------------------------------------------------------------
# step: inspect-gtf
# --------------------------------------------------------------------------
def parse_gtf_attributes(attr: str) -> dict:
    out = {}
    for m in re.finditer(r'(\S+)\s+"([^"]*)"', attr):
        out[m.group(1)] = m.group(2)
    return out


def read_fasta_names(path: Path, limit: int = 100000) -> list[str]:
    opener = gzip.open if str(path).endswith(".gz") else open
    names = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                names.append(line[1:].split()[0])
                if len(names) >= limit:
                    break
    return names


def step_inspect_gtf(cfg: Config, args) -> None:
    """Inspect the gEVE GTF and check that its sequence names match the genome.

    TopHat2 and Cufflinks silently return nothing when the chromosome names
    differ by even one character, so always run this before the real analysis.
    """
    if not cfg.geve_gtf.exists():
        raise PipelineError(f"GTF not found: {cfg.geve_gtf}")

    feature_counts = collections.Counter()
    seqnames = collections.Counter()
    attr_keys = collections.Counter()
    env_hits = 0
    env_re = re.compile(cfg.env_pattern, re.IGNORECASE)
    examples = []

    with open(cfg.geve_gtf, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            seqnames[f[0]] += 1
            feature_counts[f[2]] += 1
            attrs = parse_gtf_attributes(f[8])
            attr_keys.update(attrs.keys())
            if env_re.search(f[8]):
                env_hits += 1
                if len(examples) < 5:
                    examples.append(attrs.get("gene_id", f[8][:80]))

    LOG.info("=== %s ===", cfg.geve_gtf)
    LOG.info("Feature types : %s",
             ", ".join(f"{k}={v}" for k, v in feature_counts.most_common()))
    LOG.info("Attribute keys: %s", ", ".join(sorted(attr_keys)))
    LOG.info("Sequence names: %d  e.g. %s", len(seqnames),
             ", ".join(list(seqnames)[:8]))
    LOG.info("Lines matching env pattern '%s': %d", cfg.env_pattern, env_hits)
    for e in examples:
        LOG.info("    e.g. %s", e)
    if env_hits == 0:
        LOG.warning("No lines matched the env pattern. "
                    "Adjust --env-pattern to match your GTF.")

    # Compare chromosome names
    if cfg.genome_fasta.exists():
        fa = set(read_fasta_names(cfg.genome_fasta))
        gtf = set(seqnames)
        shared = fa & gtf
        LOG.info("=== Chromosome name check ===")
        LOG.info("FASTA: %d sequences / GTF: %d sequences / shared: %d",
                 len(fa), len(gtf), len(shared))
        if not shared:
            LOG.error("No sequence names in common! "
                      "A 'chr1' vs '1' naming mismatch is the usual cause.")
            LOG.error("  FASTA e.g.: %s", ", ".join(sorted(fa)[:5]))
            LOG.error("  GTF   e.g.: %s", ", ".join(sorted(gtf)[:5]))
            LOG.error("  Use --fix-chrom-names to write a converted copy of the GTF.")
        else:
            only_gtf = gtf - fa
            if only_gtf:
                LOG.warning("%d sequence names occur only in the GTF (e.g. %s). "
                            "Those regions will not be quantified.",
                            len(only_gtf), ", ".join(sorted(only_gtf)[:5]))
    else:
        LOG.warning("Genome FASTA not found; skipping the chromosome name check.")

    if getattr(args, "fix_chrom_names", False) and cfg.genome_fasta.exists():
        _fix_chrom_names(cfg, set(read_fasta_names(cfg.genome_fasta)), args.dry_run)


def _fix_chrom_names(cfg: Config, fasta_names: set[str], dry_run: bool) -> None:
    """Write a copy of the GTF with 'chr' added to or removed from its
    sequence names so that they match the genome FASTA."""
    out = cfg.geve_gtf.with_name(cfg.geve_gtf.stem + ".renamed.gtf")
    add_chr = any(n.startswith("chr") for n in fasta_names)
    LOG.info("Converting GTF sequence names (%s 'chr') -> %s",
             "adding" if add_chr else "removing", out)
    if dry_run:
        return
    n = 0
    with open(cfg.geve_gtf, encoding="utf-8", errors="replace") as fin, \
         open(out, "w", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            f = line.split("\t", 1)
            if len(f) != 2:
                continue
            name = f[0]
            new = ("chr" + name) if add_chr and not name.startswith("chr") \
                else (name[3:] if (not add_chr and name.startswith("chr")) else name)
            if new in fasta_names:
                n += 1
            fout.write(new + "\t" + f[1])
    LOG.info("Done. %d lines now match the FASTA. "
             "Point geve_gtf in your config at %s.", n, out)


# --------------------------------------------------------------------------
# step: fetch
# --------------------------------------------------------------------------
def fastq_files_for(cfg: Config, srr: str) -> list[Path]:
    """Return the already-downloaded fastq files (two if paired, one if single)."""
    pe = [cfg.dir_fastq / f"{srr}_1.fastq.gz", cfg.dir_fastq / f"{srr}_2.fastq.gz"]
    se = [cfg.dir_fastq / f"{srr}.fastq.gz"]
    if all(p.exists() for p in pe):
        return pe
    if se[0].exists():
        return se
    return []


def step_fetch(cfg: Config, args) -> None:
    """Download fastq files from the SRA. Layout is detected automatically."""
    require("fasterq-dump")
    cfg.dir_fastq.mkdir(parents=True, exist_ok=True)
    cfg.dir_sra.mkdir(parents=True, exist_ok=True)

    compressor = ["pigz", "-p", str(cfg.threads)] if tool_path("pigz") else ["gzip"]

    for srr in cfg.srr_ids:
        marker = done_marker(cfg.dir_fastq / srr)
        if is_done(marker, args.force):
            LOG.info("[skip] %s already downloaded (%s)", srr,
                     ", ".join(p.name for p in fastq_files_for(cfg, srr)))
            continue

        if tool_path("prefetch"):
            run(["prefetch", "--max-size", "100G", "-O", cfg.dir_sra, srr],
                dry_run=args.dry_run)

        cmd = ["fasterq-dump", "--split-3", "--skip-technical",
               "--threads", str(cfg.threads),
               "--outdir", str(cfg.dir_fastq)]
        if cfg.sra_tmpdir:
            cmd += ["--temp", cfg.sra_tmpdir]
        sra_path = cfg.dir_sra / srr
        cmd.append(str(sra_path) if sra_path.exists() else srr)
        run(cmd, dry_run=args.dry_run)

        if not args.dry_run:
            produced = sorted(cfg.dir_fastq.glob(f"{srr}*.fastq"))
            if not produced:
                raise PipelineError(f"{srr}: no fastq files were produced.")
            LOG.info("%s: got %d fastq file(s), compressing", srr, len(produced))
            for p in produced:
                run(compressor + [str(p)], dry_run=False)
        mark_done(marker, args.dry_run)

    if not args.dry_run:
        for srr in cfg.srr_ids:
            files = fastq_files_for(cfg, srr)
            LOG.info("%s: %s (%s)", srr,
                     "paired-end" if len(files) == 2 else "single-end",
                     ", ".join(p.name for p in files))


# --------------------------------------------------------------------------
# step: trim  (Cutadapt)
# --------------------------------------------------------------------------
def trimmed_files_for(cfg: Config, srr: str) -> list[Path]:
    raw = fastq_files_for(cfg, srr)
    if len(raw) == 2:
        return [cfg.dir_trim / f"{srr}_1.trimmed.fastq.gz",
                cfg.dir_trim / f"{srr}_2.trimmed.fastq.gz"]
    return [cfg.dir_trim / f"{srr}.trimmed.fastq.gz"]


def step_trim(cfg: Config, args) -> None:
    """Trim adapter sequences and low-quality bases with Cutadapt."""
    require("cutadapt")
    cfg.dir_trim.mkdir(parents=True, exist_ok=True)
    cfg.dir_log.mkdir(parents=True, exist_ok=True)

    for srr in cfg.srr_ids:
        raw = fastq_files_for(cfg, srr)
        if not raw:
            raise PipelineError(f"{srr}: no fastq found. Run the 'fetch' step first.")
        out = trimmed_files_for(cfg, srr)
        marker = done_marker(cfg.dir_trim / srr)
        if is_done(marker, args.force):
            LOG.info("[skip] %s already trimmed", srr)
            continue

        cmd = ["cutadapt",
               "-q", str(cfg.quality_cutoff),
               "-m", str(cfg.min_length),
               "-a", cfg.adapter_r1]
        if len(raw) == 2:
            cmd += ["-A", cfg.adapter_r2,
                    "-o", str(out[0]), "-p", str(out[1]),
                    str(raw[0]), str(raw[1])]
        else:
            cmd += ["-o", str(out[0]), str(raw[0])]

        run(cmd, dry_run=args.dry_run,
            stdout_path=cfg.dir_log / f"cutadapt.{srr}.log")
        mark_done(marker, args.dry_run)


# --------------------------------------------------------------------------
# step: index  (bowtie2-build)
# --------------------------------------------------------------------------
def step_index(cfg: Config, args) -> None:
    """Build the Bowtie2 index, skipping it if one already exists."""
    idx = cfg.index_prefix()
    if (Path(f"{idx}.1.bt2").exists() or Path(f"{idx}.1.bt2l").exists()) and not args.force:
        LOG.info("[skip] bowtie2 index already exists: %s.*.bt2", idx)
    else:
        require("bowtie2-build")
        if not cfg.genome_fasta.exists():
            raise PipelineError(f"Genome FASTA not found: {cfg.genome_fasta}")
        idx.parent.mkdir(parents=True, exist_ok=True)
        LOG.info("Building the bowtie2 index (1-2 hours for the mouse genome)")
        run(["bowtie2-build", "--threads", str(cfg.threads),
             str(cfg.genome_fasta), str(idx)], dry_run=args.dry_run)

    # TopHat2 expects <index_prefix>.fa alongside the index
    fa_link = Path(f"{idx}.fa")
    if not fa_link.exists() and cfg.genome_fasta.exists() and not args.dry_run:
        LOG.info("Creating %s for TopHat2 (link to %s)", fa_link, cfg.genome_fasta)
        try:
            os.symlink(cfg.genome_fasta.resolve(), fa_link)
        except OSError:
            shutil.copy2(cfg.genome_fasta, fa_link)


# --------------------------------------------------------------------------
# step: align  (TopHat2)
# --------------------------------------------------------------------------
def step_align(cfg: Config, args) -> None:
    """Align reads to the genome with TopHat2."""
    require("tophat2")
    require("bowtie2")
    idx = cfg.index_prefix()
    if not args.dry_run and not (Path(f"{idx}.1.bt2").exists()
                                 or Path(f"{idx}.1.bt2l").exists()):
        raise PipelineError(f"No bowtie2 index at {idx}. Run the 'index' step first.")

    for srr in cfg.srr_ids:
        trimmed = trimmed_files_for(cfg, srr)
        if not all(p.exists() for p in trimmed) and not args.dry_run:
            raise PipelineError(f"{srr}: no trimmed fastq. Run the 'trim' step first.")
        odir = cfg.dir_align / srr
        marker = done_marker(odir)
        if is_done(marker, args.force):
            LOG.info("[skip] %s already aligned (%s)", srr, odir / "accepted_hits.bam")
            continue

        cmd = ["tophat2",
               "-o", str(odir),
               "-p", str(cfg.threads),
               "-g", str(cfg.tophat_max_multihits),
               "--library-type", cfg.tophat_library_type]
        if len(trimmed) == 2:
            cmd += ["-r", str(cfg.tophat_mate_inner_dist),
                    "--mate-std-dev", str(cfg.tophat_mate_std_dev)]
        cmd.append(str(idx))
        cmd += [str(p) for p in trimmed]

        run(cmd, dry_run=args.dry_run)

        bam = odir / "accepted_hits.bam"
        if not args.dry_run:
            if not bam.exists():
                raise PipelineError(f"{srr}: {bam} was not produced.")
            if tool_path("samtools"):
                run(["samtools", "index", str(bam)], dry_run=False)
                run(["samtools", "flagstat", str(bam)], dry_run=False,
                    stdout_path=cfg.dir_log / f"flagstat.{srr}.txt")
        mark_done(marker, args.dry_run)


# --------------------------------------------------------------------------
# step: quant  (Cufflinks)
# --------------------------------------------------------------------------
def step_quant(cfg: Config, args) -> None:
    """Compute FPKM values against the gEVE annotation with Cufflinks.

    Using -G means no novel transcripts are assembled: abundance is estimated
    only for the coordinates in the supplied GTF, which matches the described
    procedure of quantifying with the gEVE annotation.
    """
    require("cufflinks")
    if not cfg.geve_gtf.exists():
        raise PipelineError(f"gEVE GTF not found: {cfg.geve_gtf}")

    for srr in cfg.srr_ids:
        bam = cfg.dir_align / srr / "accepted_hits.bam"
        if not bam.exists() and not args.dry_run:
            raise PipelineError(f"{srr}: {bam} not found. Run the 'align' step first.")
        odir = cfg.dir_quant / srr
        marker = done_marker(odir)
        if is_done(marker, args.force):
            LOG.info("[skip] %s already quantified", srr)
            continue

        cmd = ["cufflinks",
               "-o", str(odir),
               "-p", str(cfg.threads),
               "-G", str(cfg.geve_gtf),
               "--library-type", cfg.cufflinks_library_type]
        if cfg.cufflinks_multi_read_correct:
            cmd.append("-u")        # redistribute multi-mapping reads over repeats
        if cfg.cufflinks_frag_bias_correct and cfg.genome_fasta.exists():
            cmd += ["-b", str(cfg.genome_fasta)]
        cmd.append(str(bam))

        run(cmd, dry_run=args.dry_run,
            stdout_path=cfg.dir_log / f"cufflinks.{srr}.log")
        mark_done(marker, args.dry_run)


# --------------------------------------------------------------------------
# step: summarize
# --------------------------------------------------------------------------
def read_fpkm_tracking(path: Path) -> dict[str, dict]:
    """Read a genes.fpkm_tracking file into {tracking_id: {...}}."""
    rows = {}
    with open(path, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        ix = {name: i for i, name in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            tid = f[ix["tracking_id"]]
            rows[tid] = {
                "locus":  f[ix.get("locus", 0)],
                "fpkm":   float(f[ix["FPKM"]]),
                "lo":     float(f[ix["FPKM_conf_lo"]]),
                "hi":     float(f[ix["FPKM_conf_hi"]]),
                "status": f[ix["FPKM_status"]],
            }
    return rows


def step_summarize(cfg: Config, args) -> None:
    """Combine the per-sample FPKM values and extract the env candidates."""
    cfg.dir_result.mkdir(parents=True, exist_ok=True)

    per_sample = {}
    for srr in cfg.srr_ids:
        p = cfg.dir_quant / srr / "genes.fpkm_tracking"
        if not p.exists():
            raise PipelineError(f"{p} not found. Run the 'quant' step first.")
        per_sample[srr] = read_fpkm_tracking(p)
        LOG.info("%s: read %d genes", srr, len(per_sample[srr]))

    all_ids = sorted(set().union(*(set(d) for d in per_sample.values())))
    env_re = re.compile(cfg.env_pattern, re.IGNORECASE)

    matrix = cfg.dir_result / "geve_fpkm_matrix.tsv"
    env_tab = cfg.dir_result / "geve_env_expressed.tsv"

    if args.dry_run:
        LOG.info("[dry-run] would write %s and %s", matrix, env_tab)
        return

    cols = list(per_sample)
    with open(matrix, "w", encoding="utf-8") as fh:
        fh.write("gene_id\tlocus\t" + "\t".join(f"FPKM_{c}" for c in cols)
                 + "\t" + "\t".join(f"status_{c}" for c in cols)
                 + "\tmean_FPKM\tis_env\n")
        env_rows = []
        for gid in all_ids:
            vals, stats, locus = [], [], ""
            for c in cols:
                r = per_sample[c].get(gid)
                vals.append(r["fpkm"] if r else 0.0)
                stats.append(r["status"] if r else "NA")
                locus = locus or (r["locus"] if r else "")
            mean = sum(vals) / len(vals)
            is_env = bool(env_re.search(gid))
            fh.write(f"{gid}\t{locus}\t"
                     + "\t".join(f"{v:.4f}" for v in vals) + "\t"
                     + "\t".join(stats) + f"\t{mean:.4f}\t{int(is_env)}\n")
            if is_env and mean > 0:
                env_rows.append((gid, locus, vals, stats, mean))

    env_rows.sort(key=lambda r: r[4], reverse=True)
    with open(env_tab, "w", encoding="utf-8") as fh:
        fh.write("gene_id\tlocus\t" + "\t".join(f"FPKM_{c}" for c in cols)
                 + "\t" + "\t".join(f"status_{c}" for c in cols) + "\tmean_FPKM\n")
        for gid, locus, vals, stats, mean in env_rows:
            fh.write(f"{gid}\t{locus}\t"
                     + "\t".join(f"{v:.4f}" for v in vals) + "\t"
                     + "\t".join(stats) + f"\t{mean:.4f}\n")

    LOG.info("All EVEs: %d -> %s", len(all_ids), matrix)
    LOG.info("Expressed env candidates: %d -> %s", len(env_rows), env_tab)
    for gid, locus, vals, stats, mean in env_rows[:15]:
        LOG.info("  %-40s %-28s mean FPKM = %8.2f", gid[:40], locus[:28], mean)
    if not env_rows:
        LOG.warning("No env candidates found. "
                    "Check --env-pattern and the output of inspect-gtf.")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
STEPS = {
    "check":       step_check,
    "inspect-gtf": step_inspect_gtf,
    "fetch":       step_fetch,
    "trim":        step_trim,
    "index":       step_index,
    "align":       step_align,
    "quant":       step_quant,
    "summarize":   step_summarize,
}

ALL_ORDER = ["check", "fetch", "trim", "index", "align", "quant", "summarize"]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description="Quantify ERV-derived env gene expression in mouse placenta RNA-seq",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent("""\
            steps:
              check        verify the required tools and reference files
              inspect-gtf  inspect the gEVE GTF and check chromosome names
              fetch        download fastq files from the SRA
              trim         Cutadapt
              index        bowtie2-build
              align        TopHat2
              quant        Cufflinks
              summarize    combine FPKM values and extract env candidates
              all          run everything from check onwards
        """))
    ap.add_argument("step", choices=list(STEPS) + ["all"])
    ap.add_argument("--config", help="JSON file with the settings")
    ap.add_argument("--outdir")
    ap.add_argument("--genome-fasta")
    ap.add_argument("--geve-gtf")
    ap.add_argument("--bowtie2-index")
    ap.add_argument("--srr-ids", nargs="+", dest="srr_ids")
    ap.add_argument("--env-pattern")
    ap.add_argument("--threads", type=int)
    ap.add_argument("--force", action="store_true",
                    help="redo steps that are already marked as done")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the commands without running them")
    ap.add_argument("--fix-chrom-names", action="store_true",
                    help="inspect-gtf: write a GTF whose sequence names match the FASTA")
    ap.add_argument("--print-config", action="store_true")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)

    overrides = {k: getattr(args, k) for k in
                 ["outdir", "genome_fasta", "geve_gtf", "bowtie2_index",
                  "srr_ids", "env_pattern", "threads"]}
    cfg = load_config(args.config, overrides)

    setup_logging(cfg.dir_log / "pipeline.log", args.verbose)

    if args.print_config:
        print(cfg.dump())
        return 0

    LOG.info("Working directory: %s", cfg.outdir.resolve())
    LOG.info("Runs: %s", ", ".join(cfg.srr_ids))

    steps = ALL_ORDER if args.step == "all" else [args.step]
    try:
        for name in steps:
            LOG.info("")
            LOG.info("########## %s ##########", name)
            STEPS[name](cfg, args)
    except PipelineError as e:
        LOG.error("%s", e)
        return 1
    except KeyboardInterrupt:
        LOG.error("Interrupted.")
        return 130

    LOG.info("")
    LOG.info("Finished. Results: %s", cfg.dir_result)
    return 0


if __name__ == "__main__":
    sys.exit(main())
