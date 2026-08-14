#!/usr/bin/env python3
"""Rebuild Wave containers after a GINFINITY release and pin them in the modules.

Updates every `environment.yml` / `environment.gpu.yml`, Wave-freezes the CPU
and GPU images (Docker + Singularity), then writes the new URLs into the
`BUILD_RNA_GRAPHS` and `EMBED_RNA_GRAPHS` `container` directives.

Example:

    python3 scripts/bump_ginfinity_containers.py --version 1.0.2

CPU builds take a few minutes; GPU builds often take 10–20 minutes each.
Requires the `wave` CLI on PATH (https://seqera.io/wave/).
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
CHANNEL = "nicolas.aira"
PACKAGE = "ginfinity"

BUILD_ENV = REPO_ROOT / "modules/build_rna_graphs/environment.yml"
BUILD_MAIN = REPO_ROOT / "modules/build_rna_graphs/main.nf"
BUILD_META = REPO_ROOT / "modules/build_rna_graphs/meta.yml"
EMBED_ENV = REPO_ROOT / "modules/embed_rna_graphs/environment.yml"
EMBED_GPU_ENV = REPO_ROOT / "modules/embed_rna_graphs/environment.gpu.yml"
EMBED_MAIN = REPO_ROOT / "modules/embed_rna_graphs/main.nf"
AGENTS = REPO_ROOT / "AGENTS.md"

GIN_DEP_RE = re.compile(r"^(\s*-\s*nicolas\.aira::ginfinity=)\S+\s*$", re.MULTILINE)
PYTORCH_DEP_RE = re.compile(r"^(\s*-\s*conda-forge::pytorch-gpu=)\S+\s*$", re.MULTILINE)
CUDA_DEP_RE = re.compile(r"^(\s*-\s*conda-forge::cuda-version=)\S+\s*$", re.MULTILINE)

BUILD_CONTAINER_RE = re.compile(
    r"(container \"\$\{ workflow\.containerEngine == 'singularity' "
    r"&& !task\.ext\.singularity_pull_docker_container \?\n"
    r"        ')([^']+)(' :\n        ')([^']+)(' \}\")"
)
EMBED_CONTAINER_RE = re.compile(
    r"(container \"\$\{ workflow\.containerEngine in \['singularity', 'apptainer'\] "
    r"&& !task\.ext\.singularity_pull_docker_container \?\n"
    r"        \(task\.accelerator \? ')([^']+)(' : ')([^']+)('\) :\n"
    r"        \(task\.accelerator \? ')([^']+)(' : ')([^']+)('\) \}\")"
)

VERSION_RE = re.compile(r"^[0-9]+(?:\.[0-9]+)*(?:[a-zA-Z0-9._+-]*)?$")


@dataclass(frozen=True)
class WaveImage:
    container_image: str
    build_id: str
    cached: bool

    @property
    def docker_ref(self) -> str:
        image = self.container_image.strip()
        for prefix in ("https://", "http://", "docker://"):
            if image.startswith(prefix):
                image = image[len(prefix) :]
        return image


@dataclass(frozen=True)
class Pins:
    cpu_docker: str
    cpu_sif: str
    cpu_docker_build_id: str
    cpu_sif_build_id: str
    cpu_oras: str
    gpu_docker: str | None = None
    gpu_sif: str | None = None
    gpu_docker_build_id: str | None = None
    gpu_sif_build_id: str | None = None
    gpu_oras: str | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--version",
        required=True,
        help="New nicolas.aira::ginfinity version, e.g. 1.0.2",
    )
    parser.add_argument(
        "--pytorch-gpu",
        help="Override conda-forge::pytorch-gpu pin in environment.gpu.yml",
    )
    parser.add_argument(
        "--cuda-version",
        help="Override conda-forge::cuda-version pin in environment.gpu.yml",
    )
    parser.add_argument(
        "--await",
        dest="await_timeout",
        default="25m",
        help="Wave --await timeout (default: 25m)",
    )
    parser.add_argument(
        "--skip-gpu",
        action="store_true",
        help="Only rebuild CPU images; leave environment.gpu.yml and GPU URLs alone",
    )
    parser.add_argument(
        "--skip-anaconda-check",
        action="store_true",
        help="Do not verify the version exists on anaconda.org before freezing",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Rewrite nothing and do not call Wave; print the planned work",
    )
    parser.add_argument(
        "--cpu-docker",
        help="Skip Wave CPU Docker freeze and use this image reference",
    )
    parser.add_argument(
        "--cpu-sif",
        help="Skip Wave CPU Singularity freeze and use this HTTPS SIF URL",
    )
    parser.add_argument(
        "--gpu-docker",
        help="Skip Wave GPU Docker freeze and use this image reference",
    )
    parser.add_argument(
        "--gpu-sif",
        help="Skip Wave GPU Singularity freeze and use this HTTPS SIF URL",
    )
    return parser.parse_args()


def die(message: str, code: int = 1) -> None:
    print(f"error: {message}", file=sys.stderr)
    raise SystemExit(code)


def info(message: str) -> None:
    print(message, flush=True)


def require_version(version: str) -> str:
    version = version.strip().lstrip("v")
    if not VERSION_RE.match(version):
        die(f"invalid ginfinity version {version!r}")
    return version


def current_pin(text: str, pattern: re.Pattern[str], label: str) -> str:
    match = pattern.search(text)
    if not match:
        die(f"could not find {label} in environment file")
    return match.group(0).split("=", 1)[1].strip()


def rewrite_env(
    text: str,
    version: str,
    pytorch_gpu: str | None = None,
    cuda_version: str | None = None,
) -> str:
    if not GIN_DEP_RE.search(text):
        die("environment file has no nicolas.aira::ginfinity= pin")
    text = GIN_DEP_RE.sub(rf"\g<1>{version}", text)
    if pytorch_gpu is not None:
        if not PYTORCH_DEP_RE.search(text):
            die("environment file has no conda-forge::pytorch-gpu= pin")
        text = PYTORCH_DEP_RE.sub(rf"\g<1>{pytorch_gpu}", text)
    if cuda_version is not None:
        if not CUDA_DEP_RE.search(text):
            die("environment file has no conda-forge::cuda-version= pin")
        text = CUDA_DEP_RE.sub(rf"\g<1>{cuda_version}", text)
    return text


def conda_channels(text: str) -> str:
    channels: list[str] = []
    in_channels = False
    for raw in text.splitlines():
        line = raw.rstrip()
        if line == "channels:":
            in_channels = True
            continue
        if in_channels:
            if line.startswith("  - "):
                channels.append(line[4:].strip())
                continue
            break
    if not channels:
        die("environment file lists no conda channels")
    return ",".join(channels)


def parse_wave_json(stdout: str) -> dict:
    stdout = stdout.strip()
    if not stdout:
        die("Wave produced no stdout")
    try:
        return json.loads(stdout)
    except json.JSONDecodeError:
        start = stdout.rfind("{")
        if start < 0:
            die(f"Wave stdout was not JSON:\n{stdout}")
        return json.loads(stdout[start:])


def run_wave(args: list[str], *, dry_run: bool) -> dict | None:
    cmd = ["wave", *args]
    info("+ " + subprocess.list2cmdline(cmd))
    if dry_run:
        return None
    proc = subprocess.run(cmd, check=False, text=True, capture_output=True)
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip()
        die(f"Wave failed ({proc.returncode})\n{tail}")
    payload = parse_wave_json(proc.stdout)
    if payload.get("succeeded") is False:
        die(f"Wave reported failure: {json.dumps(payload, indent=2)}")
    return payload


def wave_freeze(env_file: Path, *, singularity: bool, await_timeout: str, dry_run: bool) -> WaveImage | None:
    channels = conda_channels(env_file.read_text())
    args = [
        "--conda-file",
        str(env_file),
        "--conda-channels",
        channels,
        "--freeze",
        "--await",
        await_timeout,
        "-o",
        "json",
    ]
    if singularity:
        args.append("--singularity")
    payload = run_wave(args, dry_run=dry_run)
    if payload is None:
        return None
    image = payload.get("containerImage") or payload.get("targetImage")
    build_id = payload.get("buildId") or ""
    if not image:
        die(f"Wave JSON missing containerImage:\n{json.dumps(payload, indent=2)}")
    return WaveImage(
        container_image=str(image),
        build_id=str(build_id),
        cached=bool(payload.get("cached")),
    )


def inspect_sif_digest(oras_image: str) -> str:
    payload = run_wave(["--inspect", "-i", oras_image, "-o", "json"], dry_run=False)
    assert payload is not None
    container = payload.get("container", payload)
    layers = container.get("manifest", {}).get("layers", [])
    if not layers:
        die(f"Wave inspect returned no layers for {oras_image}")
    for layer in layers:
        media = str(layer.get("mediaType", "")).lower()
        if "sif" in media:
            return str(layer["digest"])
    return str(max(layers, key=lambda layer: int(layer.get("size") or 0))["digest"])


def sif_https_url(digest: str) -> str:
    digest = digest.removeprefix("sha256:")
    if len(digest) != 64 or any(ch not in "0123456789abcdef" for ch in digest):
        die(f"unexpected SIF digest {digest!r}")
    return (
        "https://community-cr-prod.seqera.io/docker/registry/v2/blobs/"
        f"sha256/{digest[:2]}/{digest}/data"
    )


def wait_for_sif(url: str, timeout_s: int = 180) -> int:
    deadline = time.time() + timeout_s
    last_error = "no response"
    while time.time() < deadline:
        request = urllib.request.Request(url, method="HEAD")
        try:
            with urllib.request.urlopen(request, timeout=30) as response:
                length = int(response.headers.get("Content-Length") or 0)
                if response.status == 200 and length > 1_000_000:
                    return length
                last_error = f"HTTP {response.status} length={length}"
        except (urllib.error.URLError, TimeoutError, ValueError) as exc:
            last_error = str(exc)
        time.sleep(8)
    die(f"SIF URL not ready after {timeout_s}s ({last_error}): {url}")
    raise AssertionError("unreachable")


def check_anaconda(version: str) -> None:
    url = f"https://api.anaconda.org/package/{CHANNEL}/{PACKAGE}"
    info(f"Checking {url} for version {version}")
    request = urllib.request.Request(url, headers={"Accept": "application/json"})
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            payload = json.loads(response.read().decode())
    except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
        die(f"could not query anaconda.org: {exc}")
    versions = {str(item) for item in payload.get("versions", [])}
    if version not in versions:
        files = payload.get("files") or []
        file_versions = {str(item.get("version")) for item in files if isinstance(item, dict)}
        if version not in file_versions:
            available = ", ".join(sorted(versions)[-8:]) or "none"
            die(f"{CHANNEL}::{PACKAGE}={version} is not on anaconda.org (latest: {available})")
    info(f"Found {CHANNEL}::{PACKAGE}={version} on anaconda.org")


def replace_build_containers(text: str, cpu_sif: str, cpu_docker: str) -> str:
    if not BUILD_CONTAINER_RE.search(text):
        die(f"could not parse container directive in {BUILD_MAIN.relative_to(REPO_ROOT)}")
    return BUILD_CONTAINER_RE.sub(
        rf"\g<1>{cpu_sif}\g<3>{cpu_docker}\g<5>",
        text,
        count=1,
    )


def replace_embed_containers(
    text: str,
    gpu_sif: str,
    cpu_sif: str,
    gpu_docker: str,
    cpu_docker: str,
) -> str:
    if not EMBED_CONTAINER_RE.search(text):
        die(f"could not parse container directive in {EMBED_MAIN.relative_to(REPO_ROOT)}")
    return EMBED_CONTAINER_RE.sub(
        rf"\g<1>{gpu_sif}\g<3>{cpu_sif}\g<5>{gpu_docker}\g<7>{cpu_docker}\g<9>",
        text,
        count=1,
    )


def replace_meta(text: str, pins: Pins) -> str:
    replacements = [
        (
            r"(containers:\n  docker:\n    linux/amd64:\n      name: )\S+",
            rf"\g<1>{pins.cpu_docker}",
            "docker name",
        ),
        (
            r"(containers:\n  docker:\n    linux/amd64:\n      name: .+\n      build_id: )\S+",
            rf"\g<1>{pins.cpu_docker_build_id}",
            "docker build_id",
        ),
        (
            r"(singularity:\n    linux/amd64:\n      name: )\S+",
            rf"\g<1>{pins.cpu_oras}",
            "singularity name",
        ),
        (
            r"(singularity:\n    linux/amd64:\n      name: .+\n      build_id: )\S+",
            rf"\g<1>{pins.cpu_sif_build_id}",
            "singularity build_id",
        ),
        (
            r"(singularity:\n    linux/amd64:\n      name: .+\n      build_id: .+\n      https: )\S+",
            rf"\g<1>{pins.cpu_sif}",
            "singularity https",
        ),
    ]
    for pattern, repl, label in replacements:
        text, count = re.subn(pattern, repl, text, count=1)
        if count != 1:
            die(f"could not update {BUILD_META.relative_to(REPO_ROOT)} {label}")
    return text


def replace_agents(
    text: str,
    version: str,
    pins: Pins,
    pytorch_gpu: str,
    cuda_version: str,
) -> str:
    text, n = re.subn(
        r"nicolas\.aira::ginfinity==\S+",
        f"nicolas.aira::ginfinity=={version}",
        text,
        count=1,
    )
    if n != 1:
        info("warning: AGENTS.md ginfinity version pin not found")
    text, n = re.subn(
        r"CPU Docker: `[^`]+`",
        f"CPU Docker: `{pins.cpu_docker}`",
        text,
        count=1,
    )
    if n != 1:
        info("warning: AGENTS.md CPU Docker pin not found")
    if pins.gpu_docker:
        text, n = re.subn(
            r"GPU Docker: `[^`]+`(?: \(`pytorch-gpu=[^`]+`, `cuda-version=[^`]+`\))?",
            f"GPU Docker: `{pins.gpu_docker}` (`pytorch-gpu={pytorch_gpu}`, `cuda-version={cuda_version}`)",
            text,
            count=1,
        )
        if n != 1:
            info("warning: AGENTS.md GPU Docker pin not found")
    if pins.gpu_sif:
        digest = pins.gpu_sif.rstrip("/").split("/")[-2]
        text, n = re.subn(
            r"GPU Singularity: community-cr-prod SIF `[0-9a-f]+`",
            f"GPU Singularity: community-cr-prod SIF `{digest}`",
            text,
            count=1,
        )
        if n != 1:
            info("warning: AGENTS.md GPU Singularity pin not found")
    return text


def write_text(path: Path, text: str) -> None:
    path.write_text(text)
    info(f"Updated {path.relative_to(REPO_ROOT)}")


def freeze_sif(env_file: Path, await_timeout: str, dry_run: bool) -> tuple[WaveImage | None, str | None]:
    image = wave_freeze(env_file, singularity=True, await_timeout=await_timeout, dry_run=dry_run)
    if image is None:
        return None, None
    digest = inspect_sif_digest(image.container_image)
    url = sif_https_url(digest)
    size = wait_for_sif(url)
    info(f"SIF ready ({size} bytes): {url}")
    return image, url


def provided_urls(args: argparse.Namespace) -> bool:
    supplied = [args.cpu_docker, args.cpu_sif, args.gpu_docker, args.gpu_sif]
    if args.skip_gpu:
        supplied = supplied[:2]
    present = [item is not None for item in supplied]
    if any(present) and not all(present):
        die("provide all container URL overrides for the selected path, or none")
    return all(present) and any(present)


def main() -> int:
    args = parse_args()
    version = require_version(args.version)

    if shutil.which("wave") is None and not args.dry_run and not provided_urls(args):
        die("wave CLI not found on PATH")

    for path in (BUILD_ENV, BUILD_MAIN, EMBED_ENV, EMBED_GPU_ENV, EMBED_MAIN):
        if not path.is_file():
            die(f"missing {path}")

    if not args.skip_anaconda_check and not args.dry_run:
        check_anaconda(version)

    build_env = BUILD_ENV.read_text()
    embed_env = EMBED_ENV.read_text()
    gpu_env = EMBED_GPU_ENV.read_text()
    pytorch_gpu = args.pytorch_gpu or current_pin(gpu_env, PYTORCH_DEP_RE, "pytorch-gpu")
    cuda_version = args.cuda_version or current_pin(gpu_env, CUDA_DEP_RE, "cuda-version")

    new_build_env = rewrite_env(build_env, version)
    new_embed_env = rewrite_env(embed_env, version)
    new_gpu_env = gpu_env
    if not args.skip_gpu:
        new_gpu_env = rewrite_env(
            gpu_env,
            version,
            pytorch_gpu=pytorch_gpu,
            cuda_version=cuda_version,
        )

    info(f"GINFINITY {current_pin(build_env, GIN_DEP_RE, 'ginfinity')} -> {version}")
    if not args.skip_gpu:
        info(f"GPU pins: pytorch-gpu={pytorch_gpu} cuda-version={cuda_version}")

    if args.dry_run:
        info(f"Would update {BUILD_ENV.relative_to(REPO_ROOT)}")
        info(f"Would update {EMBED_ENV.relative_to(REPO_ROOT)}")
        if not args.skip_gpu:
            info(f"Would update {EMBED_GPU_ENV.relative_to(REPO_ROOT)}")
        if not BUILD_CONTAINER_RE.search(BUILD_MAIN.read_text()):
            die(f"dry-run: cannot parse {BUILD_MAIN.relative_to(REPO_ROOT)}")
        if not EMBED_CONTAINER_RE.search(EMBED_MAIN.read_text()):
            die(f"dry-run: cannot parse {EMBED_MAIN.relative_to(REPO_ROOT)}")
        info("Container directives parse OK")
        wave_freeze(BUILD_ENV, singularity=False, await_timeout=args.await_timeout, dry_run=True)
        wave_freeze(BUILD_ENV, singularity=True, await_timeout=args.await_timeout, dry_run=True)
        if not args.skip_gpu:
            wave_freeze(EMBED_GPU_ENV, singularity=False, await_timeout=args.await_timeout, dry_run=True)
            wave_freeze(EMBED_GPU_ENV, singularity=True, await_timeout=args.await_timeout, dry_run=True)
        info("Dry run only; no files written")
        return 0

    write_text(BUILD_ENV, new_build_env)
    write_text(EMBED_ENV, new_embed_env)
    if not args.skip_gpu:
        write_text(EMBED_GPU_ENV, new_gpu_env)

    if provided_urls(args):
        pins = Pins(
            cpu_docker=args.cpu_docker,
            cpu_sif=args.cpu_sif,
            cpu_docker_build_id="",
            cpu_sif_build_id="",
            cpu_oras="",
            gpu_docker=None if args.skip_gpu else args.gpu_docker,
            gpu_sif=None if args.skip_gpu else args.gpu_sif,
        )
    else:
        info("Freezing CPU Docker image")
        cpu_docker_img = wave_freeze(
            BUILD_ENV, singularity=False, await_timeout=args.await_timeout, dry_run=False
        )
        assert cpu_docker_img is not None
        info(f"CPU Docker: {cpu_docker_img.docker_ref} (cached={cpu_docker_img.cached})")

        info("Freezing CPU Singularity image")
        cpu_sif_img, cpu_sif = freeze_sif(BUILD_ENV, args.await_timeout, dry_run=False)
        assert cpu_sif_img is not None and cpu_sif is not None

        gpu_docker_img = None
        gpu_sif_img = None
        gpu_sif = None
        if not args.skip_gpu:
            info("Freezing GPU Docker image")
            gpu_docker_img = wave_freeze(
                EMBED_GPU_ENV, singularity=False, await_timeout=args.await_timeout, dry_run=False
            )
            assert gpu_docker_img is not None
            info(f"GPU Docker: {gpu_docker_img.docker_ref} (cached={gpu_docker_img.cached})")

            info("Freezing GPU Singularity image")
            gpu_sif_img, gpu_sif = freeze_sif(EMBED_GPU_ENV, args.await_timeout, dry_run=False)
            assert gpu_sif_img is not None and gpu_sif is not None

        pins = Pins(
            cpu_docker=cpu_docker_img.docker_ref,
            cpu_sif=cpu_sif,
            cpu_docker_build_id=cpu_docker_img.build_id,
            cpu_sif_build_id=cpu_sif_img.build_id,
            cpu_oras=cpu_sif_img.container_image,
            gpu_docker=None if gpu_docker_img is None else gpu_docker_img.docker_ref,
            gpu_sif=gpu_sif,
            gpu_docker_build_id=None if gpu_docker_img is None else gpu_docker_img.build_id,
            gpu_sif_build_id=None if gpu_sif_img is None else gpu_sif_img.build_id,
            gpu_oras=None if gpu_sif_img is None else gpu_sif_img.container_image,
        )

    write_text(BUILD_MAIN, replace_build_containers(BUILD_MAIN.read_text(), pins.cpu_sif, pins.cpu_docker))
    embed_text = EMBED_MAIN.read_text()
    if args.skip_gpu:
        match = EMBED_CONTAINER_RE.search(embed_text)
        if not match:
            die(f"could not parse container directive in {EMBED_MAIN.relative_to(REPO_ROOT)}")
        gpu_sif = match.group(2)
        gpu_docker = match.group(6)
    else:
        assert pins.gpu_sif and pins.gpu_docker
        gpu_sif = pins.gpu_sif
        gpu_docker = pins.gpu_docker
    write_text(
        EMBED_MAIN,
        replace_embed_containers(embed_text, gpu_sif, pins.cpu_sif, gpu_docker, pins.cpu_docker),
    )

    if BUILD_META.is_file() and pins.cpu_docker_build_id:
        write_text(BUILD_META, replace_meta(BUILD_META.read_text(), pins))

    if AGENTS.is_file() and not provided_urls(args):
        write_text(
            AGENTS,
            replace_agents(AGENTS.read_text(), version, pins, pytorch_gpu, cuda_version),
        )

    info("Done")
    info(f"  CPU Docker      {pins.cpu_docker}")
    info(f"  CPU Singularity {pins.cpu_sif}")
    if pins.gpu_docker:
        info(f"  GPU Docker      {pins.gpu_docker}")
    if pins.gpu_sif:
        info(f"  GPU Singularity {pins.gpu_sif}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        die("interrupted", 130)
