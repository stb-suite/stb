#!/usr/bin/env python3
"""Checks that the suite has ONE version everywhere (stb/_version.py).

- `stb.__version__` is `RELEASE` on a release branch (DEV = False), or
  `RELEASE.dev<N>` on the development branch, N = commits since BASE_TAG;
  `STB_VERSION` overrides it.
- Every console command reports exactly that version in its `--version`, and the
  interactive menu's own VERSION is the same object. No tool defines a version
  of its own any more.
- The banners all carry the same year.
- Building the package reads the version from the git checkout (and honours
  STB_VERSION), and an installed copy, which has no checkout, reports the number
  its metadata was built with. This is checked on a throw-away git repository
  copied from `stb-suite/`, in both modes (release, and dev with a tag and one
  commit after it), so it neither touches nor depends on the real
  `stb-suite/build/` or the history here.

Usage: check_versions.py <work-dir>     (exits 1 on any failure)
Needs: the stb-suite package and its `ml` extra installed (the ML tools refuse
to start without MACE); git; setuptools.
"""
import os
import re
import shutil
import subprocess
import sys
import tomllib
import zipfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "stb-suite" / "src" / "stb"
failures = []


def fail(message):
    failures.append(message)
    print(f"[FAIL] {message}")


def ok(message):
    print(f"[ OK ] {message}")


def run(args, **kwargs):
    return subprocess.run(args, capture_output=True, text=True, stdin=subprocess.DEVNULL, **kwargs)


def _version_module():
    import stb._version as v
    return v


def check_runtime_version():
    import stb
    v = _version_module()
    version = stb.__version__
    if not v.DEV:
        (ok if version == v.RELEASE else fail)(f"stb.__version__ = {version} (release branch, RELEASE = {v.RELEASE})")
    else:
        count = run(["git", "rev-list", "--count", f"{v.BASE_TAG}..HEAD"], cwd=ROOT).stdout.strip()
        expected = f"{v.RELEASE}.dev{count}"
        (ok if version == expected else fail)(
            f"stb.__version__ = {version} (dev branch, expected {expected}: {count} commits since {v.BASE_TAG})")
    env = {**os.environ, "STB_VERSION": "9.9.4242"}
    got = run([sys.executable, "-c", "import stb; print(stb.__version__)"], env=env).stdout.strip()
    (ok if got == "9.9.4242" else fail)(f"STB_VERSION overrides the version (got {got})")
    return version


def check_commands(version):
    scripts = tomllib.loads((ROOT / "stb-suite" / "pyproject.toml").read_text())["project"]["scripts"]
    menu = "stb-suite"

    def version_of(item):
        name, target = item
        module, function = target.split(":")
        code = (f"import sys; sys.argv = [{name!r}, '--version']; "
                f"from {module} import {function} as main; main()")
        env = {**os.environ, "PYTHONWARNINGS": "ignore", "NO_COLOR": "1"}
        result = run([sys.executable, "-c", code], env=env, timeout=300)
        lines = result.stdout.strip().splitlines()
        return name, result.returncode, (lines[-1] if lines else ""), result.stderr.strip()[-200:]

    with ThreadPoolExecutor(max_workers=8) as pool:
        results = list(pool.map(version_of, [(n, t) for n, t in scripts.items() if n != menu]))
    bad = [(n, rc, line, err) for n, rc, line, err in results
           if rc != 0 or not line.startswith("stb") or not line.endswith(f" {version}")]
    for name, rc, line, err in bad:
        fail(f"{name} --version -> rc={rc}, {line!r} {err}")
    if not bad:
        ok(f"all {len(results)} commands print '... {version}' for --version")
    import stb.stb_suite
    (ok if stb.stb_suite.VERSION == version else fail)(f"{menu}'s VERSION is {stb.stb_suite.VERSION}")


def check_sources():
    literal = re.compile(r'^VERSION\s*=\s*["\']', re.M)
    own = [p.name for p in sorted(SRC.glob("*.py")) if literal.search(p.read_text(encoding="utf-8"))]
    (fail if own else ok)(f"no tool defines its own VERSION literal" + (f": {own}" if own else ""))
    years = set()
    for p in SRC.glob("*.py"):
        years |= set(re.findall(r"University of Brasilia - (\d{4})", p.read_text(encoding="utf-8")))
    (ok if len(years) == 1 else fail)(f"the banners carry one year ({sorted(years)})")


def check_build(work):
    work = Path(work)
    v = _version_module()
    git = ["git", "-c", "user.name=t", "-c", "user.email=t@t"]

    def make_repo(name, dev):
        repo = work / name
        if repo.exists():
            shutil.rmtree(repo)
        shutil.copytree(ROOT / "stb-suite", repo / "stb-suite",
                        ignore=shutil.ignore_patterns("build", "dist", "__pycache__", "*.egg-info"))
        vf = repo / "stb-suite" / "src" / "stb" / "_version.py"
        text = vf.read_text(encoding="utf-8")
        text = re.sub(r"^DEV = (True|False)$", f"DEV = {dev}", text, flags=re.M)
        vf.write_text(text, encoding="utf-8")
        run(["git", "init", "-q"], cwd=repo)
        run(["git", "add", "-A"], cwd=repo)
        run(git + ["commit", "-q", "-m", "one"], cwd=repo)
        run(git + ["tag", v.BASE_TAG], cwd=repo)
        run(git + ["commit", "-q", "--allow-empty", "-m", "two"], cwd=repo)
        return repo

    def wheel(repo, env_extra):
        out = work / "wheels"
        shutil.rmtree(out, ignore_errors=True)
        env = {k: val for k, val in os.environ.items() if k != "STB_VERSION"} | env_extra
        result = run([sys.executable, "-m", "pip", "wheel", "--no-deps", "--no-build-isolation",
                      "-q", "-w", str(out), str(repo / "stb-suite")], env=env, cwd=work)
        found = sorted(out.glob("stb_suite-*.whl"))
        return (found[0] if found else None), result.stderr[-300:]

    release_repo = make_repo("repo_release", False)
    built, err = wheel(release_repo, {})
    if built is None:
        return fail(f"building the package failed: {err}")
    version = built.name.split("-")[1]
    (ok if version == v.RELEASE else fail)(f"a release-branch build is version {version} (expected {v.RELEASE})")
    installed = work / "installed"
    shutil.rmtree(installed, ignore_errors=True)
    zipfile.ZipFile(built).extractall(installed)
    got = run([sys.executable, "-c", "import stb; print(stb.__version__, stb.__file__)"],
              env={**os.environ, "PYTHONPATH": str(installed)}, cwd=work).stdout.split()
    (ok if got and got[0] == version and str(installed) in got[-1] else fail)(
        f"an installed copy (no checkout) reports its metadata version ({got})")
    override, err = wheel(release_repo, {"STB_VERSION": "9.9.999"})
    (ok if override is not None and override.name.split("-")[1] == "9.9.999" else fail)(
        f"STB_VERSION sets the built version ({override.name if override else err})")
    dev_built, err = wheel(make_repo("repo_dev", True), {})
    expected = f"{v.RELEASE}.dev1"
    got_dev = dev_built.name.split("-")[1] if dev_built else err
    (ok if got_dev == expected else fail)(
        f"a dev-branch build one commit after {v.BASE_TAG} is version {got_dev} (expected {expected})")


def main():
    work = sys.argv[1]
    os.makedirs(work, exist_ok=True)
    version = check_runtime_version()
    if version:
        check_commands(version)
    check_sources()
    check_build(work)
    print(f"{len(failures)} failure(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
