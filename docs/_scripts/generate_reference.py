#!/usr/bin/env python3
"""Generate the command reference pages in docs/reference/ from a SeekDeep binary.

Runs `SeekDeep` to list the commands, runs `SeekDeep <command> --help` for each one, and writes
one page per command group (SeekDeep, Pmo, SeekDeepServer, SeekDeepUtils) plus an index, so the
reference always matches what the program accepts. Hidden aliases (e.g. clusterDown for qluster)
are read from the runners' addFunc(...) calls in src/.

Run it from anywhere after building SeekDeep, and commit the pages it writes:

    python3 docs/_scripts/generate_reference.py                 # uses ./bin/SeekDeep
    python3 docs/_scripts/generate_reference.py --seekdeep /path/to/SeekDeep
"""

import argparse
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

ANSI = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
COMMAND_LINE = re.compile(r"^\s*\d+\)\s+(\S+)\s*$")
OPTION_TAIL = re.compile(
    r"^(?P<desc>.*?)(?:; (?P<required>required))?; default=(?P<default>.*); \((?P<type>[^()]*)\)$"
)
ADD_FUNC = re.compile(r'addFunc\(\s*"(?P<name>[^"]+)"\s*,\s*(?P<func>[A-Za-z0-9_:]+)\s*,\s*(?P<hidden>true|false)\s*\)')
TOP_LEVEL_HEADERS = {"Required Options", "Optional Options"}
# every command takes these, so they're mentioned once at the top of each page instead of in every table
UNIVERSAL_FLAGS = {"--getFlags", "-h", "--help"}

REPO = Path(__file__).resolve().parents[2]
OUT_DIR = REPO / "docs" / "reference"

GROUP_BLURBS = {
    "SeekDeep": "The main pipeline programs: extraction, per-sample clustering and population clustering.",
    "Pmo": "Building, merging and validating [Portable Microhaplotype Object (PMO)](https://plasmogenepi.github.io/PMO_Docs/) files.",
    "SeekDeepServer": "Local web servers for viewing results.",
    "SeekDeepUtils": "Setup, QC, benchmarking, variant calling and other utilities.",
}


def run(cmd):
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, errors="replace")
    return ANSI.sub("", proc.stdout)


def list_commands(seekdeep):
    """Return (version, {group: [command, ...]}) from the program's own listing."""
    out = run([seekdeep])
    version = ""
    groups = {}
    current = None
    for line in out.splitlines():
        if line.startswith("Version "):
            version = line.split(" ", 1)[1].strip()
            continue
        m = COMMAND_LINE.match(line)
        if m:
            if current is None:
                sys.exit(f"command listed before any group: {line!r}")
            groups[current].append(m.group(1))
            continue
        stripped = line.strip()
        if stripped and " " not in stripped:
            current = stripped
            groups[current] = []
    groups = {g: cmds for g, cmds in groups.items() if cmds}
    if not groups:
        sys.exit(f"could not find any commands in the output of {seekdeep}")
    return version, groups


def find_aliases():
    """Map a visible command name to its hidden aliases, from addFunc(name, func, hidden) calls."""
    by_func = defaultdict(lambda: {"visible": [], "hidden": []})
    for src in (REPO / "src").rglob("*Runner.cpp"):
        for m in ADD_FUNC.finditer(src.read_text(errors="replace")):
            key = (src, m.group("func"))
            by_func[key]["hidden" if m.group("hidden") == "true" else "visible"].append(m.group("name"))
    aliases = defaultdict(list)
    for names in by_func.values():
        for visible in names["visible"]:
            aliases[visible].extend(names["hidden"])
    return aliases


def parse_help(text):
    """Split --help output into a description and a list of (top, group, [options])."""
    lines = [l.rstrip() for l in text.splitlines()]
    while lines and not lines[0].strip():
        lines.pop(0)
    if lines and lines[0].startswith("SeekDeep "):
        lines.pop(0)
    description = []
    sections = []
    top = None
    group = None
    for line in lines:
        if "\t" in line:
            flags, rest = line.split("\t", 1)
            m = OPTION_TAIL.match(rest.strip())
            if m:
                opt = {
                    "flags": [f.strip() for f in flags.strip().split(",") if f.strip()],
                    "desc": m.group("desc").strip(),
                    "required": bool(m.group("required")),
                    "default": m.group("default").strip(),
                    "type": m.group("type").strip(),
                }
            else:
                opt = {"flags": [f.strip() for f in flags.strip().split(",") if f.strip()],
                       "desc": rest.strip(), "required": False, "default": "", "type": ""}
            if set(opt["flags"]) & UNIVERSAL_FLAGS:
                continue
            if not sections or sections[-1][0] != top or sections[-1][1] != group:
                sections.append((top, group, []))
            sections[-1][2].append(opt)
            continue
        stripped = line.strip()
        if not stripped:
            continue
        if stripped in TOP_LEVEL_HEADERS:
            top = stripped
            group = None
        elif top is None:
            description.append(stripped)
        else:
            group = stripped
    return " ".join(description), sections


# defaults that are paths on the machine that built SeekDeep (server resource directories) say
# nothing useful to a reader and leak the build machine's paths, so they're replaced
BUILD_PATH = re.compile(r"^(/Users/|/home/|/opt/|/usr/local/|/private/)")


def clean_default(text):
    if BUILD_PATH.match(text) and not text.startswith("/tmp"):
        return "*(set when SeekDeep is built)*"
    return code_cell(text)


def md_cell(text):
    return text.replace("|", "\\|").replace("<", "\\<").replace(">", "\\>")


def short_description(text, limit=160):
    """First sentence of a description, for the index table."""
    if not text:
        return "*(no description yet)*"
    first = re.split(r"(?<=[.;])\s", text, maxsplit=1)[0]
    if len(first) > limit:
        first = first[:limit].rsplit(" ", 1)[0] + "..."
    return md_cell(first)


def code_cell(text):
    if not text:
        return ""
    return "`" + text.replace("|", "\\|") + "`"


def write_group_page(group, commands, helps, aliases, version):
    slug = group.lower()
    out = [
        "---",
        f'title: "{group} commands"',
        "toc-depth: 2",
        "---",
        "",
        f"<!-- Generated by docs/_scripts/generate_reference.py from SeekDeep {version}. Do not edit by hand. -->",
        "",
        GROUP_BLURBS.get(group, ""),
        "",
        f"Generated from the `--help` output of SeekDeep {version}. Run `SeekDeep <command> --help` to see the same "
        "for your installed version. Commands and flags are not case sensitive, and every command also takes "
        "`--help`/`-h` and `--getFlags`.",
        "",
    ]
    for cmd in commands:
        description, sections = helps[cmd]
        out += [f"## {cmd} {{#{cmd}}}", ""]
        if description:
            out += [md_cell(description), ""]
        if aliases.get(cmd):
            out += ["Also callable as " + ", ".join(f"`{a}`" for a in aliases[cmd]) + ".", ""]
        out += ["```bash", f"SeekDeep {cmd} --help", "```", ""]
        if not sections:
            out += ["*This command lists no options.*", ""]
            continue
        last_top = None
        for top, sub, opts in sections:
            if top != last_top and top:
                out += [f"**{top}**", ""]
                last_top = top
            if sub:
                out += [f"*{sub}*", ""]
            out += ["| Flag | Description | Default | Type |", "|---|---|---|---|"]
            for o in opts:
                flags = "<br>".join(code_cell(f) for f in o["flags"])
                desc = md_cell(o["desc"]) + (" **(required)**" if o["required"] else "")
                out.append(f"| {flags} | {desc} | {clean_default(o['default'])} | {md_cell(o['type'])} |")
            out.append("")
    (OUT_DIR / f"{slug}.qmd").write_text("\n".join(out))
    return slug


def write_index(groups, helps, slugs, version):
    out = [
        "---",
        'title: "Command reference"',
        "---",
        "",
        f"<!-- Generated by docs/_scripts/generate_reference.py from SeekDeep {version}. Do not edit by hand. -->",
        "",
        f"Every command in SeekDeep {version}, with its flags as printed by `--help`. Run `SeekDeep` with no "
        "arguments to list the commands of your installed version.",
        "",
    ]
    for group, commands in groups.items():
        out += [f"## [{group}]({slugs[group]}.qmd)", "", GROUP_BLURBS.get(group, ""), "",
                "| Command | What it does |", "|---|---|"]
        for cmd in commands:
            out.append(f"| [{cmd}]({slugs[group]}.qmd#{cmd}) | {short_description(helps[cmd][0])} |")
        out.append("")
    (OUT_DIR / "index.qmd").write_text("\n".join(out))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--seekdeep", default=str(REPO / "bin" / "SeekDeep"), help="SeekDeep binary to document")
    args = ap.parse_args()

    version, groups = list_commands(args.seekdeep)
    aliases = find_aliases()
    helps = {}
    for commands in groups.values():
        for cmd in commands:
            helps[cmd] = parse_help(run([args.seekdeep, cmd, "--help"]))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    slugs = {g: write_group_page(g, cmds, helps, aliases, version) for g, cmds in groups.items()}
    write_index(groups, helps, slugs, version)
    n = sum(len(c) for c in groups.values())
    print(f"wrote {n} commands in {len(groups)} groups from SeekDeep {version} to {OUT_DIR}")


if __name__ == "__main__":
    main()
