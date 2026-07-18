"""
fue.bugs — lightweight in-repo bug tracker.

Every bug is one Markdown file under ``bugs/`` at the repository root, with a
small YAML-style frontmatter block and a free-form body.  This module parses
those files, validates them, lists/filters them, renders the index, and creates
new ones.  The workflow is:

    1.  file a report        →  ``fue-bug new`` (or copy TEMPLATE.md)
    2.  fix references the id →  commit ``fix(component): BUG-NNNN ...``
    3.  close it             →  set ``status: fixed`` and ``fixed_in: <version>``

No third-party dependencies: the frontmatter parser is intentionally tiny and
covers only the subset the schema uses (``key: value``, inline ``[a, b]`` lists,
and ``- item`` block lists).
"""

from __future__ import annotations

import os
import re
import datetime as _dt
from dataclasses import dataclass, field
from pathlib import Path

# ── Schema ────────────────────────────────────────────────────────────────────

STATUSES   = ("open", "in-progress", "fixed", "wontfix", "duplicate")
SEVERITIES = ("low", "medium", "high", "critical")

# Required / optional frontmatter keys.
_REQUIRED = ("id", "title", "status", "severity", "component", "reported")
_LIST_KEYS = ("tags", "references")

_ID_RE = re.compile(r"^BUG-\d{4}$")


@dataclass
class Bug:
    """A single bug report (frontmatter + markdown body)."""
    id: str
    title: str
    status: str
    severity: str
    component: str
    reported: str
    reporter: str = ""
    found_in: str = ""
    fixed_in: str = ""
    tags: list = field(default_factory=list)
    references: list = field(default_factory=list)
    body: str = ""
    path: Path | None = None

    # -- validation --------------------------------------------------------
    def problems(self):
        """Return a list of schema violations (empty ⇒ valid)."""
        errs = []
        if not _ID_RE.match(self.id or ""):
            errs.append(f"id '{self.id}' is not of the form BUG-NNNN")
        if not self.title:
            errs.append("title is empty")
        if self.status not in STATUSES:
            errs.append(f"status '{self.status}' not in {STATUSES}")
        if self.severity not in SEVERITIES:
            errs.append(f"severity '{self.severity}' not in {SEVERITIES}")
        if not self.component:
            errs.append("component is empty")
        if not _valid_date(self.reported):
            errs.append(f"reported '{self.reported}' is not YYYY-MM-DD")
        if self.status == "fixed" and not self.fixed_in:
            errs.append("status is 'fixed' but fixed_in is empty")
        return errs

    @property
    def is_open(self):
        return self.status in ("open", "in-progress")


# ── Frontmatter parsing / rendering ───────────────────────────────────────────

def _valid_date(s):
    try:
        _dt.date.fromisoformat(str(s))
        return True
    except (ValueError, TypeError):
        return False


def _split_frontmatter(text):
    """Return (frontmatter_str, body_str).  Frontmatter is between the first
    two ``---`` fences at the very top of the file."""
    if not text.startswith("---"):
        raise ValueError("file does not start with a '---' frontmatter fence")
    parts = text.split("\n")
    if parts[0].strip() != "---":
        raise ValueError("first line must be exactly '---'")
    end = None
    for i in range(1, len(parts)):
        if parts[i].strip() == "---":
            end = i
            break
    if end is None:
        raise ValueError("frontmatter closing '---' fence not found")
    return "\n".join(parts[1:end]), "\n".join(parts[end + 1:]).lstrip("\n")


def _parse_frontmatter(fm):
    """Tiny YAML-subset parser: scalars, inline ``[a, b]`` and ``- item`` lists."""
    data = {}
    key = None
    for raw in fm.split("\n"):
        line = raw.rstrip()
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        if line.startswith(("  - ", "- ", "\t- ")):        # block-list item
            if key is None:
                continue
            if not isinstance(data.get(key), list):         # promote scalar → list
                data[key] = []
            data[key].append(_strip(line.split("-", 1)[1]))
            continue
        if ":" not in line:
            continue
        k, v = line.split(":", 1)
        key = k.strip()
        v = v.strip()
        if v == "":
            data[key] = ""                # empty: scalar "" unless a "- " follows
        elif v.startswith("[") and v.endswith("]"):
            inner = v[1:-1].strip()
            data[key] = [_strip(x) for x in inner.split(",")] if inner else []
        else:
            data[key] = _strip(v)
    return data


def _strip(s):
    s = s.strip()
    if len(s) >= 2 and s[0] == s[-1] and s[0] in ("'", '"'):
        s = s[1:-1]
    return s


def _bug_from_text(text, path=None):
    fm, body = _split_frontmatter(text)
    d = _parse_frontmatter(fm)
    kwargs = dict(
        id=str(d.get("id", "")),
        title=str(d.get("title", "")),
        status=str(d.get("status", "")),
        severity=str(d.get("severity", "")),
        component=str(d.get("component", "")),
        reported=str(d.get("reported", "")),
        reporter=str(d.get("reporter", "")),
        found_in=str(d.get("found_in", "")),
        fixed_in=str(d.get("fixed_in", "")),
        tags=list(d.get("tags", []) or []),
        references=list(d.get("references", []) or []),
        body=body,
        path=Path(path) if path else None,
    )
    return Bug(**kwargs)


def render_frontmatter(bug):
    """Serialise a Bug's frontmatter back to text (stable key order)."""
    lines = ["---"]
    for k in ("id", "title", "status", "severity", "component",
              "found_in", "fixed_in", "reported", "reporter"):
        lines.append(f"{k}: {getattr(bug, k) or ''}")
    for k in _LIST_KEYS:
        vals = getattr(bug, k)
        if vals:
            lines.append(f"{k}:")
            lines.extend(f"  - {v}" for v in vals)
        else:
            lines.append(f"{k}: []")
    lines.append("---")
    return "\n".join(lines)


# ── Repository discovery / loading ────────────────────────────────────────────

def find_bugs_dir(start=None):
    """Locate the ``bugs/`` directory by walking up from *start* (cwd by
    default) and from this package's location.  Returns a Path or None."""
    seeds = [Path(start).resolve() if start else Path.cwd()]
    seeds.append(Path(__file__).resolve().parent)   # inside the installed pkg
    for seed in seeds:
        for d in (seed, *seed.parents):
            cand = d / "bugs"
            if cand.is_dir() and any(cand.glob("BUG-*.md")):
                return cand
    return None


def load_bug(path):
    return _bug_from_text(Path(path).read_text(encoding="utf-8"), path=path)


def list_bugs(bugs_dir=None, status=None, component=None):
    """Return all bugs (sorted by id), optionally filtered by status/component."""
    bugs_dir = Path(bugs_dir) if bugs_dir else find_bugs_dir()
    if bugs_dir is None:
        return []
    out = []
    for p in sorted(Path(bugs_dir).glob("BUG-*.md")):
        try:
            b = load_bug(p)
        except ValueError:
            continue
        if status and b.status != status:
            continue
        if component and b.component != component:
            continue
        out.append(b)
    return sorted(out, key=lambda b: b.id)


def validate_all(bugs_dir=None):
    """Return {bug_id_or_path: [problems]} for every malformed report, plus
    duplicate-id detection.  Empty dict ⇒ all good."""
    bugs_dir = Path(bugs_dir) if bugs_dir else find_bugs_dir()
    report = {}
    seen = {}
    if bugs_dir is None:
        return {"<bugs dir>": ["no bugs/ directory with BUG-*.md found"]}
    for p in sorted(Path(bugs_dir).glob("BUG-*.md")):
        try:
            b = load_bug(p)
        except ValueError as exc:
            report[p.name] = [f"parse error: {exc}"]
            continue
        errs = b.problems()
        if b.id in seen:
            errs.append(f"duplicate id (also in {seen[b.id]})")
        else:
            seen[b.id] = p.name
        if errs:
            report[b.id or p.name] = errs
    return report


def next_id(bugs_dir=None):
    bugs = list_bugs(bugs_dir)
    n = 0
    for b in bugs:
        m = re.match(r"BUG-(\d{4})$", b.id)
        if m:
            n = max(n, int(m.group(1)))
    return f"BUG-{n + 1:04d}"


# ── Index rendering ───────────────────────────────────────────────────────────

def render_index(bugs_dir=None):
    """Render bugs/README.md content: a status-grouped table of all reports."""
    bugs = list_bugs(bugs_dir)
    lines = [
        "# Bug reports",
        "",
        "In-repo bug tracker for **fue**.  One Markdown file per bug "
        "(`BUG-NNNN-slug.md`); this index is generated by `fue-bug index`.",
        "",
        "New report: `fue-bug new` (or copy `TEMPLATE.md`).  "
        "Validate: `fue-bug check`.  A fix commit references the id, e.g. "
        "`fix(forecast): BUG-0001 …`.",
        "",
    ]
    order = {s: i for i, s in enumerate(("open", "in-progress", "fixed",
                                         "wontfix", "duplicate"))}
    n_open = sum(1 for b in bugs if b.is_open)
    lines.append(f"**{len(bugs)} report(s), {n_open} open.**")
    lines.append("")
    lines.append("| id | status | sev | component | title | fixed in |")
    lines.append("|----|--------|-----|-----------|-------|----------|")
    for b in sorted(bugs, key=lambda b: (order.get(b.status, 9), b.id)):
        link = b.path.name if b.path else b.id
        lines.append(
            f"| [{b.id}]({link}) | {b.status} | {b.severity} | "
            f"{b.component} | {b.title} | {b.fixed_in or '—'} |")
    lines.append("")
    return "\n".join(lines)


def write_index(bugs_dir=None):
    bugs_dir = Path(bugs_dir) if bugs_dir else find_bugs_dir()
    if bugs_dir is None:
        raise FileNotFoundError("no bugs/ directory found")
    (bugs_dir / "README.md").write_text(render_index(bugs_dir) + "\n",
                                        encoding="utf-8")
    return bugs_dir / "README.md"


# ── Creation ──────────────────────────────────────────────────────────────────

def _slug(title):
    s = re.sub(r"[^a-z0-9]+", "-", title.lower()).strip("-")
    return s[:48] or "untitled"


def new_bug(title, *, component, severity="medium", found_in="",
            reporter="", tags=None, bugs_dir=None):
    """Create a new BUG-NNNN-slug.md from the template.  Returns its Path."""
    bugs_dir = Path(bugs_dir) if bugs_dir else find_bugs_dir()
    if bugs_dir is None:
        raise FileNotFoundError("no bugs/ directory found (create bugs/ first)")
    bid = next_id(bugs_dir)
    bug = Bug(
        id=bid, title=title, status="open", severity=severity,
        component=component, reported=_dt.date.today().isoformat(),
        reporter=reporter, found_in=found_in, tags=list(tags or []),
    )
    body = _TEMPLATE_BODY
    text = render_frontmatter(bug) + "\n\n" + body + "\n"
    path = bugs_dir / f"{bid}-{_slug(title)}.md"
    path.write_text(text, encoding="utf-8")
    return path


_TEMPLATE_BODY = """\
## Summary

One paragraph: what is wrong and its user-visible effect.

## Impact

Who/what is affected, and how badly.

## Reproduction

Minimal steps or code to reproduce.

## Root cause

Where in the code, and why.

## Fix

The proposed or applied change.

## Validation

How the fix is checked (tests, references)."""
