"""
Tests for the in-repo bug tracker (fue.bugs) and its reports.

These run in CI so that every committed bug report stays schema-valid, ids stay
unique, and the generated index (bugs/README.md) stays in sync.
"""

import os
from pathlib import Path

import pytest

from fue import bugs

# Locate bugs/ relative to this test file (repo root = two levels up).
BUGS_DIR = Path(__file__).resolve().parents[1] / "bugs"


def test_bugs_dir_exists():
    assert BUGS_DIR.is_dir(), f"missing {BUGS_DIR}"
    assert list(BUGS_DIR.glob("BUG-*.md")), "no BUG-*.md reports"


def test_all_reports_valid():
    problems = bugs.validate_all(BUGS_DIR)
    assert not problems, f"invalid bug reports: {problems}"


def test_ids_unique_and_well_formed():
    items = bugs.list_bugs(BUGS_DIR)
    ids = [b.id for b in items]
    assert len(ids) == len(set(ids)), f"duplicate ids: {ids}"
    for b in items:
        assert bugs._ID_RE.match(b.id), b.id
        assert b.status in bugs.STATUSES
        assert b.severity in bugs.SEVERITIES


def test_bug_0001_present_and_about_drift():
    items = {b.id: b for b in bugs.list_bugs(BUGS_DIR)}
    assert "BUG-0001" in items, "the forecast-drift bug (BUG-0001) must exist"
    b = items["BUG-0001"]
    assert b.component == "forecast"
    assert "drift" in (b.title + " " + b.body).lower()


def test_index_is_in_sync():
    readme = BUGS_DIR / "README.md"
    assert readme.is_file(), "run 'fue-bug index' to create bugs/README.md"
    content = readme.read_text(encoding="utf-8")
    # every report must appear in the index
    for b in bugs.list_bugs(BUGS_DIR):
        assert b.id in content, f"{b.id} missing from bugs/README.md — run 'fue-bug index'"


def test_frontmatter_roundtrip():
    b = bugs.list_bugs(BUGS_DIR)[0]
    fm = bugs.render_frontmatter(b)
    reparsed = bugs._bug_from_text(fm + "\n\nbody")
    assert reparsed.id == b.id
    assert reparsed.title == b.title
    assert reparsed.tags == b.tags
    assert reparsed.references == b.references


def test_new_bug_lifecycle(tmp_path):
    # a scratch bugs/ dir with one report so find/next_id work in isolation
    d = tmp_path / "bugs"
    d.mkdir()
    (d / "BUG-0001-seed.md").write_text(
        (BUGS_DIR / "BUG-0001-forecast-mean-drift.md").read_text(encoding="utf-8"),
        encoding="utf-8")
    assert bugs.next_id(d) == "BUG-0002"
    p = bugs.new_bug("A brand new problem", component="report",
                     severity="low", bugs_dir=d)
    assert p.exists() and p.name.startswith("BUG-0002-")
    assert not bugs.validate_all(d)                 # the new one is valid
    assert bugs.next_id(d) == "BUG-0003"


def test_status_filter():
    open_bugs = bugs.list_bugs(BUGS_DIR, status="open")
    for b in open_bugs:
        assert b.status == "open"
