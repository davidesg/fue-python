"""
fue-bug — command-line front end for the in-repo bug tracker (fue.bugs).

    fue-bug list [--status S] [--component C]   list reports (open by default: all)
    fue-bug show BUG-NNNN                        print a report
    fue-bug new "title" --component forecast     create a new report
    fue-bug index                               regenerate bugs/README.md
    fue-bug check                               validate all reports (CI-friendly)
"""

from __future__ import annotations

import argparse
import sys

from . import bugs as _bugs


def _cmd_list(args):
    items = _bugs.list_bugs(status=args.status, component=args.component)
    if not items:
        print("no bug reports found.")
        return 0
    for b in items:
        flag = " " if b.status != "open" else "*"
        print(f"{flag}{b.id}  [{b.status:11s}] {b.severity:8s} "
              f"{b.component:10s} {b.title}")
    n_open = sum(1 for b in items if b.is_open)
    print(f"\n{len(items)} report(s), {n_open} open.")
    return 0


def _cmd_show(args):
    for b in _bugs.list_bugs():
        if b.id == args.id:
            print(_bugs.render_frontmatter(b))
            print()
            print(b.body)
            return 0
    print(f"fue-bug: no report with id {args.id}", file=sys.stderr)
    return 1


def _cmd_new(args):
    path = _bugs.new_bug(
        args.title, component=args.component, severity=args.severity,
        found_in=args.found_in, reporter=args.reporter,
        tags=args.tag or [])
    print(f"created {path}")
    print("edit it, then run 'fue-bug index' to refresh bugs/README.md")
    return 0


def _cmd_index(args):
    path = _bugs.write_index()
    print(f"wrote {path}")
    return 0


def _cmd_check(args):
    report = _bugs.validate_all()
    if not report:
        n = len(_bugs.list_bugs())
        print(f"OK — {n} report(s), all valid.")
        return 0
    print("INVALID bug reports:", file=sys.stderr)
    for who, errs in report.items():
        for e in errs:
            print(f"  {who}: {e}", file=sys.stderr)
    return 1


def main(argv=None):
    p = argparse.ArgumentParser(prog="fue-bug",
                                description="fue in-repo bug tracker")
    sub = p.add_subparsers(dest="cmd")

    pl = sub.add_parser("list", help="list bug reports")
    pl.add_argument("--status", choices=_bugs.STATUSES)
    pl.add_argument("--component")
    pl.set_defaults(func=_cmd_list)

    ps = sub.add_parser("show", help="print a report")
    ps.add_argument("id")
    ps.set_defaults(func=_cmd_show)

    pn = sub.add_parser("new", help="create a new report")
    pn.add_argument("title")
    pn.add_argument("--component", required=True)
    pn.add_argument("--severity", choices=_bugs.SEVERITIES, default="medium")
    pn.add_argument("--found-in", dest="found_in", default="")
    pn.add_argument("--reporter", default="")
    pn.add_argument("--tag", action="append")
    pn.set_defaults(func=_cmd_new)

    pi = sub.add_parser("index", help="regenerate bugs/README.md")
    pi.set_defaults(func=_cmd_index)

    pc = sub.add_parser("check", help="validate all reports")
    pc.set_defaults(func=_cmd_check)

    args = p.parse_args(argv)
    if not getattr(args, "func", None):
        p.print_help()
        return 1
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
