#!/usr/bin/env python3
"""Find, and optionally fix, commits that list an AI tool as an author or co-author.

X-PSI's policy on generative AI (docs/source/ai_policy.rst) asks contributors
to record AI assistance with an ``Assisted-by:`` commit trailer rather than a
``Co-authored-by:`` trailer, because GitHub counts co-authors as contributors.

Usage:

    python3 .github/scripts/ai_coauthor.py check <base> <head>
        Report every commit in <base>..<head> that lists an AI tool as
        co-author or author. Exits with status 1 if any are found. Used by the
        "AI co-author check" GitHub Actions workflow.

    python3 .github/scripts/ai_coauthor.py fix
        Rewrite the current commit (HEAD): replace AI ``Co-authored-by:``
        trailers with ``Assisted-by:``, and make you the author if an AI tool
        is. To fix every commit on a branch, run it from ``git rebase --exec``
        (see the message printed by ``check``).
"""

import os
import re
import subprocess
import sys

# Email addresses used by AI tools in commit trailers or as commit authors.
AI_EMAIL = re.compile(
    r"(@anthropic\.com"
    r"|@cursor\.(sh|com)"
    r"|@openai\.com"
    r"|^copilot@github\.com"
    r"|\+copilot@users\.noreply\.github\.com)$",
    re.IGNORECASE,
)

# Names used by AI tools. A plain "Claude" is not matched, so that people
# called Claude are not flagged; the email address catches those cases.
AI_NAME = re.compile(
    r"^(claude (code|opus|sonnet|haiku|fable|mythos)\b.*"
    r"|(github )?copilot( agent)?"
    r"|cursor( agent)?)$",
    re.IGNORECASE,
)

TRAILER = re.compile(
    r"^co-authored-by:[ \t]*(?P<name>.*?)[ \t]*<(?P<email>[^>]*)>[ \t]*$",
    re.IGNORECASE | re.MULTILINE,
)

POLICY_URL = "https://xpsi-group.github.io/xpsi/ai_policy.html"


def is_ai(name, email):
    return bool(AI_EMAIL.search(email.strip()) or AI_NAME.match(name.strip()))


def git(*args, stdin=None):
    return subprocess.run(
        ["git", *args], check=True, text=True, capture_output=True, input=stdin
    ).stdout


def commit_info(rev):
    sep = "\x1f"
    out = git("log", "-1", f"--format=%an{sep}%ae{sep}%s{sep}%B", rev)
    name, email, subject, body = out.split(sep, 3)
    return name, email, subject, body


def ai_trailers(body):
    return [m for m in TRAILER.finditer(body) if is_ai(m["name"], m["email"])]


def escape(text):
    # Escape data for GitHub Actions workflow commands.
    return text.replace("%", "%25").replace("\r", "%0D").replace("\n", "%0A")


def check(base, head):
    commits = git("rev-list", "--reverse", f"{base}..{head}").split()
    problems = []
    for sha in commits:
        name, email, subject, body = commit_info(sha)
        issues = [
            f"lists `{m['name']} <{m['email']}>` as co-author"
            for m in ai_trailers(body)
        ]
        if is_ai(name, email):
            issues.append(f"has `{name} <{email}>` as its author")
        if issues:
            problems.append((sha[:7], subject, issues))

    if not problems:
        print(f"Checked {len(commits)} commit(s): no AI tools listed as "
              "co-author or author.")
        return 0

    lines = [
        "## AI tools listed as co-author or author",
        "",
        "X-PSI's [policy on generative AI](" + POLICY_URL + ") asks you to "
        "record AI assistance with an `Assisted-by:` line in the commit "
        "message, not `Co-authored-by:`, because GitHub counts co-authors "
        "as contributors to X-PSI. Some tools (for example Claude Code) add "
        "a `Co-authored-by:` line automatically.",
        "",
        "These commits need changing:",
        "",
    ]
    for short, subject, issues in problems:
        for issue in issues:
            lines.append(f"- `{short}` ({subject}): {issue}")
    lines += [
        "",
        "**To fix all commits on your branch**, run these commands in your "
        "clone, with your branch checked out (`upstream` is the "
        "`xpsi-group/xpsi` remote):",
        "",
        "```bash",
        "git fetch upstream",
        "git rebase upstream/main --exec "
        "\"python3 .github/scripts/ai_coauthor.py fix\"",
        "git push --force-with-lease",
        "```",
        "",
        "This replaces each AI `Co-authored-by:` line with `Assisted-by:`, "
        "and makes you the author of any commit authored by an AI tool. "
        "Rewriting history needs a force push; that is fine here.",
        "",
        "**If only your latest commit is affected**, you can instead run "
        "`git commit --amend`, change the line by hand to "
        "`Assisted-by: <tool>:<model>`, and push with "
        "`git push --force-with-lease`.",
        "",
        "To stop it happening again, see your tool's settings (in Claude "
        "Code, set `\"attribution\": {\"commit\": \"\"}` in "
        "`~/.claude/settings.json`).",
    ]
    report = "\n".join(lines)
    print(report)

    for short, subject, issues in problems:
        msg = f"Commit {short} " + "; ".join(issues).replace("`", "")
        print(f"::error title=AI tool listed as co-author::{escape(msg)}")

    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary:
        with open(summary, "a", encoding="utf-8") as fh:
            fh.write(report + "\n")
    return 1


def fix():
    name, email, _, body = commit_info("HEAD")

    def to_assisted(m):
        return f"Assisted-by: {m['name'] or 'AI tool'}"

    new_body = TRAILER.sub(
        lambda m: to_assisted(m) if is_ai(m["name"], m["email"]) else m[0],
        body,
    )
    reset_author = is_ai(name, email)
    if new_body == body and not reset_author:
        return 0

    cmd = ["commit", "--amend", "--no-verify", "--cleanup=whitespace",
           "-F", "-"]
    if reset_author:
        cmd.append("--reset-author")
    git(*cmd, stdin=new_body)
    print(f"Fixed commit: {git('log', '-1', '--format=%h %s').strip()}")
    return 0


def main(argv):
    if len(argv) == 4 and argv[1] == "check":
        return check(argv[2], argv[3])
    if len(argv) == 2 and argv[1] == "fix":
        return fix()
    print(__doc__, file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv))
