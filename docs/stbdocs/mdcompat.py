"""Makes GitHub-flavoured Markdown render the same way under Python-Markdown.

The example READMEs are written for GitHub (CommonMark) but the documentation
site renders them with Python-Markdown (MkDocs). The two disagree on how deep a
nested list is indented and on whether a list may follow a paragraph line
directly; `github_markdown_to_python_markdown` bridges both. test/docs/
compares the result against a CommonMark renderer for every README.
"""

import re

# A list item marker with the whitespace after it: `- `, `* `, `+ `, `1. `, `1) `.
_LIST_START = re.compile(r"^(\s*)([-*+]|\d+[.)])(\s+)\S")
# Lines that are never plain paragraph text: heading, table row, blockquote, HTML, rule.
_NOT_PARAGRAPH = re.compile(r"^\s*(?:#|\||>|<|(?:[-*_]\s*){3,}$)")
_QUOTE_LINE = re.compile(r"^ {0,3}>")
_FENCE = ("```", "~~~")


_CODE_SPAN = re.compile(r"(`+)(?!`)(.+?)(?<!`)\1(?!`)")


def _fix_escaped_pipes(line):
    """In a GitHub table row `\\|` is a literal pipe, including inside a code
    span; Python-Markdown would show the backslash there. A bare `|` inside a
    code span does not split a cell in Python-Markdown, so inside code the
    backslash is dropped, and outside code the pipe becomes an entity."""
    out, last = [], 0
    for match in _CODE_SPAN.finditer(line):
        out.append(line[last:match.start()].replace("\\|", "&#124;"))
        out.append(match.group(0).replace("\\|", "|"))
        last = match.end()
    out.append(line[last:].replace("\\|", "&#124;"))
    return "".join(out)


def _indent(line):
    return len(line) - len(line.lstrip())


def _content_offset(match):
    """Column where a list item's text starts (what its continuation lines and
    nested lists must be indented to, per CommonMark)."""
    lead, marker, gap = match.group(1), match.group(2), match.group(3)
    return len(lead) + len(marker) + (len(gap) if len(gap) <= 4 else 1)


def _normalize_lists(text):
    """Rewrites the list structure of GitHub-flavoured Markdown (CommonMark) so
    that Python-Markdown (MkDocs) renders it the same way. Two differences:

    - Nesting. CommonMark nests a list under an item when it is indented to
      the item's text (2 spaces under `- `, 3 under `1. `); Python-Markdown
      only nests at 4 spaces per level. Every list level, and every
      continuation line inside an item, is re-indented to 4 spaces per level.
    - Interrupting a paragraph. CommonMark lets a top-level list follow a
      paragraph line directly; Python-Markdown then folds it into the
      paragraph as literal text, so a blank line is inserted first. (An
      ordered list interrupts a paragraph only if it starts at 1.)

    Table rows also get their escaped pipes (`\\|`) fixed. Fenced code blocks
    are shifted along with the item they belong to and otherwise left alone."""
    out, stack = [], []          # stack: content offsets of the open list items
    prev_kind, prev_blank = None, True
    fence, fence_delta = None, 0
    for raw in text.splitlines(keepends=True):
        stripped = raw.strip()
        if fence is None and stripped.startswith("|") and "\\|" in raw:
            raw = _fix_escaped_pipes(raw)
        if fence is not None:
            if fence_delta >= 0:
                out.append(" " * fence_delta + raw if stripped else raw)
            else:
                out.append(raw[min(-fence_delta, _indent(raw)):] if stripped else raw)
            if stripped.startswith(fence):
                fence = None
                prev_kind, prev_blank = "special", False
            continue
        if not stripped:
            out.append(raw)
            prev_blank = True
            continue

        indent = _indent(raw)
        in_list_before = bool(stack)
        while stack and indent < stack[-1]:
            stack.pop()
        depth = len(stack)
        match = _LIST_START.match(raw)

        if match:
            interrupts = depth == 0 and prev_kind == "para" \
                and (not match.group(2)[0].isdigit() or match.group(2)[:-1] == "1")
            if interrupts:
                out.append("\n")
            out.append(" " * (4 * depth) + raw.lstrip())
            stack.append(_content_offset(match))
            prev_kind = "list"
        else:
            new_indent = indent if depth == 0 else 4 * depth + max(indent - stack[-1], 0)
            out.append(" " * new_indent + raw.lstrip())
            if stripped.startswith(_FENCE):
                fence, fence_delta = stripped[:3], new_indent - indent
                prev_kind = "special"
            elif depth > 0 or (in_list_before and not prev_blank):
                prev_kind = "list"    # item content, or a lazy continuation of one
            elif _NOT_PARAGRAPH.match(raw):
                prev_kind = "special"
            else:
                prev_kind = "para"
        prev_blank = False
    return "".join(out)


def github_markdown_to_python_markdown(text):
    """Applies `_normalize_lists` to the document and, separately, to the body
    of every blockquote (a list after a paragraph line inside a `>` block needs
    the same fix). Fenced code blocks, including any `>` inside them, are
    left untouched."""
    out, run, quoting, in_fence = [], [], False, False

    def flush():
        if quoting:
            inner = [re.sub(r"^ {0,3}> ?", "", line, count=1) for line in run]
            for line in _normalize_lists("".join(inner)).splitlines(keepends=True):
                out.append("> " + line if line.strip() else ">\n")
        else:
            out.append(_normalize_lists("".join(run)))

    for line in text.splitlines(keepends=True):
        is_quote = not in_fence and bool(_QUOTE_LINE.match(line))
        if run and is_quote != quoting:
            flush()
            run = []
        quoting = is_quote
        run.append(line)
        if not is_quote and line.lstrip().startswith(_FENCE):
            in_fence = not in_fence
    if run:
        flush()
    return "".join(out)
