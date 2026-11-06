"""Check a caller-supplied release policy without publishing its accounting."""

import ast
import io
from pathlib import Path
import re
import subprocess
import tokenize

from .verifier import _read_json


def _python_lines(text):
    excluded = set()
    for node in ast.walk(ast.parse(text)):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            if ast.get_docstring(node, clean=False) is not None:
                doc = node.body[0]
                excluded.update(range(doc.lineno, doc.end_lineno + 1))
    syntax = {tokenize.COMMENT, tokenize.NL, tokenize.NEWLINE, tokenize.ENDMARKER,
              tokenize.INDENT, tokenize.DEDENT, tokenize.ENCODING}
    lines = set()
    for token in tokenize.generate_tokens(io.StringIO(text).readline):
        if token.type not in syntax:
            lines.update(range(token.start[0], token.end[0] + 1))
    return len(lines - excluded)


def private_policy_passes(path: Path, commit: str) -> bool:
    policy = _read_json(path, 65536)
    if policy["schema"] != 1:
        raise ValueError("unsupported release policy")
    required_lines, required_commits = policy["production_lines"], policy["meaningful_commits"]
    if not all(type(value) is int and value > 0 for value in (required_lines, required_commits)):
        raise ValueError("invalid release policy thresholds")
    tracked = subprocess.check_output(["git", "ls-files", "-z"], text=True).split("\0")
    count = 0
    for name in tracked:
        source = Path(name)
        suffix = source.suffix
        prefixes = {".py": "python/diachron/", ".rs": "rust/diachron-core/src/",
                    ".ts": "ui/src/", ".css": "ui/src/"}
        if suffix not in prefixes or not name.startswith(prefixes[suffix]):
            continue
        text = source.read_text(encoding="utf-8")
        if suffix == ".py":
            count += _python_lines(text)
        else:
            if suffix == ".rs":
                text = text.split("#[cfg(test)]", 1)[0]
            text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
            count += sum(bool(line.strip()) and not line.lstrip().startswith("//") for line in text.splitlines())
    commits = subprocess.check_output(["git", "log", "--format=%T", commit], text=True).splitlines()
    return count >= required_lines and len(set(commits)) >= required_commits
