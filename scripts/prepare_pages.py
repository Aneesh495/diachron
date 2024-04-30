"""Add navigation to copies of the checksum-verified standalone release files."""
import sys
from pathlib import Path


def prepare(site: Path) -> None:
    preview = site / "Diachron_Preview.html"
    thesis = site / "Diachron_Thesis.html"
    css = """<style>
.diachron-links { display: flex; flex-wrap: wrap; gap: 1rem; padding: .8rem 1rem;
  color: #e0edf7; background: #142030; font: 600 .85rem/1.4 system-ui, sans-serif; }
.diachron-links a { color: #96ddff; }
@media (max-width: 48rem) {
  pre { overflow-x: auto; max-width: 100%; }
  :not(pre) > code { overflow-wrap: anywhere; }
  figure { margin-inline: 0; max-width: 100%; }
  img { max-width: 100%; height: auto; }
  table { display: block; max-width: 100%; overflow-x: auto; }
}
</style>"""

    for source, target, href, label in [
        (preview, site / "index.html", "Diachron_Thesis.html", "Research report"),
        (thesis, thesis, "./", "Open the workbench"),
    ]:
        original = source.read_text(encoding="utf-8")
        body = original.find("<body")
        start = original.index(">", body) + 1 if body >= 0 else original.index("<h1")
        nav = (f'<nav class="diachron-links" aria-label="Diachron pages"><a href="{href}">{label}</a>'
               '<a href="https://github.com/Aneesh495/diachron">Source and reproducibility</a></nav>')
        output = original[:start] + nav + original[start:]
        if "</head>" in output:
            output = output.replace("</head>", css + "</head>", 1)
        else:
            output = output[:start] + css + output[start:]
        target.write_text(output, encoding="utf-8")

    # The published thesis copy adds navigation. Its release checksum belongs to
    # the original release asset, which remains available there.
    (site / "SHA256SUMS").unlink()


if __name__ == "__main__":
    prepare(Path(sys.argv[1]))
