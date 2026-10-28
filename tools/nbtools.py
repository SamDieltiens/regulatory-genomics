"""Small helpers to build and execute the teaching notebooks.

Cells tagged "skip-execution" contain long-running commands (full download,
full alignment). These are run separately with the same commands and are not
re-executed when the notebook is executed.
"""
import re
import sys
import nbformat as nbf
from nbclient import NotebookClient

SKIP = "skip-execution"
ANSI = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")


def md(src):
    return nbf.v4.new_markdown_cell(src.strip("\n"))


def code(src, skip=False):
    c = nbf.v4.new_code_cell(src.strip("\n"))
    if skip:
        c.metadata["tags"] = [SKIP]
    return c


def new_notebook(cells, kernel, display_name, language):
    nb = nbf.v4.new_notebook()
    nb.cells = cells
    nb.metadata["kernelspec"] = {
        "name": kernel, "display_name": display_name, "language": language}
    return nb


def sync_markdown(path, nb):
    """Copy the markdown cells of nb into the executed notebook at path, keeping its outputs."""
    old = nbf.read(path, as_version=4)
    assert len(old.cells) == len(nb.cells), "cell structure changed, re-execute instead"
    for o, n in zip(old.cells, nb.cells):
        assert o.cell_type == n.cell_type, "cell structure changed, re-execute instead"
        if n.cell_type == "markdown":
            o.source = n.source
        else:
            assert o.source == n.source, "code changed, re-execute instead"
    return old


def execute(nb, cwd, timeout=7200):
    """Execute all cells that are not tagged with SKIP, in order."""
    client = NotebookClient(nb, kernel_name=nb.metadata["kernelspec"]["name"],
                            timeout=timeout, allow_errors=True, resources={"metadata": {"path": cwd}})
    with client.setup_kernel():
        n = 0
        for i, cell in enumerate(nb.cells):
            if cell.cell_type != "code":
                continue
            n += 1
            if SKIP in cell.metadata.get("tags", []):
                cell.outputs = []
                continue
            print(f"[{n}] {cell.source.splitlines()[0][:70]}", file=sys.stderr, flush=True)
            client.execute_cell(cell, i)
            for o in cell.outputs:  # strip terminal colour codes (ls/grep aliases)
                if o.output_type == "stream":
                    o.text = ANSI.sub("", o.text)
    return nb
