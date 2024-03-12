# Saved-draw workbench

The workbench opens a bounded export of saved draws and their evidence. Its
model banner reports the convergence outcome for the full fit. Tree summaries
and location displays use a deterministic subset of at most 120 draws per
model. Catalog comparisons count exact sampled clades over every retained draw.

Export model bundles from the same frozen study. Supply the pinned lexical
source and catalog archive for the research tour and taxonomy comparisons.

```sh
.venv/bin/diachron export-workbench \
  --bundle tree_only=runs/inference/tree_only/tree_only \
  --bundle dated=runs/inference/dated/dated \
  --bundle geography=runs/inference/geography/geography \
  --bundle borrowing=runs/inference/borrowing/borrowing \
  --source runs/inputs/iecor --catalog runs/inputs/glottolog53.zip \
  --out runs/workbench/workbench.json
npm --prefix ui ci --ignore-scripts
npm --prefix ui run build
python3 scripts/bundle_workbench_html.py \
  --data runs/workbench/workbench.json --out runs/workbench/Diachron.html
```

Open the HTML directly. It contains the compiled interface and evidence snapshot.
The browser reads the embedded JSON and uses no network requests for exploration.
External source links remain available for deliberate inspection.

The tour uses matrix states rebuilt from pinned lexical bytes, including source
proposed cognacy and the recorded fruit loan. The Knossos case compares two
published historical assignments of an archive date. Mycenaean Greek is outside
this frozen lexical cohort. The workbench shows approximate calendar markers
alongside source citations. The IE-CoR terminal date likelihood remains defined
by its own input. [Hood's earlier-date argument](https://doi.org/10.1017/S0003598X00029550)
[Palmer's find-place analysis](https://doi.org/10.1017/S0003598X00036024)

Compare geographic draws by chain to examine differences behind the pooled
cloud. The protocol 3 fits leave the number of historical modes unresolved.
The chain filter reveals that dependence without changing the underlying draws.
Catalog groups use only the members present in the study and distinguish the
IE-CoR source classification from Glottolog 5.3. Their frequencies describe
the saved sample.

With Puppeteer Chromium available, run the browser check against the standalone
file. The report includes the exact HTML digest and the desktop and mobile
screenshots.

The observed four-model snapshot passed [browser verification](results/workbench-full-browser.json).
The [desktop](screenshots/workbench-observed-full.png) and
[mobile](screenshots/workbench-observed-full-mobile.png) screenshots show the
same offline export.

Add a verified regional bundle as `--bundle regional=PATH` to compare it
with the planar geographic fit. Regional marks show sampled category
frequencies at centroids in the selected time window. The legend and marker
details identify the aggregation. [Regional browser replay](results/workbench-regional-browser.json)
[Regional screenshot](screenshots/workbench-regional.png)

```sh
node ui/qa/workbench.mjs runs/workbench/Diachron.html runs/workbench/qa
```
