import assert from "node:assert/strict";
import { fileURLToPath } from "node:url";
import { build } from "esbuild";

const built = await build({
  entryPoints: [fileURLToPath(new URL("../src/posterior.ts", import.meta.url))],
  bundle: true, platform: "node", format: "esm", write: false,
});
const source = Buffer.from(built.outputFiles[0].contents).toString("base64");
const { PosteriorSummary } = await import(`data:text/javascript;base64,${source}`);
const majority = { newick: "((a:1,b:1):1,(c:1,d:1):1);", root_age: 2 };
const first = { newick: "((a:1,c:1):1,(b:1,d:1):1);", root_age: 2 };
const summary = new PosteriorSummary({ tip_labels: ["a", "b", "c", "d"], chains: [[first, majority, majority]] });
assert.equal(summary.maxCredibilityTree(), summary.trees[1]);
assert.equal(summary.supportForTaxa(new Set([0, 1])), 2 / 3);
console.log("Maximum clade credibility chooses the supported draw.");
