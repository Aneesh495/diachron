import assert from "node:assert/strict";
import { mkdir, readFile, writeFile } from "node:fs/promises";
import { createHash } from "node:crypto";
import { resolve } from "node:path";
import { pathToFileURL } from "node:url";
import puppeteer from "puppeteer";

const [input, directory] = process.argv.slice(2);
if (!input || !directory) throw new Error("Supply a standalone HTML file and an output directory.");
const filename = resolve(input);
const url = pathToFileURL(filename).href;
await mkdir(directory, { recursive: true });
const browser = await puppeteer.launch({ headless: true, args: ["--no-sandbox"] });
const failures = [];
const requests = [];
try {
  async function open(viewport) {
    const page = await browser.newPage();
    page.on("pageerror", error => failures.push(String(error)));
    page.on("request", request => requests.push(request.url()));
    await page.setViewport(viewport);
    await page.goto(url, { waitUntil: "load" });
    await page.waitForSelector(".research-tour");
    return page;
  }
  const page = await open({ width: 1440, height: 920 });
  const modelSelector = 'select[aria-label="Posterior model"]';
  const regional = await page.$(`${modelSelector} option[value="regional"]`) !== null;
  const modelCount = regional ? 5 : 4;
  assert.match(await page.$eval(".wb-mixing", el => el.textContent), /Diagnostic draws|Recorded scalar checks passed/);
  assert.equal(await page.$$(".wb-panel").then(nodes => nodes.length), 4);
  const next = ".tour-navigation button:last-child";
  const nextStep = async (title) => {
    await page.click(next);
    assert.match(await page.$eval(".research-tour h2", el => el.textContent), title);
  };
  await page.click(".research-tour summary");
  await page.click(".research-tour div:last-child button");
  let readout = await page.$eval(".matrix-readout", el => el.textContent);
  assert.match(readout, /Unknown/);
  assert.match(readout, /Coded absence/);
  await nextStep(/alternate cognacy/);
  await page.click(".research-tour div:last-child button");
  assert.equal(await page.$$(".matrix-readout li").then(nodes => nodes.length), 2);
  await nextStep(/loan annotation/);
  await page.click(".research-tour div:last-child button");
  readout = await page.$eval(".matrix-readout", el => el.textContent);
  assert.match(readout, /4364/);
  assert.match(readout, /227/);
  await nextStep(/Knossos/);
  const dates = await page.$$(".research-tour p button");
  assert.equal(dates.length, 2);
  await dates[0].click();
  assert.match(await page.$eval(".timeline-readout", el => el.textContent), /3399/);
  await dates[1].click();
  assert.match(await page.$eval(".timeline-readout", el => el.textContent), /3199/);
  await nextStep(/location modes/);
  await page.click(".research-tour div:last-child button");
  assert.equal(await page.$eval(modelSelector, el => el.value), "geography");
  const firstChain = await page.$eval(".map-view circle[data-chain]", el => el.dataset.chain);
  await page.select('select[aria-label="Location draw chain"]', firstChain);
  const visible = await page.$$eval(".map-view circle[data-chain]", nodes => nodes.map(el => el.dataset.chain));
  assert(visible.length > 0 && visible.every(chain => chain === firstChain));
  await page.select('select[aria-label="Location draw chain"]', "all");
  await page.select('select[aria-label="Geographic projection"]', "globe");
  await page.waitForSelector('[data-globe="visible"]');
  await page.click('button[aria-label="Rotate east"]');
  await page.click('button[aria-label="Rotate west"]');
  await page.click(".research-tour summary");
  await page.click(".taxonomy-comparison summary");
  assert.match(await page.$eval(".taxonomy-comparison ul", el => el.textContent), /4800 draws/);
  assert.equal(await page.$$(".taxonomy-comparison li").then(nodes => nodes.length), modelCount);
  await page.click(".taxonomy-comparison button");
  await page.waitForSelector('.map-view rect[stroke="#182739"]');
  await page.click(".taxonomy-comparison summary");
  await page.click(".tree-controls button:nth-child(2)");
  await page.click(".model-comparison summary");
  assert.equal(await page.$$(".model-comparison tbody tr").then(nodes => nodes.length), modelCount);
  await page.click(".model-comparison summary");
  await page.select('select[aria-label="Geographic projection"]', "map");
  await page.click(".attestation-inspector summary");
  assert.match(await page.$eval(".attestation-inspector a", el => el.href), /\/HD\d{6}$/);
  await page.click(".attestation-inspector button");
  assert.match(await page.$eval(".timeline-readout", el => el.textContent), /before 2000 CE/);
  assert.match(await page.$eval(".attestation-inspector p", el => el.textContent), /before 1950 CE/);
  assert.equal(await page.evaluate(() => document.activeElement.getAttribute("aria-label")), "Select an attestation");
  await page.waitForSelector('.wb-map path[fill="#e6a23c"]');
  await page.click(".attestation-inspector summary");
  await page.screenshot({ path: resolve(directory, "workbench-desktop.png") });

  if (regional) {
    await page.click(".tree-controls button:nth-child(2)");
    await page.select(modelSelector, "regional");
    await page.evaluate(() => {
      const data = JSON.parse(document.getElementById("workbench-data").textContent);
      const age = Math.min(...data.bundles.regional.rootLocations.map(point => point.age));
      const slider = document.querySelector('input[aria-label="Years before present"]');
      slider.value = String(Math.round(age));
      slider.dispatchEvent(new Event("input"));
    });
    await page.waitForSelector("circle[data-location-region]");
    assert.match(await page.$eval("circle[data-location-region] title", el => el.textContent), /region centroid\. \d+ of \d+ draws/);
    assert.match(await page.$eval(".map-view", el => el.textContent), /sampled region centroids/);
    await page.screenshot({ path: resolve(directory, "workbench-regional.png") });
  }

  const mobile = await open({ width: 390, height: 844, isMobile: true, hasTouch: true });
  await mobile.click(".research-tour summary");
  assert.equal(await mobile.evaluate(() => document.documentElement.scrollWidth), 390);
  const timeline = 'input[aria-label="Years before present"]';
  await mobile.focus(timeline);
  const before = await mobile.$eval(".timeline-readout", el => el.textContent);
  await mobile.keyboard.press("ArrowRight");
  assert.notEqual(await mobile.$eval(".timeline-readout", el => el.textContent), before);
  await mobile.screenshot({ path: resolve(directory, "workbench-mobile.png"), fullPage: true });
  assert.deepEqual(failures, []);
  assert(requests.every(request => request.startsWith("file:") || request.startsWith("data:")));
  const html = await readFile(filename);
  await writeFile(resolve(directory, "workbench-qa.json"), JSON.stringify({
    passed: true, browser: await browser.version(), desktop: [1440, 920], mobile: [390, 844],
    html_sha256: createHash("sha256").update(html).digest("hex"),
    verified: ["offline standalone HTML", "four linked panels", "five evidence examples",
      "two sourced calendar alternatives", "location chain filter", "catalog group frequencies and selection",
      "globe rotation", "four-model comparison", "EDH epoch conversion and focus", "mobile keyboard timeline",
      ...(regional ? ["regional centroid frequencies on saved draws"] : [])],
    network_requests: 0,
  }, null, 2) + "\n");
} finally {
  await browser.close();
}
