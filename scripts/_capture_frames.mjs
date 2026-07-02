// Headless frame capture for an animated GIF of a reglscatterpy save_html plot.
// Drives the plot's public API (headless can't fake a mouse drag) and screenshots
// each frame. Params via a JSON file argv[3].
//   node _capture_frames.mjs <html> <outdir> <params.json>
import { createRequire } from "module";
import { mkdirSync } from "fs";
const require = createRequire("/home/goguxor/Desktop/reglScatterplotR/js/");
const puppeteer = require("puppeteer-core");

const [html, outdir, paramsPath] = process.argv.slice(2);
const P = require(paramsPath);                 // {w,h,select,holdStart,zoomIn,holdSel,zoomOut,holdEnd,padding}
mkdirSync(outdir, { recursive: true });

const b = await puppeteer.launch({ executablePath: "/usr/bin/chromium", headless: "new",
  args: ["--use-gl=angle","--use-angle=swiftshader","--no-sandbox","--ignore-gpu-blocklist","--force-device-scale-factor=1"] });
const page = await b.newPage();
await page.setViewport({ width: P.w, height: P.h, deviceScaleFactor: 2 });
await page.goto("file://" + html, { waitUntil: "networkidle0" });
await new Promise(r => setTimeout(r, 3000));   // widget init + first draw

const getPlot = `(() => { const reg = window.__myScatterplotRegistry; return reg && Array.from(reg.values())[0]; })()`;

// initial + target cameraView (target = frame the selected points, or a center box)
const cams = await page.evaluate((sel, pad) => {
  const reg = window.__myScatterplotRegistry; const e = Array.from(reg.values())[0];
  const cam0 = Array.from(e.plot.get("cameraView"));
  if (sel && sel.length) e.plot.zoomToPoints(sel, { transition: false, padding: pad });
  else e.plot.zoomToArea({ x: -0.35, y: -0.35, width: 0.7, height: 0.7 }, { transition: false });
  const cam1 = Array.from(e.plot.get("cameraView"));
  e.plot.set({ cameraView: cam0 });            // reset
  return { cam0, cam1 };
}, P.select || null, P.padding ?? 0.3);

const lerp = (a, b, t) => a.map((v, i) => v + (b[i] - v) * t);
const ease = t => t < .5 ? 2*t*t : 1 - Math.pow(-2*t+2, 2)/2;   // easeInOutQuad
let f = 0;
const shot = async () => { await page.screenshot({ path: `${outdir}/f${String(f++).padStart(3,"0")}.png` }); };

async function setCam(c, selected) {
  await page.evaluate((cam, sel) => {
    const e = Array.from(window.__myScatterplotRegistry.values())[0];
    e.plot.set({ cameraView: cam });
    if (sel === true && window.__rs_sel) e.plot.select(window.__rs_sel, { preventEvent: true });
    if (sel === false) e.plot.deselect({ preventEvent: true });
  }, c, selected);
  await new Promise(r => setTimeout(r, 60));   // let swiftshader paint
}

// stash selection indices in-page
if (P.select) await page.evaluate(s => { window.__rs_sel = s; }, P.select);

for (let i = 0; i < P.holdStart; i++) { await setCam(cams.cam0); await shot(); }
for (let i = 1; i <= P.zoomIn; i++)   { await setCam(lerp(cams.cam0, cams.cam1, ease(i/P.zoomIn))); await shot(); }
if (P.select) await setCam(cams.cam1, true);
for (let i = 0; i < P.holdSel; i++)   { await setCam(cams.cam1, i===0 ? true : undefined); await shot(); }
if (P.select) await page.evaluate(() => Array.from(window.__myScatterplotRegistry.values())[0].plot.deselect({preventEvent:true}));
for (let i = 1; i <= P.zoomOut; i++)  { await setCam(lerp(cams.cam1, cams.cam0, ease(i/P.zoomOut))); await shot(); }
for (let i = 0; i < P.holdEnd; i++)   { await setCam(cams.cam0); await shot(); }

console.log("frames:", f);
await b.close();
