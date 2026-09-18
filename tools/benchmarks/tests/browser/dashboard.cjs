const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const http = require('node:http');
const { chromium } = require('playwright');

(async () => {
  const root = path.resolve(process.argv[2]);
  const server = http.createServer((req, res) => {
    const target = path.resolve(root, '.' + decodeURIComponent(req.url.split('?')[0]));
    if (!target.startsWith(root + path.sep) || !fs.existsSync(target) || !fs.statSync(target).isFile()) { res.writeHead(404); res.end(); return; }
    res.setHeader('Content-Type', target.endsWith('.js') ? 'text/javascript' : target.endsWith('.css') ? 'text/css' : target.endsWith('.json') ? 'application/json' : 'text/html');
    res.end(fs.readFileSync(target));
  });
  await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
  let browser;
  try {
    browser = await chromium.launch({ headless: true });
    const page = await browser.newPage();
    const errors = [], failures = [];
    page.on('pageerror', e => errors.push(e.message));
    page.on('console', m => { if (m.type() === 'error') errors.push(m.text()); });
    page.on('requestfailed', r => failures.push(r.url()));
    const open = async name => {
      await page.goto(`http://127.0.0.1:${server.address().port}/${name}/index.html`);
      await page.waitForFunction(() => document.querySelector('#loadStatus').dataset.state);
    };
    await open('success');
    assert.equal(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth), true);
    await page.setViewportSize({width:390,height:844});
    assert.equal(await page.evaluate(() => document.documentElement.scrollWidth <= innerWidth), true);
    await page.setViewportSize({width:1280,height:900});
    assert.equal(await page.locator('#loadStatus').getAttribute('data-state'), 'loaded');
    assert.equal(await page.locator('#latestTable tbody tr').count(), 4);
    assert.match(await page.locator('#latestTable').innerText(), /3\.0000 s \(n=2\)/);
    assert.match(await page.locator('#buildMeasurements').innerText(), /n=1 build/);
    assert.equal(await page.locator('#bvhChart circle').count(), 4);
    await page.selectOption('#benchmarkFilter', 'bvh_many_spheres');
    assert.equal(await page.locator('#latestTable tbody tr').count(), 2);
    await page.selectOption('#configFilter', 'o3_native_no_simd');
    assert.equal(await page.locator('#latestTable tbody tr').count(), 1);
    await open('mixed');
    assert.match(await page.locator('#latestAttempt').innerText(), /attempt 2.*failed/);
    assert.match(await page.locator('#latestComplete').innerText(), /attempt 1.*complete/);
    assert.match(await page.locator('#selectionStatus').innerText(), /attempt 1/);
    assert.match(await page.locator('#diagnostics').innerText(), /Both builds failed/);
    await page.selectOption('#runFilter', { index: 0 });
    assert.match(await page.locator('#latestTable').innerText(), /build_failed/);
    assert.equal(await page.locator('#latestTable tbody').innerText().then(t => t.includes('3.0000 s')), false);
    await open('partial');
    assert.match(await page.locator('#latestComplete').innerText(), /No complete valid/);
    assert.match(await page.locator('#latestTable').innerText(), /3\.0000 s/);
    assert.match(await page.locator('#latestTable').innerText(), /build_failed/);
    await open('failed');
    assert.match(await page.locator('#renderChart').innerText(), /No successful/);
    await page.locator('#diagnostics details summary').first().click();
    assert.match(await page.locator('#diagnostics pre').first().innerText(), /<test> & diagnostics <\/script>/);
    await open('legacy');
    assert.match(await page.locator('#bvhChart').innerText(), /BVH timing not recorded/);
    assert.equal(await page.locator('#renderChart circle').count(), 4);
    await open('empty');
    assert.equal(await page.locator('#loadStatus').getAttribute('data-state'), 'empty');
    assert.equal(await page.locator('#latestTable tbody tr').count(), 0);
    assert.deepEqual(errors, []); assert.deepEqual(failures, []);
    // Exercise a real failed fetch, independently of missing benchmark metrics.
    await page.route('**/empty/index.html', async route => {
      const html = fs.readFileSync(path.join(root, 'empty/index.html'), 'utf8').replace(/(<script type="application\/json" id="benchmark-data">)[\s\S]*?(<\/script>)/, '$1$2');
      await route.fulfill({ contentType: 'text/html', body: html });
    });
    await page.route('**/data/dashboard.json', route => route.fulfill({ status: 503, body: 'unavailable' }));
    await open('empty');
    assert.equal(await page.locator('#loadStatus').getAttribute('data-state'), 'error');
    assert.match(await page.locator('#loadStatus').innerText(), /Dashboard load error: HTTP 503/);
    assert.equal(errors.filter(e => !e.includes('503')).length, 0);
    console.log('Browser smoke passed: successful, failed, partial, legacy, empty, filters, latest valid comparison, escaped diagnostics, HTTP error.');
  } finally {
    if (browser) await browser.close();
    await new Promise(resolve => server.close(resolve));
  }
})().catch(error => { console.error(error); process.exitCode = 1; });
