#!/usr/bin/env node
/**
 * Packages BaitBench for distribution.
 *
 * Produces dist-release/ containing:
 *   macOS:   BaitBench.app  (ad-hoc signed)  +  baitbench  (CLI binary)
 *   Linux:   BaitBench.AppImage              +  baitbench
 *   Windows: BaitBench_*.msi / setup.exe     +  baitbench.exe
 *
 * Usage (called by `make release`):
 *   node scripts/package-release.mjs
 */

import { execSync, spawnSync } from 'child_process';
import {
  cpSync, mkdirSync, existsSync, readdirSync,
  copyFileSync, chmodSync,
} from 'fs';
import { resolve, dirname, join, extname, basename } from 'path';
import { fileURLToPath } from 'url';
import { platform } from 'os';

const __dirname = dirname(fileURLToPath(import.meta.url));
const repoRoot   = resolve(__dirname, '../..');
const guiDir     = resolve(__dirname, '..');
const bundleRoot = resolve(guiDir, 'src-tauri/target/release/bundle');
const outDir     = resolve(guiDir, 'dist-release');

const os = platform(); // 'darwin' | 'linux' | 'win32'

// ── Helpers ─────────────────────────────────────────────────────────────────

function run(cmd, opts = {}) {
  console.log(`  $ ${cmd}`);
  execSync(cmd, { stdio: 'inherit', ...opts });
}

function findFirst(dir, predicate) {
  if (!existsSync(dir)) return null;
  for (const entry of readdirSync(dir, { withFileTypes: true })) {
    const full = join(dir, entry.name);
    if (predicate(entry, full)) return full;
    if (entry.isDirectory()) {
      const found = findFirst(full, predicate);
      if (found) return found;
    }
  }
  return null;
}

// ── Prepare output dir ───────────────────────────────────────────────────────

mkdirSync(outDir, { recursive: true });
console.log(`\n📦  Packaging release into ${outDir}\n`);

// ── Copy CLI binary ───────────────────────────────────────────────────────────

const cliSrc  = resolve(repoRoot, os === 'win32' ? 'target/release/baitbench.exe' : 'target/release/baitbench');
const cliDest = resolve(outDir,   os === 'win32' ? 'baitbench.exe' : 'baitbench');

if (!existsSync(cliSrc)) {
  console.error(`CLI binary not found at ${cliSrc}. Run 'cargo build --release' first.`);
  process.exit(1);
}
copyFileSync(cliSrc, cliDest);
if (os !== 'win32') chmodSync(cliDest, 0o755);
console.log(`✓  CLI binary  →  ${cliDest}`);

// ── macOS: copy .app and ad-hoc sign ─────────────────────────────────────────

if (os === 'darwin') {
  const appSrc = findFirst(join(bundleRoot, 'macos'), (e) => e.name.endsWith('.app') && e.isDirectory());
  if (!appSrc) {
    console.error('macOS .app bundle not found. Did `npm run tauri:build` succeed?');
    process.exit(1);
  }
  const appDest = join(outDir, basename(appSrc));
  cpSync(appSrc, appDest, { recursive: true });
  console.log(`✓  App bundle  →  ${appDest}`);

  // Ad-hoc sign: removes the quarantine "damaged app" error for unsigned distributions.
  // This does NOT provide notarization — users will still see an Gatekeeper prompt
  // on first launch (right-click → Open), but the app won't appear as corrupted.
  const codesign = spawnSync('codesign', [
    '--force', '--deep', '--sign', '-',
    '--options', 'runtime',
    appDest,
  ], { stdio: 'inherit' });
  if (codesign.status === 0) {
    console.log('✓  Ad-hoc signed');
  } else {
    console.warn('⚠  codesign not available or failed — app will not be ad-hoc signed');
  }

  // Also copy the .dmg if present
  const dmgSrc = findFirst(join(bundleRoot, 'dmg'), (e) => e.name.endsWith('.dmg'));
  if (dmgSrc) {
    const dmgDest = join(outDir, basename(dmgSrc));
    copyFileSync(dmgSrc, dmgDest);
    console.log(`✓  DMG installer  →  ${dmgDest}`);
  }
}

// ── Linux: copy AppImage ──────────────────────────────────────────────────────

if (os === 'linux') {
  const appimage = findFirst(join(bundleRoot, 'appimage'), (e) => e.name.endsWith('.AppImage'));
  if (appimage) {
    const dest = join(outDir, basename(appimage));
    copyFileSync(appimage, dest);
    chmodSync(dest, 0o755);
    console.log(`✓  AppImage  →  ${dest}`);
  } else {
    console.warn('⚠  AppImage not found — check the tauri build output');
  }
}

// ── Windows: copy installer ───────────────────────────────────────────────────

if (os === 'win32') {
  for (const subdir of ['msi', 'nsis']) {
    const installer = findFirst(join(bundleRoot, subdir), (e) =>
      e.name.endsWith('.msi') || e.name.endsWith('.exe')
    );
    if (installer) {
      const dest = join(outDir, basename(installer));
      copyFileSync(installer, dest);
      console.log(`✓  Installer  →  ${dest}`);
    }
  }
}

// ── Summary ───────────────────────────────────────────────────────────────────

console.log(`
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  Release artifacts in:  ${outDir}/
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

Contents:`);
for (const f of readdirSync(outDir)) {
  console.log(`  ${f}`);
}

if (os === 'darwin') {
  console.log(`
macOS distribution notes
────────────────────────
• BaitBench.app  — drag to /Applications and double-click to launch.
• baitbench      — standalone CLI; copy to /usr/local/bin or ~/bin.

Gatekeeper (first launch):
  Because this build is not notarized, macOS will block the first open.
  Users should right-click → Open (or run the one-time command below):
    xattr -d com.apple.quarantine BaitBench.app

For fully seamless distribution (no prompt), sign and notarize with an
Apple Developer account:
    xcrun notarytool submit BaitBench.dmg --apple-id ... --notarize
`);
}
