# BaitBench

BaitBench is an in-silico probe capture simulation tool for evaluating how well a probeset performs. It answers questions like:

- Does the probeset capture all target sequences?
- Does it reject background (distractor) sequences?
- Can it discriminate between organisms within the target panel?
- How does performance change at different target abundances (CT values)?
- What sequencing depth is needed for adequate genome coverage?

---

## Documentation

BaitBench's documentation is organized around four kinds of content:

| Section | Purpose |
|---------|---------|
| [**Tutorials**](tutorials/index.md) | Step-by-step walkthroughs for newcomers. Start here. |
| [**How-To Guides**](how-to/index.md) | Goal-oriented guides for specific tasks (input prep, parameter tuning, interpreting results). |
| [**Reference**](reference/index.md) | Complete, precise information: every subcommand, flag, output column, and file format. |
| [**Explanation**](explanation/index.md) | Background and concepts: the thermodynamic model, classification system, pipeline design. |

---

## Quick Start

<style>
.qs-wrap { border: 2px solid #b0bec5; border-radius: 8px; overflow: hidden; margin: 1.5em 0; }
.qs-tabs { display: flex; background: #dde4ec; border-bottom: 2px solid #b0bec5; }
.qs-tab {
  padding: 16px 36px; cursor: pointer; border: none; background: none;
  font-size: 1.4rem; font-weight: 700; color: #455a64;
  border-right: 1px solid #b0bec5;
}
.qs-tab:last-child { border-right: none; }
.qs-tab:hover { background: #cdd6e0; color: #1a252f; }
.qs-tab.active { background: #fff; color: #1565c0; border-bottom: 4px solid #1565c0; margin-bottom: -2px; }
.qs-panel { display: none; padding: 28px 32px; background: #fff; line-height: 1.6; }
.qs-panel.active { display: block; }
.qs-panel p { margin: 0.8em 0; }
.qs-panel h4 { margin: 1.6em 0 0.5em; font-size: 1.3rem; font-weight: 700; border-bottom: 1px solid #eee; padding-bottom: 4px; }
.qs-panel hr { border: none; border-top: 1px solid #ddd; margin: 1.6em 0; }
.qs-panel ol, .qs-panel ul { padding-left: 1.6em; margin: 0.6em 0; }
.qs-panel li { margin: 0.4em 0; }
.qs-panel pre { background: #f4f6f8; border: 1px solid #d0d7de; border-radius: 6px; padding: 16px 18px; overflow-x: auto; margin: 0.8em 0; }
.qs-panel pre, .qs-panel pre * { font-family: 'Menlo', 'Consolas', 'Monaco', monospace; color: #24292f !important; background: none !important; }
.qs-panel code { font-family: 'Menlo', 'Consolas', 'Monaco', monospace; color: #24292f !important; background: #f0f3f5 !important; border: 1px solid #d0d7de; border-radius: 4px; padding: 1px 5px; }
</style>

<div class="qs-wrap">
<div class="qs-tabs">
  <button class="qs-tab active" onclick="qsTab(event,'cli')">Command Line</button>
  <button class="qs-tab" onclick="qsTab(event,'mac')">Mac GUI</button>
  <button class="qs-tab" onclick="qsTab(event,'win')">Windows GUI</button>
</div>

<div id="qs-cli" class="qs-panel active">
<pre><code># 1. Clone the repository
git clone https://github.com/niel-infante/BaitBench.git
cd BaitBench

# 2. Install dependencies (requires conda)
conda env create -f environment.yml
conda activate baitbench

# 3. Build (requires the Rust toolchain — https://rustup.rs)
cargo build --release

# 4. Run a simulation
./target/release/baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --outdir results/</code></pre>
<p>See <a href="tutorials/first-run.md">Your First Simulation</a> for a complete walkthrough.</p>
</div>

<div id="qs-mac" class="qs-panel">
<p><strong>Download:</strong> <a href="https://github.com/niel-infante/BaitBench/releases/latest">BaitBench-macOS.dmg — latest release</a> (Apple Silicon only)</p>
<p>The app is currently <strong>unsigned and unnotarized</strong>, so macOS will block it on first launch. Steps to allow it are below.</p>
<hr>
<h4>Step 1 — Install BaitBench</h4>
<ol>
  <li>Open your <strong>Downloads</strong> folder and double-click <strong>BaitBench-macOS.dmg</strong>. A window will open showing the BaitBench icon and an Applications shortcut.</li>
  <li>Drag the <strong>BaitBench</strong> icon onto the <strong>Applications</strong> folder shortcut in that window to install it.</li>
  <li>You can now eject the disk image by pressing the Eject button next to it in the Finder sidebar, or simply close the window.</li>
</ol>
<hr>
<h4>Step 2 — Allow BaitBench to open</h4>
<p>Because the app is unsigned, macOS will block it the first time. The steps differ slightly by macOS version.</p>
<p><strong>macOS 15 Sequoia and later</strong></p>
<ol>
  <li>Double-click <strong>BaitBench</strong> in your Applications folder. You will see <em>"BaitBench can't be opened because Apple cannot check it for malicious software."</em> Click <strong>Done</strong> (not Move to Trash).</li>
  <li>Open <strong>System Settings → Privacy &amp; Security</strong>.</li>
  <li>Scroll down to the Security section. You should see <em>"BaitBench was blocked from use because it is not from an identified developer."</em></li>
  <li>Click <strong>Open Anyway</strong> and authenticate with Touch ID or your password.</li>
  <li>Double-click <strong>BaitBench</strong> again — it will now open. You only need to do this once.</li>
</ol>
<p><strong>macOS 13 Ventura and macOS 14 Sonoma</strong></p>
<ol>
  <li><strong>Right-click</strong> (or Control-click) <strong>BaitBench</strong> in your Applications folder and choose <strong>Open</strong>.</li>
  <li>Click <strong>Open</strong> in the dialog that appears. You only need to do this once.</li>
  <li>If right-click does not work, follow the Sequoia steps above — they work on all versions.</li>
</ol>
<p><strong>Optional: Terminal shortcut (any macOS version)</strong></p>
<p>If you are comfortable with the Terminal app, you can remove the quarantine flag directly instead of going through System Settings. Open Terminal and run:</p>
<pre><code>xattr -dr com.apple.quarantine /Applications/BaitBench.app</code></pre>
<p>Then double-click <strong>BaitBench</strong> normally. You only need to do this once.</p>
<hr>
<h4>Step 3 — First launch setup</h4>
<p>The first time BaitBench opens it will check for the pipeline tools it needs. You will be offered two options:</p>
<ul>
  <li><strong>Set up automatically</strong> — BaitBench will install <a href="https://www.anaconda.com/docs/getting-started/installation">conda</a> if needed and create the <code>baitbench</code> environment for you. This takes 5–10 minutes.</li>
  <li><strong>Use an existing environment</strong> — if you already have a <code>baitbench</code> conda environment, BaitBench will detect it automatically or let you browse to it.</li>
</ul>
<p>After setup completes, BaitBench is ready to use.</p>
</div>

<div id="qs-win" class="qs-panel">
<p><strong>Download:</strong> <a href="https://github.com/niel-infante/BaitBench/releases/latest">BaitBench-Windows.msi — latest release</a></p>
<p>The installer is currently <strong>unsigned</strong>, so Windows SmartScreen will warn you on first run.</p>
<hr>
<h4>Step 1 — Install BaitBench</h4>
<ol>
  <li>Open your <strong>Downloads</strong> folder and double-click <strong>BaitBench-Windows.msi</strong>.</li>
  <li>You will see a SmartScreen warning: <em>"Windows protected your PC"</em>. Click <strong>More info</strong>, then click <strong>Run anyway</strong>.</li>
  <li>Follow the installer prompts. BaitBench will be added to your Start menu.</li>
</ol>
<hr>
<h4>Step 2 — First launch setup</h4>
<p>The first time BaitBench opens it will check for the pipeline tools it needs. You will be offered two options:</p>
<ul>
  <li><strong>Set up automatically</strong> — BaitBench will install <a href="https://www.anaconda.com/docs/getting-started/installation">conda</a> if needed and create the <code>baitbench</code> environment for you. This takes 5–10 minutes.</li>
  <li><strong>Use an existing environment</strong> — if you already have a <code>baitbench</code> conda environment, BaitBench will detect it automatically or let you browse to it.</li>
</ul>
<p>After setup completes, BaitBench is ready to use.</p>
</div>

</div>

<script>
function qsTab(e, id) {
  document.querySelectorAll('.qs-tab').forEach(function(t){ t.classList.remove('active'); });
  document.querySelectorAll('.qs-panel').forEach(function(p){ p.classList.remove('active'); });
  e.target.classList.add('active');
  document.getElementById('qs-' + id).classList.add('active');
}
</script>
