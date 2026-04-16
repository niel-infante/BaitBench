<script lang="ts">
  import { invoke } from '@tauri-apps/api/core';
  import { currentView, pipelineStatus, logLines, reportPath, condaEnvPath } from '../stores';
  import FilePicker from '../components/FilePicker.svelte';
  import AdvancedOptions from '../components/AdvancedOptions.svelte';
  import type { PipelineConfig, ToolDef } from '../types';

  // ── Tool definitions ─────────────────────────────────────────────────────
  const TOOLS: ToolDef[] = [
    { id: 'run',            category: 'Simulation',  label: 'Run Pipeline',         description: 'Full end-to-end pipeline (prepare → simulate → sequence → map → metrics → report)' },
    { id: 'coverage-curve', category: 'Simulation',  label: 'Coverage Curve',       description: 'Sweep across CT / capture-fraction values and plot coverage depth curves' },
    { id: 'build-probes',   category: 'Probes',      label: 'Build Probes',         description: 'Design a probe set from target sequences' },
    { id: 'assess-probes',  category: 'Probes',      label: 'Assess Probes',        description: 'Combined probe coverage + cross-reactivity assessment' },
    { id: 'xreact',         category: 'Analysis',    label: 'Cross-Reactivity',     description: 'Check probe cross-reactivity against genomes or other probes' },
    { id: 'panel-qc',       category: 'Analysis',    label: 'Panel QC',             description: 'Assess target panel discriminability between species' },
    { id: 'identify',       category: 'Analysis',    label: 'Identify Species',     description: 'Call species present/absent/ambiguous from detection patterns' },
  ];

  const CATEGORIES = [...new Set(TOOLS.map(t => t.category))];

  let selectedTool = TOOLS[0].id;
  let runError = '';
  let launching = false;

  $: tool = TOOLS.find(t => t.id === selectedTool)!;

  // ── Shared state ──────────────────────────────────────────────────────────
  let outdir = '';
  let threads = '4';
  let outputPrefix = '';
  let reportMode = 'full';
  let cleanup = false;

  // ── run ───────────────────────────────────────────────────────────────────
  let r_targets = '';
  let r_probes = '';
  let r_distractors = '';
  let r_sample = '';
  let r_sampleIsFile = true;
  let r_sampleInline = '';
  let r_distractorMode: 'fraction' | 'ct' = 'fraction';
  let r_distractorFraction = '0.9';
  let r_ct = '25';
  let r_ctBaseline = '20';
  let r_ctBaselineFraction = '0.01';
  let r_simulateMode = 'thermodynamic';
  let r_numFragments = '10000';
  let r_captureFraction = '0.5';
  let r_hybTemp = '70';
  let r_readLength = '120';
  let r_seed = '';
  let r_genomes = '';
  let r_sampleTargetMap = '';
  let r_groups = '';
  let r_hostFasta = '';
  let r_identify = false;
  let r_runName = '';
  // fragment params
  let r_fragLenMean = '175';
  let r_fragLenMin = '150';
  let r_fragLenMax = '200';
  let r_minimapPreset = 'sr';

  // ── build-probes ──────────────────────────────────────────────────────────
  let bp_targets = '';
  let bp_method = 'tile';
  let bp_probeLength = '120';
  let bp_step = '-60';
  let bp_catchStride = '60';
  let bp_catchMismatches = '5';
  let bp_syottiMismatches = '40';
  let bp_syottiSeedLen = '20';
  let bp_minGc = '0.20';
  let bp_maxGc = '0.80';
  let bp_dustThreshold = '2.0';
  let bp_skipAssess = false;
  let bp_genomes = '';

  // ── assess-probes ─────────────────────────────────────────────────────────
  let ap_targets = '';
  let ap_probes = '';
  let ap_genomes = '';
  let ap_threshold = '80';
  let ap_allIndividual = false;

  // ── xreact ───────────────────────────────────────────────────────────────
  let xr_probes = '';
  let xr_against = '';
  let xr_self = false;
  let xr_threshold = '80';

  // ── panel-qc ─────────────────────────────────────────────────────────────
  let pq_targets = '';
  let pq_sampleTargetMap = '';
  let pq_identityThreshold = '90';

  // ── identify ─────────────────────────────────────────────────────────────
  let id_detectedDetail = '';
  let id_sampleTargetMap = '';
  let id_targets = '';
  let id_identityThreshold = '90';
  let id_minUniqueTargets = '1';

  // ── coverage-curve ────────────────────────────────────────────────────────
  let cc_targets = '';
  let cc_probes = '';
  let cc_distractors = '';
  let cc_sample = '';
  let cc_ctValues = '20 25 30';
  let cc_captureFractionFixed = '0.5';
  let cc_numFragments = '10000';
  let cc_simulateMode = 'thermodynamic';

  // ── Build config and run ──────────────────────────────────────────────────
  function buildArgs(): Record<string, string> | null {
    const a: Record<string, string> = {};
    const add = (flag: string, val: string | number) => { const s = String(val); if (s !== '' && s !== 'undefined' && s !== 'null') a[flag] = s; };
    const flag = (f: string, on: boolean) => { if (on) a[f] = ''; };

    // Shared
    if (outdir) add('--outdir', outdir);
    if (outputPrefix) add('--output-prefix', outputPrefix);
    if (reportMode !== 'full') add('--report', reportMode);
    flag('--cleanup', cleanup);

    switch (selectedTool) {
      case 'run': {
        if (!r_targets || !r_probes || !r_distractors) {
          runError = 'Targets, probes, and at least one distractor file are required.';
          return null;
        }
        add('--targets', r_targets);
        add('--probes', r_probes);
        // Multi-value: split by tab into separate args representation
        add('--distractors', r_distractors); // tab-separated handled in Rust
        // sample
        const sampleVal = r_sampleIsFile ? r_sample : r_sampleInline;
        if (sampleVal) add('--sample', sampleVal);
        // distractor mode
        if (r_distractorMode === 'fraction') {
          add('--distractor-fraction', r_distractorFraction);
        } else {
          add('--ct', r_ct);
          if (r_ctBaseline !== '20') add('--ct-baseline', r_ctBaseline);
          if (r_ctBaselineFraction !== '0.01') add('--ct-baseline-fraction', r_ctBaselineFraction);
        }
        add('--simulate-mode', r_simulateMode);
        add('--num-fragments', r_numFragments);
        add('--capture-fraction', r_captureFraction);
        if (r_simulateMode === 'thermodynamic') add('--hybridization-temperature', r_hybTemp);
        add('--read-length', r_readLength);
        if (r_seed) add('--seed', r_seed);
        add('--threads', threads);
        // optional
        if (r_genomes) add('--genomes', r_genomes);
        if (r_sampleTargetMap) add('--sample-target-map', r_sampleTargetMap);
        if (r_groups) add('--groups', r_groups);
        if (r_hostFasta) add('--host-fasta', r_hostFasta);
        flag('--identify', r_identify);
        if (r_runName) add('--run-name', r_runName);
        // Advanced fragment params
        add('--fragment-length-mean', r_fragLenMean);
        add('--fragment-length-min', r_fragLenMin);
        add('--fragment-length-max', r_fragLenMax);
        add('--minimap-preset', r_minimapPreset);
        break;
      }

      case 'build-probes': {
        if (!bp_targets) { runError = 'Targets file is required.'; return null; }
        add('--targets', bp_targets);
        add('--method', bp_method);
        add('--probe-length', bp_probeLength);
        if (bp_method === 'tile') add('--step', bp_step);
        if (bp_method === 'catch-lite' || bp_method === 'catch') {
          add('--catch-probe-stride', bp_catchStride);
          add('--catch-mismatches', bp_catchMismatches);
        }
        if (bp_method === 'syotti-lite') {
          add('--syotti-mismatches', bp_syottiMismatches);
          add('--syotti-seed-len', bp_syottiSeedLen);
        }
        add('--min-gc', bp_minGc);
        add('--max-gc', bp_maxGc);
        add('--dust-threshold', bp_dustThreshold);
        if (bp_genomes) add('--genomes', bp_genomes);
        flag('--skip-assess', bp_skipAssess);
        add('--threads', threads);
        break;
      }

      case 'assess-probes': {
        if (!ap_targets || !ap_probes) { runError = 'Targets and probes files are required.'; return null; }
        add('--targets', ap_targets);
        add('--probes', ap_probes);
        if (ap_genomes) add('--genomes', ap_genomes);
        add('--threshold', ap_threshold);
        flag('--all-individual-targets', ap_allIndividual);
        add('--threads', threads);
        break;
      }

      case 'xreact': {
        if (!xr_probes) { runError = 'Probes file is required.'; return null; }
        add('--probes', xr_probes);
        if (xr_against) add('--against', xr_against);
        flag('--self', xr_self);
        add('--threshold', xr_threshold);
        if (!xr_against && !xr_self) {
          runError = 'Specify at least one of: genome FASTA (Against) or Self cross-reactivity check.';
          return null;
        }
        break;
      }

      case 'panel-qc': {
        if (!pq_targets || !pq_sampleTargetMap) { runError = 'Targets and sample-target-map are required.'; return null; }
        add('--targets', pq_targets);
        add('--sample-target-map', pq_sampleTargetMap);
        add('--identity-threshold', pq_identityThreshold);
        break;
      }

      case 'identify': {
        if (!id_detectedDetail || !id_sampleTargetMap) { runError = 'detected_detail.tsv and sample-target-map are required.'; return null; }
        add('--detected-detail', id_detectedDetail);
        add('--sample-target-map', id_sampleTargetMap);
        if (id_targets) add('--targets', id_targets);
        add('--identity-threshold', id_identityThreshold);
        add('--min-unique-targets', id_minUniqueTargets);
        break;
      }

      case 'coverage-curve': {
        if (!cc_targets || !cc_probes || !cc_distractors) { runError = 'Targets, probes, and distractors are required.'; return null; }
        add('--targets', cc_targets);
        add('--probes', cc_probes);
        add('--distractors', cc_distractors);
        if (cc_sample) add('--sample', cc_sample);
        // CT values as a space-separated list passed as individual args
        if (cc_ctValues) {
          const vals = cc_ctValues.trim().split(/\s+/);
          a['--ct-values'] = vals.join('\t');
        }
        add('--capture-fraction', cc_captureFractionFixed);
        add('--num-fragments', cc_numFragments);
        add('--simulate-mode', cc_simulateMode);
        add('--threads', threads);
        break;
      }
    }

    return a;
  }

  async function runTool() {
    runError = '';
    const args = buildArgs();
    if (!args) return;

    const config: PipelineConfig = {
      tool: selectedTool,
      args,
      conda_env: $condaEnvPath,
    };

    launching = true;
    logLines.set([]);
    reportPath.set(null);
    pipelineStatus.set('running');

    // Infer report path from outdir for the "run" tool
    if (outdir && selectedTool === 'run') {
      reportPath.set(`${outdir}/report.html`);
    }

    try {
      await invoke('run_pipeline', { config });
      currentView.set('log');
    } catch (e) {
      runError = String(e);
      pipelineStatus.set('idle');
    } finally {
      launching = false;
    }
  }

  function goToSetup() {
    currentView.set('setup');
  }

  $: isRunning = $pipelineStatus === 'running';
</script>

<div class="run-view">
  <!-- Header -->
  <div class="header">
    <span class="logo">BaitBench</span>
    <div class="spacer"></div>
    <button class="btn-ghost small" on:click={goToSetup}>⚙ Change Environment</button>
    {#if isRunning}
      <button class="btn-ghost small active-run" on:click={() => currentView.set('log')}>
        ● View Running Pipeline
      </button>
    {/if}
  </div>

  <div class="body">
    <!-- Tool selector sidebar -->
    <nav class="sidebar">
      {#each CATEGORIES as cat}
        <div class="cat-label">{cat}</div>
        {#each TOOLS.filter(t => t.category === cat) as t}
          <button
            class="tool-btn"
            class:active={selectedTool === t.id}
            on:click={() => { selectedTool = t.id; runError = ''; }}
          >{t.label}</button>
        {/each}
      {/each}
    </nav>

    <!-- Form area -->
    <div class="form-area">
      <div class="tool-header">
        <h2>{tool.label}</h2>
        <p class="tool-desc">{tool.description}</p>
      </div>

      <div class="form-scroll">
        <!-- ── run ─────────────────────────────────────── -->
        {#if selectedTool === 'run'}
          <section class="form-section">
            <h3>Inputs</h3>
            <FilePicker label="Targets FASTA" bind:value={r_targets} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <FilePicker label="Probes FASTA" bind:value={r_probes} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta'] }]} />
            <FilePicker label="Distractor FASTA(s)" bind:value={r_distractors} multiple required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
          </section>

          <section class="form-section">
            <h3>Sample Manifest</h3>
            <div class="toggle-row">
              <button class="toggle-btn" class:active={r_sampleIsFile}
                on:click={() => r_sampleIsFile = true}>File (TSV)</button>
              <button class="toggle-btn" class:active={!r_sampleIsFile}
                on:click={() => r_sampleIsFile = false}>Inline IDs</button>
            </div>
            {#if r_sampleIsFile}
              <FilePicker label="Sample manifest" bind:value={r_sample}
                filters={[{ name: 'TSV', extensions: ['tsv', 'txt'] }]} />
            {:else}
              <label class="field-label" for="sample-inline">
                Space-separated target IDs (optionally followed by weight)
              </label>
              <input id="sample-inline" class="text-input" type="text"
                bind:value={r_sampleInline}
                placeholder="target_1 target_2 3 target_3" />
            {/if}
          </section>

          <section class="form-section">
            <h3>Distractor Fraction</h3>
            <div class="toggle-row">
              <button class="toggle-btn" class:active={r_distractorMode === 'fraction'}
                on:click={() => r_distractorMode = 'fraction'}>Fraction</button>
              <button class="toggle-btn" class:active={r_distractorMode === 'ct'}
                on:click={() => r_distractorMode = 'ct'}>CT Value</button>
            </div>
            {#if r_distractorMode === 'fraction'}
              <div class="field-row">
                <label class="field-label" for="dfrac">Distractor fraction (0–1)</label>
                <input id="dfrac" class="text-input short" type="number" min="0" max="1" step="0.01"
                  bind:value={r_distractorFraction} />
              </div>
            {:else}
              <div class="field-group">
                <div class="field-row">
                  <label class="field-label" for="ct">CT value</label>
                  <input id="ct" class="text-input short" type="number" step="0.1"
                    bind:value={r_ct} />
                </div>
                <AdvancedOptions label="CT calibration">
                  <div class="field-row">
                    <label class="field-label" for="ctb">CT baseline</label>
                    <input id="ctb" class="text-input short" type="number" step="0.1"
                      bind:value={r_ctBaseline} />
                  </div>
                  <div class="field-row">
                    <label class="field-label" for="ctbf">CT baseline fraction</label>
                    <input id="ctbf" class="text-input short" type="number" step="0.001" min="0" max="1"
                      bind:value={r_ctBaselineFraction} />
                  </div>
                </AdvancedOptions>
              </div>
            {/if}
          </section>

          <section class="form-section">
            <h3>Simulation</h3>
            <div class="field-row">
              <label class="field-label" for="simmode">Simulate mode</label>
              <select id="simmode" class="select-input" bind:value={r_simulateMode}>
                <option value="thermodynamic">Thermodynamic (TNN)</option>
                <option value="simple">Simple</option>
              </select>
            </div>
            {#if r_simulateMode === 'thermodynamic'}
              <div class="field-row">
                <label class="field-label" for="hybtemp">Hybridization temperature (°C)</label>
                <input id="hybtemp" class="text-input short" type="number" step="1"
                  bind:value={r_hybTemp} />
              </div>
            {/if}
            <div class="field-row">
              <label class="field-label" for="nfrags">Number of fragments</label>
              <input id="nfrags" class="text-input short" type="number" min="100"
                bind:value={r_numFragments} />
            </div>
            <div class="field-row">
              <label class="field-label" for="capfrac">Capture fraction</label>
              <input id="capfrac" class="text-input short" type="number" min="0" max="1" step="0.01"
                bind:value={r_captureFraction} />
            </div>
            <div class="field-row">
              <label class="field-label" for="readlen">Read length (bp)</label>
              <input id="readlen" class="text-input short" type="number" min="1"
                bind:value={r_readLength} />
            </div>
          </section>

          <section class="form-section">
            <h3>Output</h3>
            <FilePicker label="Output directory" bind:value={outdir} directory required />
            <div class="field-row">
              <label class="field-label" for="report-mode">Report</label>
              <select id="report-mode" class="select-input" bind:value={reportMode}>
                <option value="full">Full HTML</option>
                <option value="rmd">RMarkdown only</option>
                <option value="none">None</option>
              </select>
            </div>
          </section>

          <AdvancedOptions label="Advanced Options">
            <div class="field-row">
              <label class="field-label" for="threads">Threads</label>
              <input id="threads" class="text-input short" type="number" min="1"
                bind:value={threads} />
            </div>
            <div class="field-row">
              <label class="field-label" for="seed">Random seed (blank = random)</label>
              <input id="seed" class="text-input short" type="text" bind:value={r_seed}
                placeholder="e.g. 42" />
            </div>
            <div class="field-row">
              <label class="field-label" for="runname">Run name</label>
              <input id="runname" class="text-input" type="text" bind:value={r_runName} />
            </div>
            <div class="field-row">
              <label class="field-label" for="outprefix">Output prefix</label>
              <input id="outprefix" class="text-input" type="text" bind:value={outputPrefix} />
            </div>
            <div class="field-row">
              <label class="field-label" for="fraglenmean">Fragment length mean (bp)</label>
              <input id="fraglenmean" class="text-input short" type="number" bind:value={r_fragLenMean} />
            </div>
            <div class="field-row">
              <label class="field-label" for="fraglenmin">Fragment length min (bp)</label>
              <input id="fraglenmin" class="text-input short" type="number" bind:value={r_fragLenMin} />
            </div>
            <div class="field-row">
              <label class="field-label" for="fraglenmax">Fragment length max (bp)</label>
              <input id="fraglenmax" class="text-input short" type="number" bind:value={r_fragLenMax} />
            </div>
            <div class="field-row">
              <label class="field-label" for="minimap">Minimap2 preset</label>
              <input id="minimap" class="text-input short" type="text" bind:value={r_minimapPreset} />
            </div>
            <label class="check-label">
              <input type="checkbox" bind:checked={cleanup} />
              Clean up intermediate files after run
            </label>
            <h4 class="subsection">Genome Mode (optional)</h4>
            <FilePicker label="Genomes FASTA" bind:value={r_genomes}
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <FilePicker label="Sample-target map" bind:value={r_sampleTargetMap}
              filters={[{ name: 'TSV', extensions: ['tsv', 'txt'] }]} />
            <h4 class="subsection">Grouping (optional)</h4>
            <FilePicker label="Target groups" bind:value={r_groups}
              filters={[{ name: 'TSV', extensions: ['tsv', 'txt'] }]} />
            <h4 class="subsection">Host Filtering (optional)</h4>
            <FilePicker label="Host FASTA" bind:value={r_hostFasta}
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <h4 class="subsection">Species Identification</h4>
            <label class="check-label">
              <input type="checkbox" bind:checked={r_identify} />
              Run species-level identification after metrics (genome mode)
            </label>
          </AdvancedOptions>

        <!-- ── build-probes ──────────────────────────── -->
        {:else if selectedTool === 'build-probes'}
          <section class="form-section">
            <h3>Inputs</h3>
            <FilePicker label="Targets FASTA" bind:value={bp_targets} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
          </section>
          <section class="form-section">
            <h3>Design Method</h3>
            <div class="field-row">
              <label class="field-label" for="bp-method">Method</label>
              <select id="bp-method" class="select-input" bind:value={bp_method}>
                <option value="tile">Tile (sliding window)</option>
                <option value="catch-lite">CATCH-lite (native Rust)</option>
                <option value="syotti-lite">Syotti-lite (native Rust)</option>
                <option value="catch">CATCH (external, requires conda pkg)</option>
              </select>
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-probelen">Probe length (bp)</label>
              <input id="bp-probelen" class="text-input short" type="number" min="1"
                bind:value={bp_probeLength} />
            </div>
            {#if bp_method === 'tile'}
              <div class="field-row">
                <label class="field-label" for="bp-step">Step (negative = overlap)</label>
                <input id="bp-step" class="text-input short" type="number"
                  bind:value={bp_step} />
              </div>
            {:else if bp_method === 'catch-lite' || bp_method === 'catch'}
              <div class="field-row">
                <label class="field-label" for="bp-catchstride">CATCH probe stride</label>
                <input id="bp-catchstride" class="text-input short" type="number"
                  bind:value={bp_catchStride} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-catchmm">CATCH mismatches</label>
                <input id="bp-catchmm" class="text-input short" type="number"
                  bind:value={bp_catchMismatches} />
              </div>
            {:else if bp_method === 'syotti-lite'}
              <div class="field-row">
                <label class="field-label" for="bp-syottimm">Syotti mismatches</label>
                <input id="bp-syottimm" class="text-input short" type="number"
                  bind:value={bp_syottiMismatches} />
              </div>
              <div class="field-row">
                <label class="field-label" for="bp-syottiseed">Syotti seed length</label>
                <input id="bp-syottiseed" class="text-input short" type="number"
                  bind:value={bp_syottiSeedLen} />
              </div>
            {/if}
          </section>
          <section class="form-section">
            <h3>Output</h3>
            <FilePicker label="Output directory" bind:value={outdir} directory required />
          </section>
          <AdvancedOptions label="Filtering & Assessment">
            <div class="field-row">
              <label class="field-label" for="bp-mingc">Min GC</label>
              <input id="bp-mingc" class="text-input short" type="number" step="0.01" min="0" max="1"
                bind:value={bp_minGc} />
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-maxgc">Max GC</label>
              <input id="bp-maxgc" class="text-input short" type="number" step="0.01" min="0" max="1"
                bind:value={bp_maxGc} />
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-dust">sDUST threshold</label>
              <input id="bp-dust" class="text-input short" type="number" step="0.1"
                bind:value={bp_dustThreshold} />
            </div>
            <div class="field-row">
              <label class="field-label" for="bp-threads">Threads</label>
              <input id="bp-threads" class="text-input short" type="number" min="1"
                bind:value={threads} />
            </div>
            <FilePicker label="Genomes FASTA (for cross-reactivity check)"
              bind:value={bp_genomes}
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <label class="check-label">
              <input type="checkbox" bind:checked={bp_skipAssess} />
              Skip assessment step
            </label>
          </AdvancedOptions>

        <!-- ── assess-probes ────────────────────────── -->
        {:else if selectedTool === 'assess-probes'}
          <section class="form-section">
            <h3>Inputs</h3>
            <FilePicker label="Targets FASTA" bind:value={ap_targets} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <FilePicker label="Probes FASTA" bind:value={ap_probes} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta'] }]} />
            <FilePicker label="Genomes FASTA (optional)" bind:value={ap_genomes}
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
          </section>
          <section class="form-section">
            <h3>Output</h3>
            <FilePicker label="Output directory" bind:value={outdir} directory required />
          </section>
          <AdvancedOptions label="Options">
            <div class="field-row">
              <label class="field-label" for="ap-threshold">Homology threshold (%)</label>
              <input id="ap-threshold" class="text-input short" type="number" min="0" max="100"
                bind:value={ap_threshold} />
            </div>
            <label class="check-label">
              <input type="checkbox" bind:checked={ap_allIndividual} />
              Compute per-target coverage individually (slower but more detailed)
            </label>
            <div class="field-row">
              <label class="field-label" for="ap-threads">Threads</label>
              <input id="ap-threads" class="text-input short" type="number" min="1"
                bind:value={threads} />
            </div>
          </AdvancedOptions>

        <!-- ── xreact ───────────────────────────────── -->
        {:else if selectedTool === 'xreact'}
          <section class="form-section">
            <h3>Inputs</h3>
            <FilePicker label="Probes FASTA" bind:value={xr_probes} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta'] }]} />
            <FilePicker label="Against FASTA (genomes to check against)"
              bind:value={xr_against}
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <label class="check-label">
              <input type="checkbox" bind:checked={xr_self} />
              Probe-vs-probe self cross-reactivity
            </label>
          </section>
          <section class="form-section">
            <h3>Output</h3>
            <FilePicker label="Output directory" bind:value={outdir} directory required />
          </section>
          <AdvancedOptions label="Options">
            <div class="field-row">
              <label class="field-label" for="xr-threshold">Homology threshold (%)</label>
              <input id="xr-threshold" class="text-input short" type="number" min="0" max="100"
                bind:value={xr_threshold} />
            </div>
          </AdvancedOptions>

        <!-- ── panel-qc ─────────────────────────────── -->
        {:else if selectedTool === 'panel-qc'}
          <section class="form-section">
            <h3>Inputs</h3>
            <FilePicker label="Targets FASTA" bind:value={pq_targets} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <FilePicker label="Sample-target map (TSV)" bind:value={pq_sampleTargetMap} required
              filters={[{ name: 'TSV', extensions: ['tsv', 'txt'] }]} />
          </section>
          <section class="form-section">
            <h3>Output</h3>
            <FilePicker label="Output directory" bind:value={outdir} directory required />
          </section>
          <AdvancedOptions label="Options">
            <div class="field-row">
              <label class="field-label" for="pq-ident">Identity threshold (%)</label>
              <input id="pq-ident" class="text-input short" type="number" min="0" max="100"
                bind:value={pq_identityThreshold} />
            </div>
          </AdvancedOptions>

        <!-- ── identify ─────────────────────────────── -->
        {:else if selectedTool === 'identify'}
          <section class="form-section">
            <h3>Inputs</h3>
            <FilePicker label="detected_detail.tsv" bind:value={id_detectedDetail} required
              filters={[{ name: 'TSV', extensions: ['tsv'] }]} />
            <FilePicker label="Sample-target map (TSV)" bind:value={id_sampleTargetMap} required
              filters={[{ name: 'TSV', extensions: ['tsv', 'txt'] }]} />
            <FilePicker label="Targets FASTA (optional, for similarity)" bind:value={id_targets}
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
          </section>
          <section class="form-section">
            <h3>Output</h3>
            <FilePicker label="Output directory" bind:value={outdir} directory required />
          </section>
          <AdvancedOptions label="Options">
            <div class="field-row">
              <label class="field-label" for="id-ident">Identity threshold (%)</label>
              <input id="id-ident" class="text-input short" type="number" min="0" max="100"
                bind:value={id_identityThreshold} />
            </div>
            <div class="field-row">
              <label class="field-label" for="id-minuniq">Min unique targets for call</label>
              <input id="id-minuniq" class="text-input short" type="number" min="1"
                bind:value={id_minUniqueTargets} />
            </div>
          </AdvancedOptions>

        <!-- ── coverage-curve ────────────────────────── -->
        {:else if selectedTool === 'coverage-curve'}
          <section class="form-section">
            <h3>Inputs</h3>
            <FilePicker label="Targets FASTA" bind:value={cc_targets} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <FilePicker label="Probes FASTA" bind:value={cc_probes} required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta'] }]} />
            <FilePicker label="Distractor FASTA(s)" bind:value={cc_distractors} multiple required
              filters={[{ name: 'FASTA', extensions: ['fa', 'fasta', 'fna'] }]} />
            <FilePicker label="Sample manifest (optional)" bind:value={cc_sample}
              filters={[{ name: 'TSV', extensions: ['tsv', 'txt'] }]} />
          </section>
          <section class="form-section">
            <h3>Parameter Sweep</h3>
            <div class="field-row">
              <label class="field-label" for="cc-ctvals">CT values (space-separated)</label>
              <input id="cc-ctvals" class="text-input" type="text"
                bind:value={cc_ctValues} placeholder="20 25 30" />
            </div>
            <div class="field-row">
              <label class="field-label" for="cc-capfrac">Capture fraction (fixed)</label>
              <input id="cc-capfrac" class="text-input short" type="number" min="0" max="1" step="0.01"
                bind:value={cc_captureFractionFixed} />
            </div>
            <div class="field-row">
              <label class="field-label" for="cc-nfrags">Number of fragments</label>
              <input id="cc-nfrags" class="text-input short" type="number" min="100"
                bind:value={cc_numFragments} />
            </div>
            <div class="field-row">
              <label class="field-label" for="cc-simmode">Simulate mode</label>
              <select id="cc-simmode" class="select-input" bind:value={cc_simulateMode}>
                <option value="thermodynamic">Thermodynamic</option>
                <option value="simple">Simple</option>
              </select>
            </div>
          </section>
          <section class="form-section">
            <h3>Output</h3>
            <FilePicker label="Output directory" bind:value={outdir} directory required />
          </section>
          <AdvancedOptions label="Advanced">
            <div class="field-row">
              <label class="field-label" for="cc-threads">Threads</label>
              <input id="cc-threads" class="text-input short" type="number" min="1"
                bind:value={threads} />
            </div>
          </AdvancedOptions>
        {/if}

        <!-- Error -->
        {#if runError}
          <div class="run-error">{runError}</div>
        {/if}

        <!-- Run button -->
        <div class="run-row">
          <button
            class="btn-run"
            disabled={isRunning || launching}
            on:click={runTool}
          >
            {#if launching}
              <span class="spinner-sm"></span> Starting…
            {:else if isRunning}
              Running…
            {:else}
              ▶ Run {tool.label}
            {/if}
          </button>
          {#if isRunning}
            <button class="btn-ghost" on:click={() => currentView.set('log')}>
              View log →
            </button>
          {/if}
        </div>
      </div>
    </div>
  </div>
</div>

<style>
  .run-view {
    display: flex;
    flex-direction: column;
    height: 100vh;
    overflow: hidden;
  }
  .header {
    display: flex;
    align-items: center;
    gap: 10px;
    padding: 8px 16px;
    border-bottom: 1px solid var(--color-border);
    background: var(--color-card);
    flex-shrink: 0;
  }
  .logo {
    font-size: 1.05rem;
    font-weight: 800;
    color: var(--color-primary);
  }
  .spacer { flex: 1; }
  .body {
    display: flex;
    flex: 1;
    overflow: hidden;
  }
  /* Sidebar */
  .sidebar {
    width: 180px;
    flex-shrink: 0;
    border-right: 1px solid var(--color-border);
    background: var(--color-card);
    padding: 10px 0;
    overflow-y: auto;
    display: flex;
    flex-direction: column;
    gap: 1px;
  }
  .cat-label {
    font-size: 0.7rem;
    font-weight: 700;
    text-transform: uppercase;
    color: var(--color-muted);
    letter-spacing: 0.04em;
    padding: 8px 12px 3px;
  }
  .tool-btn {
    width: 100%;
    text-align: left;
    padding: 7px 14px;
    font-size: 0.85rem;
    border: none;
    background: transparent;
    cursor: pointer;
    color: var(--color-text);
    border-radius: 0;
  }
  .tool-btn:hover { background: var(--color-btn-hover); }
  .tool-btn.active {
    background: var(--color-primary-light);
    color: var(--color-primary);
    font-weight: 600;
  }
  /* Form area */
  .form-area {
    flex: 1;
    display: flex;
    flex-direction: column;
    overflow: hidden;
  }
  .tool-header {
    padding: 14px 20px 10px;
    border-bottom: 1px solid var(--color-border);
    flex-shrink: 0;
  }
  .tool-header h2 { margin: 0 0 3px; font-size: 1.05rem; }
  .tool-desc { margin: 0; font-size: 0.82rem; color: var(--color-muted); }
  .form-scroll {
    flex: 1;
    overflow-y: auto;
    padding: 16px 20px 32px;
    display: flex;
    flex-direction: column;
    gap: 16px;
  }
  .form-section {
    display: flex;
    flex-direction: column;
    gap: 10px;
  }
  .form-section h3 {
    margin: 0 0 2px;
    font-size: 0.82rem;
    font-weight: 700;
    text-transform: uppercase;
    letter-spacing: 0.04em;
    color: var(--color-muted);
  }
  .field-row {
    display: flex;
    align-items: center;
    gap: 10px;
  }
  .field-row .field-label {
    flex: 1;
    font-size: 0.85rem;
    font-weight: 600;
    color: var(--color-label);
    min-width: 160px;
  }
  .field-group { display: flex; flex-direction: column; gap: 8px; }
  .field-label { font-size: 0.85rem; font-weight: 600; color: var(--color-label); }
  .text-input, .select-input {
    padding: 5px 8px;
    font-size: 0.85rem;
    border: 1px solid var(--color-border);
    border-radius: 4px;
    background: var(--color-input-bg);
    color: var(--color-text);
  }
  .text-input { flex: 1; }
  .text-input.short { width: 90px; flex: none; }
  .select-input { flex: 1; }
  .check-label {
    display: flex;
    align-items: center;
    gap: 8px;
    font-size: 0.85rem;
    color: var(--color-text);
    cursor: pointer;
  }
  .toggle-row {
    display: flex;
    border: 1px solid var(--color-border);
    border-radius: 5px;
    overflow: hidden;
    width: fit-content;
  }
  .toggle-btn {
    padding: 4px 12px;
    font-size: 0.82rem;
    border: none;
    background: var(--color-input-bg);
    cursor: pointer;
    color: var(--color-muted);
  }
  .toggle-btn.active {
    background: var(--color-primary);
    color: #fff;
    font-weight: 600;
  }
  .subsection {
    margin: 8px 0 0;
    font-size: 0.78rem;
    font-weight: 700;
    text-transform: uppercase;
    letter-spacing: 0.04em;
    color: var(--color-muted);
  }
  .run-error {
    background: #fff5f5;
    border: 1px solid #feb2b2;
    border-radius: 5px;
    padding: 8px 12px;
    font-size: 0.83rem;
    color: #c53030;
  }
  .run-row {
    display: flex;
    align-items: center;
    gap: 12px;
    margin-top: 4px;
  }
  .btn-run {
    padding: 10px 22px;
    background: var(--color-primary);
    color: #fff;
    border: none;
    border-radius: 6px;
    font-size: 0.92rem;
    font-weight: 700;
    cursor: pointer;
    display: flex;
    align-items: center;
    gap: 7px;
  }
  .btn-run:disabled { opacity: 0.5; cursor: not-allowed; }
  .btn-run:not(:disabled):hover { opacity: 0.88; }
  .spinner-sm {
    width: 12px; height: 12px;
    border: 2px solid rgba(255,255,255,0.4);
    border-top-color: #fff;
    border-radius: 50%;
    animation: spin 0.7s linear infinite;
  }
  @keyframes spin { to { transform: rotate(360deg); } }
  .btn-ghost {
    padding: 5px 10px;
    background: transparent;
    border: 1px solid var(--color-border);
    border-radius: 5px;
    font-size: 0.82rem;
    cursor: pointer;
    color: var(--color-text);
  }
  .btn-ghost:hover { background: var(--color-btn-hover); }
  .btn-ghost.small { font-size: 0.8rem; padding: 4px 8px; }
  .active-run { color: #2b6cb0; border-color: #bee3f8; }
</style>
