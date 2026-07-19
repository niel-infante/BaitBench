<script lang="ts">
  import { onMount, onDestroy } from 'svelte';
  import { invoke } from '@tauri-apps/api/core';
  import { listen, type UnlistenFn } from '@tauri-apps/api/event';
  import { open as dialogOpen } from '@tauri-apps/plugin-dialog';
  import { load } from '@tauri-apps/plugin-store';
  import { get } from 'svelte/store';
  import { currentView, condaEnvPath } from '../stores';
  import type { CondaEnv, ValidationResult, SetupCheck, SetupProgress } from '../types';

  type WizardStep = 'scanning' | 'welcome' | 'existing' | 'installing' | 'creating' | 'done' | 'error';

  let step: WizardStep = 'scanning';
  let setupCheck: SetupCheck | null = null;
  let errorMsg = '';

  let progressMessage = '';
  let progressPercent: number | null = null;
  let progressLog: string[] = [];

  let installCli = true;
  let cliPath = '';
  let cliInstalling = false;

  let envs: CondaEnv[] = [];
  let selectedPath = '';
  let customPath = '';
  let validation: ValidationResult | null = null;
  let validating = false;
  let detecting = false;

  let condaExe = '';
  let cameFromRun = false;
  let unlisten: UnlistenFn | null = null;

  onMount(async () => {
    unlisten = await listen<SetupProgress>('setup-progress', (event) => {
      const p = event.payload;
      progressMessage = p.message;
      if (p.percent !== null && p.percent !== undefined) progressPercent = p.percent;
      progressLog = [...progressLog, p.message].slice(-200);
    });
    await runScan();
  });

  onDestroy(() => { if (unlisten) unlisten(); });

  async function runScan() {
    step = 'scanning';
    errorMsg = '';
    try {
      setupCheck = await invoke<SetupCheck>('check_setup');
      if (setupCheck.baitbench_env_valid) {
        if (get(condaEnvPath)) {
          // User explicitly navigated here from RunView to change the env
          cameFromRun = true;
          await openExistingFlow();
          return;
        }
        await saveEnvAndContinue(setupCheck.baitbench_env!);
        return;
      }
      condaExe = setupCheck.conda_executable ?? '';
      step = 'welcome';
    } catch (e) {
      errorMsg = String(e);
      step = 'error';
    }
  }

  async function startAutoSetup() {
    progressLog = [];
    progressPercent = null;
    progressMessage = '';
    if (!setupCheck?.conda_found) {
      step = 'installing';
      try {
        condaExe = await invoke<string>('install_conda');
      } catch (e) {
        errorMsg = String(e);
        step = 'error';
        return;
      }
    }
    step = 'creating';
    progressLog = [];
    try {
      const envPath = await invoke<string>('create_baitbench_env', { condaExecutable: condaExe });
      await saveEnvAndContinue(envPath, false);
      step = 'done';
    } catch (e) {
      errorMsg = String(e);
      step = 'error';
    }
  }

  async function finish() {
    if (installCli) {
      cliInstalling = true;
      try {
        cliPath = await invoke<string>('install_cli_to_path');
      } catch (e) {
        cliPath = '';
        console.warn('CLI install failed:', e);
      } finally {
        cliInstalling = false;
      }
    }
    currentView.set('run');
  }

  async function openExistingFlow() {
    step = 'existing';
    detecting = true;
    try {
      envs = await invoke<CondaEnv[]>('detect_conda_envs');
    } finally {
      detecting = false;
    }
  }

  async function validate(path: string) {
    if (!path) return;
    validating = true;
    validation = null;
    try {
      validation = await invoke<ValidationResult>('validate_conda_env', { path });
    } finally {
      validating = false;
    }
  }

  async function browseForEnv() {
    const dir = await dialogOpen({ directory: true, multiple: false });
    if (typeof dir === 'string') {
      customPath = dir;
      selectedPath = '';
      await validate(dir);
    }
  }

  async function onSelectChange() {
    if (selectedPath) { customPath = ''; await validate(selectedPath); }
  }

  async function saveExisting() {
    const path = selectedPath || customPath;
    if (!path || !validation?.valid) return;
    await saveEnvAndContinue(path);
  }

  async function saveEnvAndContinue(path: string, navigate = true) {
    const store = await load('settings.json', { autoSave: true });
    await store.set('conda_env_path', path);
    await store.save();
    condaEnvPath.set(path);
    if (navigate) currentView.set('run');
  }

  $: activePath = selectedPath || customPath;
</script>

<div class="setup-view">
  <div class="hero">
    <div class="logo-text">BaitBench</div>
    <p class="tagline">Probe capture efficiency simulation</p>
  </div>

  {#if step === 'scanning'}
    <div class="card">
      <div class="spinner-row"><span class="spinner"></span> Checking for existing environment…</div>
    </div>

  {:else if step === 'welcome'}
    <div class="card">
      <h2>Welcome to BaitBench</h2>
      <p class="hint">BaitBench needs a conda environment with <code>cd-hit-est</code> and other tools. Set one up automatically or point to an existing environment.</p>
      {#if !setupCheck?.conda_found}
        <div class="info-box">No conda installation detected. Automatic setup will download <strong>Miniforge3</strong> (~80 MB) to <code>~/miniforge3</code>, then create the <strong>baitbench</strong> environment.</div>
      {:else}
        <div class="info-box ok">Conda found. Automatic setup will create a new <strong>baitbench</strong> environment.</div>
      {/if}
      <button class="btn-primary wide" on:click={startAutoSetup}>Set up automatically</button>
      <button class="btn-link" on:click={openExistingFlow}>Use an existing environment instead</button>
    </div>

  {:else if step === 'installing'}
    <div class="card wide">
      <h2>Installing Miniforge3</h2>
      {#if progressPercent !== null}
        <div class="progress-bar-wrap"><div class="progress-bar" style="width: {Math.round(progressPercent * 100)}%"></div></div>
        <div class="pct-label">{Math.round(progressPercent * 100)}%</div>
      {:else}
        <div class="spinner-row"><span class="spinner"></span> Working…</div>
      {/if}
      {#if progressMessage}<p class="hint">{progressMessage}</p>{/if}
    </div>

  {:else if step === 'creating'}
    <div class="card wide">
      <h2>Creating baitbench environment</h2>
      <p class="hint">This can take 5–10 minutes. Please leave the app open.</p>
      <div class="spinner-row"><span class="spinner"></span> {progressMessage || 'Working…'}</div>
      <div class="log-box">
        {#each progressLog as line}<div class="log-line">{line}</div>{/each}
      </div>
    </div>

  {:else if step === 'done'}
    <div class="card">
      <div class="success-icon">✓</div>
      <h2>BaitBench is ready!</h2>
      <p class="hint">The baitbench environment was created successfully.</p>
      <label class="checkbox-row">
        <input type="checkbox" bind:checked={installCli} />
        <span>Also install <code>baitbench</code> CLI to PATH</span>
      </label>
      <p class="hint indent">Lets you run <code>baitbench</code> from any terminal.{#if !installCli}<span class="muted"> (You can do this later from the app menu.)</span>{/if}</p>
      {#if cliPath}<div class="info-box ok">CLI installed to <code>{cliPath}</code></div>{/if}
      <button class="btn-primary wide" on:click={finish} disabled={cliInstalling}>
        {#if cliInstalling}<span class="spinner small"></span> Installing CLI…{:else}Open BaitBench{/if}
      </button>
    </div>

  {:else if step === 'existing'}
    <div class="card">
      <h2>Select an existing environment</h2>
      <p class="hint">Choose a conda environment that already has <code>cd-hit-est</code>.</p>
      {#if detecting}
        <div class="spinner-row"><span class="spinner"></span> Scanning for conda environments…</div>
      {:else if envs.length > 0}
        <label class="field-label" for="env-select">Detected environments</label>
        <select id="env-select" bind:value={selectedPath} on:change={onSelectChange} class="select">
          <option value="">— choose an environment —</option>
          {#each envs as env}<option value={env.path}>{env.name}</option>{/each}
        </select>
      {:else}
        <p class="muted">No conda environments detected automatically.</p>
      {/if}
      <div class="divider-row">
        <div class="divider-line"></div><span class="divider-text">or browse manually</span><div class="divider-line"></div>
      </div>
      <div class="browse-row">
        <input type="text" bind:value={customPath} placeholder="/path/to/conda/env" class="path-input"
          on:change={() => { if (customPath) { selectedPath = ''; validate(customPath); } }} />
        <button class="btn-secondary" on:click={browseForEnv}>Browse…</button>
      </div>
      {#if validating}<div class="spinner-row"><span class="spinner"></span> Validating…</div>{/if}
      {#if validation && activePath}
        <div class="validation" class:ok={validation.valid} class:fail={!validation.valid}>
          {#if validation.valid}
            <div class="val-header ok-text">✓ Environment is valid</div>
          {:else}
            <div class="val-header fail-text">✗ Missing required tools</div>
            <ul class="val-list">{#each validation.missing as m}<li>{m}</li>{/each}</ul>
          {/if}
          {#if validation.warnings.length > 0}
            <div class="warnings">
              <div class="warn-title">Optional tools missing:</div>
              <ul class="val-list warn-list">{#each validation.warnings as w}<li>{w}</li>{/each}</ul>
            </div>
          {/if}
        </div>
      {/if}
      <div class="btn-row">
        {#if cameFromRun}
          <button class="btn-secondary" on:click={() => currentView.set('run')}>Cancel</button>
        {:else}
          <button class="btn-secondary" on:click={() => (step = 'welcome')}>← Back</button>
        {/if}
        <button class="btn-primary" disabled={!activePath || !validation?.valid} on:click={saveExisting}>Save &amp; Continue</button>
      </div>
    </div>

  {:else if step === 'error'}
    <div class="card wide">
      <div class="fail-icon">✗</div>
      <h2>Setup failed</h2>
      <div class="error-msg">{errorMsg}</div>
      {#if progressLog.length > 0}
        <div class="log-box">
          {#each progressLog as line}<div class="log-line">{line}</div>{/each}
        </div>
      {/if}
      <div class="btn-row">
        <button class="btn-secondary" on:click={openExistingFlow}>Use existing env</button>
        <button class="btn-primary" on:click={runScan}>Try again</button>
      </div>
    </div>
  {/if}
</div>

<style>
  .setup-view { display:flex; flex-direction:column; align-items:center; justify-content:center; min-height:100vh; padding:24px; gap:24px; }
  .hero { text-align:center; }
  .logo-text { font-size:2.2rem; font-weight:800; color:var(--color-primary); letter-spacing:-0.5px; }
  .tagline { color:var(--color-muted); margin:4px 0 0; font-size:0.95rem; }
  .card { background:var(--color-card); border:1px solid var(--color-border); border-radius:10px; padding:28px 32px; width:100%; max-width:520px; display:flex; flex-direction:column; gap:14px; }
  .card.wide { max-width:680px; }
  h2 { margin:0 0 2px; font-size:1.15rem; }
  .hint { color:var(--color-muted); font-size:0.85rem; margin:0; }
  .hint.indent { padding-left:22px; }
  .muted { color:var(--color-muted); font-size:0.85rem; }
  code { background:var(--color-chip); border-radius:3px; padding:1px 4px; font-size:0.8rem; }
  .info-box { background:#fffbeb; border:1px solid #fbd38d; border-radius:6px; padding:10px 12px; font-size:0.84rem; color:#744210; }
  .info-box.ok { background:#f0fff4; border-color:#9ae6b4; color:#276749; }
  .success-icon, .fail-icon { font-size:2.5rem; text-align:center; line-height:1; }
  .success-icon { color:#38a169; }
  .fail-icon { color:#e53e3e; }
  .progress-bar-wrap { background:var(--color-border); border-radius:4px; height:8px; overflow:hidden; }
  .progress-bar { height:100%; background:var(--color-primary); border-radius:4px; transition:width 0.3s ease; }
  .pct-label { font-size:0.8rem; color:var(--color-muted); text-align:right; margin-top:-8px; }
  .log-box { background:#1a1a1a; border-radius:6px; padding:10px 12px; font-family:monospace; font-size:0.75rem; color:#d4d4d4; height:180px; overflow-y:auto; display:flex; flex-direction:column; gap:2px; }
  .log-line { word-break:break-all; }
  .spinner-row { display:flex; align-items:center; gap:8px; font-size:0.85rem; color:var(--color-muted); }
  .spinner { width:14px; height:14px; border:2px solid var(--color-border); border-top-color:var(--color-primary); border-radius:50%; animation:spin 0.7s linear infinite; flex-shrink:0; }
  .spinner.small { width:11px; height:11px; }
  @keyframes spin { to { transform:rotate(360deg); } }
  .field-label { font-size:0.8rem; font-weight:600; color:var(--color-label); }
  .select, .path-input { width:100%; padding:7px 10px; font-size:0.88rem; border:1px solid var(--color-border); border-radius:5px; background:var(--color-input-bg); color:var(--color-text); box-sizing:border-box; }
  .browse-row { display:flex; gap:8px; }
  .path-input { flex:1; min-width:0; }
  .divider-row { display:flex; align-items:center; gap:8px; color:var(--color-muted); font-size:0.78rem; }
  .divider-line { flex:1; height:1px; background:var(--color-border); }
  .validation { border-radius:6px; padding:10px 12px; font-size:0.85rem; display:flex; flex-direction:column; gap:6px; }
  .validation.ok { background:#f0fff4; border:1px solid #9ae6b4; }
  .validation.fail { background:#fff5f5; border:1px solid #feb2b2; }
  .val-header { font-weight:600; }
  .ok-text { color:#276749; }
  .fail-text { color:#c53030; }
  .val-list { margin:4px 0 0; padding-left:18px; display:flex; flex-direction:column; gap:2px; color:#742a2a; }
  .warnings { margin-top:4px; }
  .warn-title { font-weight:600; color:#744210; font-size:0.82rem; }
  .warn-list { color:#744210; }
  .btn-primary { padding:9px 16px; background:var(--color-primary); color:#fff; border:none; border-radius:6px; font-size:0.9rem; font-weight:600; cursor:pointer; transition:opacity 0.15s; display:flex; align-items:center; justify-content:center; gap:6px; }
  .btn-primary.wide { width:100%; }
  .btn-primary:disabled { opacity:0.4; cursor:not-allowed; }
  .btn-primary:not(:disabled):hover { opacity:0.88; }
  .btn-secondary { padding:7px 12px; background:var(--color-btn); color:var(--color-text); border:1px solid var(--color-border); border-radius:5px; font-size:0.85rem; cursor:pointer; white-space:nowrap; }
  .btn-secondary:hover { background:var(--color-btn-hover); }
  .btn-link { background:none; border:none; color:var(--color-primary); font-size:0.85rem; cursor:pointer; padding:2px 0; text-decoration:underline; text-align:left; }
  .btn-row { display:flex; gap:8px; justify-content:flex-end; }
  .checkbox-row { display:flex; align-items:center; gap:8px; font-size:0.88rem; cursor:pointer; }
  .checkbox-row input { cursor:pointer; }
  .error-msg { background:#fff5f5; border:1px solid #feb2b2; border-radius:5px; padding:8px 10px; font-size:0.82rem; color:#c53030; word-break:break-word; }
</style>
