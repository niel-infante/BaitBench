<script lang="ts">
  import { onMount } from 'svelte';
  import { invoke } from '@tauri-apps/api/core';
  import { open as dialogOpen } from '@tauri-apps/plugin-dialog';
  import { load } from '@tauri-apps/plugin-store';
  import { currentView, condaEnvPath } from '../stores';
  import type { CondaEnv, ValidationResult } from '../types';

  let envs: CondaEnv[] = [];
  let selectedPath: string = '';
  let customPath: string = '';
  let validation: ValidationResult | null = null;
  let validating = false;
  let detecting = true;
  let error = '';

  onMount(async () => {
    detecting = true;
    try {
      envs = await invoke<CondaEnv[]>('detect_conda_envs');
    } catch (e) {
      error = String(e);
    } finally {
      detecting = false;
    }
  });

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
      selectedPath = dir;
      await validate(dir);
    }
  }

  async function onSelectChange() {
    if (selectedPath) {
      customPath = '';
      await validate(selectedPath);
    }
  }

  async function saveAndContinue() {
    const path = selectedPath || customPath;
    if (!path || !validation?.valid) return;

    const store = await load('settings.json', { autoSave: true });
    await store.set('conda_env_path', path);
    await store.save();

    condaEnvPath.set(path);
    currentView.set('run');
  }

  $: activePath = selectedPath || customPath;
</script>

<div class="setup-view">
  <div class="hero">
    <div class="logo-text">BaitBench</div>
    <p class="tagline">Probe capture efficiency simulation</p>
  </div>

  <div class="card">
    <h2>Environment Setup</h2>
    <p class="hint">
      BaitBench requires a conda environment with <code>minimap2</code> and
      <code>cd-hit-est</code>. Select the environment below.
    </p>

    {#if detecting}
      <div class="spinner-row"><span class="spinner"></span> Scanning for conda environments…</div>
    {:else if envs.length > 0}
      <label class="field-label" for="env-select">Detected environments</label>
      <select
        id="env-select"
        bind:value={selectedPath}
        on:change={onSelectChange}
        class="select"
      >
        <option value="">— choose an environment —</option>
        {#each envs as env}
          <option value={env.path}>{env.name}</option>
        {/each}
      </select>
    {:else}
      <p class="muted">No conda environments detected automatically.</p>
    {/if}

    <div class="divider-row">
      <div class="divider-line"></div>
      <span class="divider-text">or browse manually</span>
      <div class="divider-line"></div>
    </div>

    <div class="browse-row">
      <input
        type="text"
        bind:value={customPath}
        placeholder="/path/to/conda/env"
        class="path-input"
        on:change={() => { if (customPath) { selectedPath = ''; validate(customPath); } }}
      />
      <button class="btn-secondary" on:click={browseForEnv}>Browse…</button>
    </div>

    {#if validating}
      <div class="spinner-row"><span class="spinner"></span> Validating…</div>
    {/if}

    {#if validation && activePath}
      <div class="validation" class:ok={validation.valid} class:fail={!validation.valid}>
        {#if validation.valid}
          <div class="val-header ok-text">✓ Environment is valid</div>
        {:else}
          <div class="val-header fail-text">✗ Missing required tools</div>
          <ul class="val-list">
            {#each validation.missing as m}
              <li>{m}</li>
            {/each}
          </ul>
        {/if}
        {#if validation.warnings.length > 0}
          <div class="warnings">
            <div class="warn-title">Warnings (optional tools missing):</div>
            <ul class="val-list warn-list">
              {#each validation.warnings as w}
                <li>{w}</li>
              {/each}
            </ul>
          </div>
        {/if}
      </div>
    {/if}

    {#if error}
      <div class="error-msg">{error}</div>
    {/if}

    <button
      class="btn-primary save-btn"
      disabled={!activePath || !validation?.valid}
      on:click={saveAndContinue}
    >
      Save &amp; Continue
    </button>
  </div>
</div>

<style>
  .setup-view {
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    min-height: 100vh;
    padding: 24px;
    gap: 24px;
  }
  .hero {
    text-align: center;
  }
  .logo-text {
    font-size: 2.2rem;
    font-weight: 800;
    color: var(--color-primary);
    letter-spacing: -0.5px;
  }
  .tagline {
    color: var(--color-muted);
    margin: 4px 0 0;
    font-size: 0.95rem;
  }
  .card {
    background: var(--color-card);
    border: 1px solid var(--color-border);
    border-radius: 10px;
    padding: 28px 32px;
    width: 100%;
    max-width: 520px;
    display: flex;
    flex-direction: column;
    gap: 14px;
  }
  h2 {
    margin: 0 0 2px;
    font-size: 1.15rem;
  }
  .hint {
    color: var(--color-muted);
    font-size: 0.85rem;
    margin: 0;
  }
  code {
    background: var(--color-chip);
    border-radius: 3px;
    padding: 1px 4px;
    font-size: 0.8rem;
  }
  .field-label {
    font-size: 0.8rem;
    font-weight: 600;
    color: var(--color-label);
  }
  .select,
  .path-input {
    width: 100%;
    padding: 7px 10px;
    font-size: 0.88rem;
    border: 1px solid var(--color-border);
    border-radius: 5px;
    background: var(--color-input-bg);
    color: var(--color-text);
    box-sizing: border-box;
  }
  .browse-row {
    display: flex;
    gap: 8px;
  }
  .path-input {
    flex: 1;
    min-width: 0;
  }
  .divider-row {
    display: flex;
    align-items: center;
    gap: 8px;
    color: var(--color-muted);
    font-size: 0.78rem;
  }
  .divider-line {
    flex: 1;
    height: 1px;
    background: var(--color-border);
  }
  .spinner-row {
    display: flex;
    align-items: center;
    gap: 8px;
    font-size: 0.85rem;
    color: var(--color-muted);
  }
  .spinner {
    width: 14px;
    height: 14px;
    border: 2px solid var(--color-border);
    border-top-color: var(--color-primary);
    border-radius: 50%;
    animation: spin 0.7s linear infinite;
    flex-shrink: 0;
  }
  @keyframes spin { to { transform: rotate(360deg); } }
  .validation {
    border-radius: 6px;
    padding: 10px 12px;
    font-size: 0.85rem;
    display: flex;
    flex-direction: column;
    gap: 6px;
  }
  .validation.ok {
    background: #f0fff4;
    border: 1px solid #9ae6b4;
  }
  .validation.fail {
    background: #fff5f5;
    border: 1px solid #feb2b2;
  }
  .val-header { font-weight: 600; }
  .ok-text { color: #276749; }
  .fail-text { color: #c53030; }
  .val-list {
    margin: 4px 0 0;
    padding-left: 18px;
    display: flex;
    flex-direction: column;
    gap: 2px;
    color: #742a2a;
  }
  .warnings { margin-top: 4px; }
  .warn-title { font-weight: 600; color: #744210; font-size: 0.82rem; }
  .warn-list { color: #744210; }
  .error-msg {
    background: #fff5f5;
    border: 1px solid #feb2b2;
    border-radius: 5px;
    padding: 8px 10px;
    font-size: 0.82rem;
    color: #c53030;
  }
  .btn-primary {
    padding: 9px 16px;
    background: var(--color-primary);
    color: #fff;
    border: none;
    border-radius: 6px;
    font-size: 0.9rem;
    font-weight: 600;
    cursor: pointer;
    transition: opacity 0.15s;
  }
  .btn-primary:disabled { opacity: 0.4; cursor: not-allowed; }
  .btn-primary:not(:disabled):hover { opacity: 0.88; }
  .btn-secondary {
    padding: 7px 12px;
    background: var(--color-btn);
    color: var(--color-text);
    border: 1px solid var(--color-border);
    border-radius: 5px;
    font-size: 0.85rem;
    cursor: pointer;
    white-space: nowrap;
  }
  .btn-secondary:hover { background: var(--color-btn-hover); }
  .save-btn { width: 100%; margin-top: 4px; }
  .muted { color: var(--color-muted); font-size: 0.85rem; margin: 0; }
</style>
