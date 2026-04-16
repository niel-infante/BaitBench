<script lang="ts">
  import { onMount, onDestroy, tick } from 'svelte';
  import { invoke } from '@tauri-apps/api/core';
  import { listen } from '@tauri-apps/api/event';
  import { currentView, pipelineStatus, logLines, reportPath } from '../stores';
  import type { LogLine } from '../types';

  let logEl: HTMLElement;
  let autoScroll = true;
  let unlisten: (() => void)[] = [];

  // Match any absolute path ending in .html in a log line.
  const HTML_PATH_RE = /([\/~][^\s'"]+\.html)/;

  function extractHtmlPath(text: string): string | null {
    const m = text.match(HTML_PATH_RE);
    return m ? m[1] : null;
  }

  onMount(async () => {
    const u1 = await listen<LogLine>('pipeline-log', (event) => {
      logLines.update((lines) => [...lines, event.payload]);
      // Pick up the report path from log output as soon as it appears.
      const found = extractHtmlPath(event.payload.text);
      if (found) reportPath.set(found);
      if (autoScroll) scrollToBottom();
    });

    const u2 = await listen('pipeline-complete', async () => {
      pipelineStatus.set('complete');
      // Auto-open the report if we found one in the log output.
      const path = $reportPath ?? null;
      // Use a small delay to let the file system finish writing.
      if (path) {
        await new Promise(r => setTimeout(r, 800));
        await invoke('open_report', { path });
      }
    });

    const u3 = await listen('pipeline-error', () => {
      pipelineStatus.set('failed');
    });

    unlisten = [u1, u2, u3];

    // Sync status on mount in case we navigated here after the run already finished.
    const status = await invoke<string>('get_pipeline_status');
    pipelineStatus.set(status as any);
  });

  onDestroy(() => {
    unlisten.forEach((u) => u());
  });

  function scrollToBottom() {
    tick().then(() => {
      if (logEl) logEl.scrollTop = logEl.scrollHeight;
    });
  }

  function onScroll() {
    if (!logEl) return;
    const atBottom = logEl.scrollHeight - logEl.scrollTop - logEl.clientHeight < 40;
    autoScroll = atBottom;
  }

  async function cancel() {
    await invoke('cancel_pipeline');
    pipelineStatus.set('failed');
  }

  function lineClass(line: LogLine): string {
    const t = line.text;
    if (line.stream === 'stderr') {
      if (/^error/i.test(t) || /\[error\]/i.test(t)) return 'line-error';
      if (/^warn/i.test(t) || /\[warn\]/i.test(t)) return 'line-warn';
      return 'line-stderr';
    }
    return 'line-stdout';
  }

  function goBack() {
    if ($pipelineStatus === 'running') {
      if (!confirm('A pipeline is still running. Go back anyway? (it will keep running in the background)')) return;
    }
    currentView.set('run');
  }

  async function viewReport() {
    if (!$reportPath) return;
    await invoke('open_report', { path: $reportPath });
  }

  function clearLog() {
    logLines.set([]);
  }
</script>

<div class="log-view">
  <div class="toolbar">
    <button class="btn-ghost back-btn" on:click={goBack}>← Back to Run</button>
    <div class="status-pill" class:running={$pipelineStatus === 'running'}
         class:complete={$pipelineStatus === 'complete'}
         class:failed={$pipelineStatus === 'failed'}
         class:idle={$pipelineStatus === 'idle'}>
      {#if $pipelineStatus === 'running'}<span class="spinner"></span>{/if}
      {$pipelineStatus.charAt(0).toUpperCase() + $pipelineStatus.slice(1)}
    </div>
    <div class="spacer"></div>
    {#if $pipelineStatus === 'running'}
      <button class="btn-danger" on:click={cancel}>Cancel</button>
    {/if}
    {#if $pipelineStatus === 'complete' && $reportPath}
      <button class="btn-primary" on:click={viewReport}>View Report</button>
    {/if}
    <button class="btn-ghost" on:click={clearLog}>Clear</button>
  </div>

  <div
    class="log-panel"
    bind:this={logEl}
    on:scroll={onScroll}
  >
    {#if $logLines.length === 0}
      <div class="empty">No output yet…</div>
    {:else}
      {#each $logLines as line}
        <div class="log-line {lineClass(line)}">{line.text}</div>
      {/each}
    {/if}
  </div>

  {#if !autoScroll}
    <button class="scroll-to-bottom" on:click={() => { autoScroll = true; scrollToBottom(); }}>
      ↓ Scroll to bottom
    </button>
  {/if}
</div>

<style>
  .log-view {
    display: flex;
    flex-direction: column;
    height: 100vh;
    overflow: hidden;
  }
  .toolbar {
    display: flex;
    align-items: center;
    gap: 8px;
    padding: 8px 14px;
    border-bottom: 1px solid var(--color-border);
    background: var(--color-card);
    flex-shrink: 0;
  }
  .spacer { flex: 1; }
  .status-pill {
    display: flex;
    align-items: center;
    gap: 5px;
    font-size: 0.78rem;
    font-weight: 600;
    padding: 3px 9px;
    border-radius: 99px;
    text-transform: capitalize;
  }
  .status-pill.idle    { background: #edf2f7; color: #4a5568; }
  .status-pill.running { background: #ebf8ff; color: #2b6cb0; }
  .status-pill.complete{ background: #f0fff4; color: #276749; }
  .status-pill.failed  { background: #fff5f5; color: #c53030; }

  .log-panel {
    flex: 1;
    overflow-y: auto;
    font-family: 'JetBrains Mono', 'Fira Mono', 'Cascadia Code', monospace;
    font-size: 0.78rem;
    line-height: 1.6;
    padding: 10px 14px;
    background: #0f1117;
    color: #d4d4d8;
  }
  .log-line { white-space: pre-wrap; word-break: break-all; }
  .line-stdout { color: #d4d4d8; }
  .line-stderr { color: #a0aec0; }
  .line-warn   { color: #f6e05e; }
  .line-error  { color: #fc8181; }

  .empty { color: #4a5568; text-align: center; margin-top: 48px; }

  .spinner {
    width: 10px; height: 10px;
    border: 2px solid #bee3f8;
    border-top-color: #2b6cb0;
    border-radius: 50%;
    animation: spin 0.7s linear infinite;
  }
  @keyframes spin { to { transform: rotate(360deg); } }

  .scroll-to-bottom {
    position: absolute;
    bottom: 20px;
    right: 20px;
    padding: 6px 12px;
    background: var(--color-primary);
    color: #fff;
    border: none;
    border-radius: 99px;
    font-size: 0.8rem;
    cursor: pointer;
    opacity: 0.9;
  }
  .scroll-to-bottom:hover { opacity: 1; }

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
  .btn-primary {
    padding: 5px 12px;
    background: var(--color-primary);
    color: #fff;
    border: none;
    border-radius: 5px;
    font-size: 0.82rem;
    font-weight: 600;
    cursor: pointer;
  }
  .btn-primary:hover { opacity: 0.88; }
  .btn-danger {
    padding: 5px 12px;
    background: #e53e3e;
    color: #fff;
    border: none;
    border-radius: 5px;
    font-size: 0.82rem;
    font-weight: 600;
    cursor: pointer;
  }
  .btn-danger:hover { background: #c53030; }
  .back-btn { color: var(--color-muted); }
</style>
