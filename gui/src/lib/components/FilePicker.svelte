<script lang="ts">
  import { open } from '@tauri-apps/plugin-dialog';

  export let label: string;
  export let value: string = '';
  export let placeholder: string = 'Select file…';
  export let directory: boolean = false;
  export let multiple: boolean = false;
  export let filters: { name: string; extensions: string[] }[] = [];
  export let required: boolean = false;
  export let id: string = label.replace(/\s+/g, '_').toLowerCase();

  async function browse() {
    const selected = await open({
      multiple,
      directory,
      filters: filters.length > 0 ? filters : undefined,
    });
    if (selected === null) return;
    if (multiple && Array.isArray(selected)) {
      value = selected.join('\t');
    } else if (typeof selected === 'string') {
      value = selected;
    }
  }
</script>

<div class="file-picker">
  <label for={id}>
    {label}
    {#if required}<span class="req">*</span>{/if}
  </label>
  <div class="input-row">
    <input
      {id}
      type="text"
      bind:value
      {placeholder}
      class:multi={multiple}
      readonly
    />
    <button type="button" on:click={browse}>Browse…</button>
  </div>
  {#if multiple && value}
    <div class="multi-paths">
      {#each value.split('\t').filter(Boolean) as p}
        <span class="path-chip">{p.split('/').at(-1)}</span>
      {/each}
    </div>
  {/if}
</div>

<style>
  .file-picker {
    display: flex;
    flex-direction: column;
    gap: 4px;
  }
  label {
    font-size: 0.85rem;
    font-weight: 600;
    color: var(--color-label);
  }
  .req {
    color: #e53e3e;
    margin-left: 2px;
  }
  .input-row {
    display: flex;
    gap: 6px;
  }
  input {
    flex: 1;
    font-size: 0.85rem;
    padding: 5px 8px;
    border: 1px solid var(--color-border);
    border-radius: 4px;
    background: var(--color-input-bg);
    color: var(--color-text);
    min-width: 0;
    cursor: default;
  }
  button {
    padding: 5px 10px;
    font-size: 0.8rem;
    border: 1px solid var(--color-border);
    border-radius: 4px;
    background: var(--color-btn);
    color: var(--color-text);
    cursor: pointer;
    white-space: nowrap;
  }
  button:hover {
    background: var(--color-btn-hover);
  }
  .multi-paths {
    display: flex;
    flex-wrap: wrap;
    gap: 4px;
    margin-top: 2px;
  }
  .path-chip {
    background: var(--color-chip);
    border-radius: 3px;
    padding: 1px 6px;
    font-size: 0.78rem;
    color: var(--color-label);
  }
</style>
