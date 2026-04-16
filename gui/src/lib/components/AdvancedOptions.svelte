<script lang="ts">
  export let label: string = 'Advanced Options';
  let open = false;
</script>

<details bind:open class="advanced">
  <summary class="advanced-toggle">
    <span class="arrow">{open ? '▾' : '▸'}</span>
    {label}
  </summary>
  <div class="advanced-body">
    <slot />
  </div>
</details>

<style>
  .advanced {
    margin-top: 8px;
    border: 1px solid var(--color-border);
    border-radius: 6px;
  }
  .advanced-toggle {
    list-style: none;
    display: flex;
    align-items: center;
    gap: 6px;
    padding: 8px 12px;
    font-size: 0.85rem;
    font-weight: 600;
    cursor: pointer;
    user-select: none;
    color: var(--color-label);
  }
  .advanced-toggle::-webkit-details-marker {
    display: none;
  }
  .arrow {
    font-size: 0.75rem;
    color: var(--color-muted);
  }
  .advanced-body {
    padding: 12px 14px;
    display: flex;
    flex-direction: column;
    gap: 12px;
    border-top: 1px solid var(--color-border);
  }
</style>
