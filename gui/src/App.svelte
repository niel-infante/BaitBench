<script lang="ts">
  import { onMount } from 'svelte';
  import { load } from '@tauri-apps/plugin-store';
  import { currentView, condaEnvPath } from './lib/stores';
  import SetupView from './lib/views/SetupView.svelte';
  import RunView from './lib/views/RunView.svelte';
  import LogView from './lib/views/LogView.svelte';

  onMount(async () => {
    try {
      const store = await load('settings.json', { autoSave: true });
      const savedEnv = await store.get<string>('conda_env_path');
      if (savedEnv) {
        condaEnvPath.set(savedEnv);
        currentView.set('run');
      } else {
        currentView.set('setup');
      }
    } catch {
      currentView.set('setup');
    }
  });

  // Keep references visible to Svelte's checker so it doesn't tree-shake the imports.
  // These are never rendered — they only satisfy the static analyser.
  const _views = [SetupView, RunView, LogView];
</script>

{#if $currentView === 'setup'}
  <SetupView />
{:else if $currentView === 'run'}
  <RunView />
{:else if $currentView === 'log'}
  <LogView />
{/if}
