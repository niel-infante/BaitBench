import { defineConfig } from 'vite';
import { svelte } from '@sveltejs/vite-plugin-svelte';
import preprocess from 'svelte-preprocess';

// https://vitejs.dev/config/
export default defineConfig(async () => ({
  plugins: [svelte({ preprocess: preprocess() })],

  // Vite options tailored for Tauri development and only applied in `tauri dev` command
  clearScreen: false,
  server: {
    port: 1420,
    strictPort: true,
    watch: {
      // Tell Vite to ignore watching `src-tauri`
      ignored: ['**/src-tauri/**'],
    },
  },
}));
