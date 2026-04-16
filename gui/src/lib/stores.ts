import { writable } from 'svelte/store';
import type { PipelineStatus, ViewName, LogLine } from './types';

export const currentView = writable<ViewName>('setup');
export const pipelineStatus = writable<PipelineStatus>('idle');
export const logLines = writable<LogLine[]>([]);
export const reportPath = writable<string | null>(null);
export const condaEnvPath = writable<string>('');
