export interface CondaEnv {
  name: string;
  path: string;
}

export interface SetupCheck {
  conda_found: boolean;
  conda_executable: string | null;
  baitbench_env: string | null;
  baitbench_env_valid: boolean;
}

export interface SetupProgress {
  step: string;
  message: string;
  percent: number | null;
}

export interface ValidationResult {
  valid: boolean;
  missing: string[];
  warnings: string[];
}

export interface LogLine {
  text: string;
  stream: 'stdout' | 'stderr';
}

export type PipelineStatus = 'idle' | 'running' | 'complete' | 'failed';

export type ViewName = 'setup' | 'run' | 'log';

/** Represents the full configuration for a pipeline invocation. */
export interface PipelineConfig {
  tool: string;
  /** flag → value. Boolean flags use "". Multi-value flags use tab-separated string. */
  args: Record<string, string>;
  conda_env: string;
}

/** One tool definition used to drive the RunView form. */
export interface ToolDef {
  id: string;
  label: string;
  category: string;
  description: string;
}
