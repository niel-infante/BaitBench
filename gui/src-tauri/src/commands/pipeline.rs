use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::env;
use std::sync::Mutex;
use tauri::{AppHandle, Emitter, Manager, State};
use tauri_plugin_notification::NotificationExt;
use tauri_plugin_shell::process::CommandEvent;
use tauri_plugin_shell::ShellExt;

#[derive(Debug, Serialize, Deserialize, Clone, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum PipelineStatus {
    Idle,
    Running,
    Complete,
    Failed,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct PipelineConfig {
    /// e.g. "run", "build-probes", "xreact", …
    pub tool: String,
    /// Flag → value pairs. Flags with no value (booleans) use "" as value.
    /// Multi-value flags (e.g. --distractors) use a tab-separated string.
    pub args: HashMap<String, String>,
    /// Full path to the conda env to prepend to PATH.
    pub conda_env: String,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct LogLine {
    pub text: String,
    pub stream: String, // "stdout" | "stderr"
}

pub struct AppState {
    pub status: Mutex<PipelineStatus>,
    // We store the pid of the child so we can kill it.
    pub child_pid: Mutex<Option<u32>>,
}

impl Default for AppState {
    fn default() -> Self {
        Self {
            status: Mutex::new(PipelineStatus::Idle),
            child_pid: Mutex::new(None),
        }
    }
}

/// Build the argv list from a PipelineConfig.
/// Each flag → one or more args in the correct order.
fn build_argv(config: &PipelineConfig) -> Vec<String> {
    let mut argv: Vec<String> = vec![config.tool.clone()];

    for (flag, value) in &config.args {
        if value.is_empty() {
            // Boolean flag
            argv.push(flag.clone());
        } else if value.contains('\t') {
            // Multi-value flag (tab-separated)
            argv.push(flag.clone());
            for v in value.split('\t') {
                if !v.is_empty() {
                    argv.push(v.to_string());
                }
            }
        } else {
            argv.push(flag.clone());
            argv.push(value.clone());
        }
    }
    argv
}

#[tauri::command]
pub async fn run_pipeline(
    config: PipelineConfig,
    app: AppHandle,
    state: State<'_, AppState>,
) -> Result<(), String> {
    // Guard: only one pipeline at a time
    {
        let status = state.status.lock().unwrap();
        if *status == PipelineStatus::Running {
            return Err("A pipeline is already running".to_string());
        }
    }

    // Set status to Running
    *state.status.lock().unwrap() = PipelineStatus::Running;

    // Build PATH with conda env prepended
    let conda_bin = format!("{}/bin", config.conda_env);
    let current_path = env::var("PATH").unwrap_or_default();
    let new_path = if current_path.is_empty() {
        conda_bin.clone()
    } else {
        format!("{}:{}", conda_bin, current_path)
    };

    // The sidecar binary IS baitbench; all args including the subcommand name are passed as argv.
    let argv = build_argv(&config);

    let shell = app.shell();
    let cmd = shell
        .sidecar("baitbench")
        .map_err(|e| e.to_string())?
        .env("PATH", &new_path)
        .args(&argv);

    let (mut rx, child) = cmd.spawn().map_err(|e| e.to_string())?;

    // Store pid for cancellation
    *state.child_pid.lock().unwrap() = Some(child.pid());

    // Stream events on a separate task so the command returns immediately
    let app_clone = app.clone();
    tokio::spawn(async move {
        let mut success = true;

        while let Some(event) = rx.recv().await {
            match event {
                CommandEvent::Stdout(line) => {
                    let text = String::from_utf8_lossy(&line).to_string();
                    let _ = app_clone.emit(
                        "pipeline-log",
                        LogLine {
                            text,
                            stream: "stdout".to_string(),
                        },
                    );
                }
                CommandEvent::Stderr(line) => {
                    let text = String::from_utf8_lossy(&line).to_string();
                    let _ = app_clone.emit(
                        "pipeline-log",
                        LogLine {
                            text,
                            stream: "stderr".to_string(),
                        },
                    );
                }
                CommandEvent::Error(err) => {
                    success = false;
                    let _ = app_clone.emit(
                        "pipeline-log",
                        LogLine {
                            text: format!("[error] {}", err),
                            stream: "stderr".to_string(),
                        },
                    );
                }
                CommandEvent::Terminated(payload) => {
                    let code = payload.code.unwrap_or(-1);
                    if code != 0 {
                        success = false;
                    }
                    break;
                }
                _ => {}
            }
        }

        // Update state and notify
        let state = app_clone.state::<AppState>();
        *state.child_pid.lock().unwrap() = None;

        if success {
            *state.status.lock().unwrap() = PipelineStatus::Complete;
            let _ = app_clone.emit("pipeline-complete", ());
            let _ = app_clone
                .notification()
                .builder()
                .title("BaitBench")
                .body("Pipeline finished successfully")
                .show();
        } else {
            *state.status.lock().unwrap() = PipelineStatus::Failed;
            let _ = app_clone.emit("pipeline-error", ());
            let _ = app_clone
                .notification()
                .builder()
                .title("BaitBench")
                .body("Pipeline finished with errors")
                .show();
        }
    });

    Ok(())
}

#[tauri::command]
pub fn get_pipeline_status(state: State<'_, AppState>) -> PipelineStatus {
    state.status.lock().unwrap().clone()
}

#[tauri::command]
pub fn cancel_pipeline(state: State<'_, AppState>) -> Result<(), String> {
    let pid_opt = state.child_pid.lock().unwrap().clone();
    if let Some(pid) = pid_opt {
        // Send SIGTERM on Unix
        #[cfg(unix)]
        unsafe {
            libc::kill(pid as libc::pid_t, libc::SIGTERM);
        }
        #[cfg(windows)]
        {
            // On Windows, use taskkill
            let _ = std::process::Command::new("taskkill")
                .args(["/PID", &pid.to_string(), "/F"])
                .output();
        }
        *state.status.lock().unwrap() = PipelineStatus::Failed;
        *state.child_pid.lock().unwrap() = None;
        Ok(())
    } else {
        Err("No pipeline is running".to_string())
    }
}
