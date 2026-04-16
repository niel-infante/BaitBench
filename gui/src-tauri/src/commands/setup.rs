use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct CondaEnv {
    pub name: String,
    pub path: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct ValidationResult {
    pub valid: bool,
    pub missing: Vec<String>,
    pub warnings: Vec<String>,
}

/// Required binaries that must be present in a valid conda env.
const REQUIRED_BINS: &[&str] = &["minimap2", "cd-hit-est"];
/// Optional binaries — absence is a warning, not an error.
const OPTIONAL_BINS: &[&str] = &["R", "pandoc", "blastn"];

fn env_from_path(path: PathBuf) -> Option<CondaEnv> {
    if !path.is_dir() {
        return None;
    }
    let bin = path.join("bin");
    // A valid conda env must have a bin/ directory
    if !bin.is_dir() {
        return None;
    }
    let name = path
        .file_name()
        .and_then(|n| n.to_str())
        .unwrap_or("unknown")
        .to_string();
    Some(CondaEnv {
        name,
        path: path.to_string_lossy().to_string(),
    })
}

#[tauri::command]
pub async fn detect_conda_envs() -> Result<Vec<CondaEnv>, String> {
    let home = dirs_next::home_dir().ok_or("Cannot determine home directory")?;

    // Candidate root conda directories
    let roots: Vec<PathBuf> = vec![
        home.join("miniconda3"),
        home.join("miniconda"),
        home.join("anaconda3"),
        home.join("anaconda"),
        home.join(".conda"),
        home.join("opt/miniconda3"),
        home.join("opt/anaconda3"),
        PathBuf::from("/opt/homebrew/Caskroom/miniconda/base"),
        PathBuf::from("/usr/local/anaconda3"),
        PathBuf::from("/opt/anaconda3"),
    ];

    let mut envs: Vec<CondaEnv> = Vec::new();

    // Check CONDA_PREFIX env var (active conda env)
    if let Ok(prefix) = std::env::var("CONDA_PREFIX") {
        let p = PathBuf::from(&prefix);
        if let Some(env) = env_from_path(p) {
            envs.push(CondaEnv {
                name: format!("{} (active)", env.name),
                path: env.path,
            });
        }
    }

    for root in &roots {
        // Base env (root itself)
        if let Some(env) = env_from_path(root.clone()) {
            if !envs.iter().any(|e| e.path == env.path) {
                envs.push(CondaEnv {
                    name: format!("{} (base)", env.name),
                    path: env.path,
                });
            }
        }

        // Named envs under envs/
        let envs_dir = root.join("envs");
        if let Ok(entries) = std::fs::read_dir(&envs_dir) {
            for entry in entries.flatten() {
                let p = entry.path();
                if let Some(env) = env_from_path(p) {
                    if !envs.iter().any(|e| e.path == env.path) {
                        envs.push(env);
                    }
                }
            }
        }
    }

    Ok(envs)
}

#[tauri::command]
pub async fn validate_conda_env(path: String) -> Result<ValidationResult, String> {
    let base = Path::new(&path);
    let bin_dir = base.join("bin");

    if !base.is_dir() {
        return Ok(ValidationResult {
            valid: false,
            missing: vec!["Environment directory does not exist".to_string()],
            warnings: vec![],
        });
    }

    let mut missing = Vec::new();
    let mut warnings = Vec::new();

    for bin in REQUIRED_BINS {
        let p = bin_dir.join(bin);
        if !p.exists() {
            missing.push(bin.to_string());
        }
    }

    for bin in OPTIONAL_BINS {
        let p = bin_dir.join(bin);
        if !p.exists() {
            warnings.push(format!("{} not found (optional — needed for reports)", bin));
        }
    }

    Ok(ValidationResult {
        valid: missing.is_empty(),
        missing,
        warnings,
    })
}
