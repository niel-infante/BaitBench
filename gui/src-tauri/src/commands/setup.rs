use futures_util::StreamExt;
use serde::{Deserialize, Serialize};
use std::path::{Path, PathBuf};
use std::process::Stdio;
use tauri::{AppHandle, Emitter, Manager};
use tokio::io::{AsyncBufReadExt, BufReader};
use tokio::process::Command as TokioCommand;

// ── Public types ──────────────────────────────────────────────────────────────

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

#[derive(Debug, Serialize, Deserialize)]
pub struct SetupCheck {
    pub conda_found: bool,
    pub conda_executable: Option<String>,
    pub baitbench_env: Option<String>,
    pub baitbench_env_valid: bool,
}

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct SetupProgress {
    pub step: String,
    pub message: String,
    pub percent: Option<f32>,
}

// ── Constants ─────────────────────────────────────────────────────────────────

const REQUIRED_BINS: &[&str] = &["minimap2", "cd-hit-est"];
const OPTIONAL_BINS: &[&str] = &["R", "pandoc", "blastn"];

// ── Shared helpers ────────────────────────────────────────────────────────────

fn emit_progress(app: &AppHandle, step: &str, message: &str, percent: Option<f32>) {
    let _ = app.emit(
        "setup-progress",
        SetupProgress {
            step: step.to_string(),
            message: message.to_string(),
            percent,
        },
    );
}

fn env_from_path(path: PathBuf) -> Option<CondaEnv> {
    if !path.is_dir() {
        return None;
    }
    // Accept env if it has bin/, Scripts/, or a conda executable directly
    let has_bin = path.join("bin").is_dir()
        || path.join("Scripts").is_dir()
        || path.join("conda.exe").exists();
    if !has_bin {
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

/// Check all standard locations for a binary inside a conda env.
fn bin_exists_in_env(env_path: &Path, bin: &str) -> bool {
    let bin_exe = if cfg!(target_os = "windows") {
        format!("{}.exe", bin)
    } else {
        bin.to_string()
    };
    let search = [
        env_path.join("bin").join(&bin_exe),
        env_path.join("bin").join(bin), // fallback without .exe
        env_path.join("Scripts").join(&bin_exe),
        env_path.join("Library").join("bin").join(&bin_exe),
        env_path.join("Library").join("usr").join("bin").join(&bin_exe),
        env_path.join("Library").join("mingw-w64").join("bin").join(&bin_exe),
    ];
    search.iter().any(|p| p.exists())
}

fn find_conda_in_root(root: &Path) -> Option<PathBuf> {
    [
        root.join("bin").join("conda"),
        root.join("condabin").join("conda"),
        root.join("Scripts").join("conda.exe"),
        root.join("conda.exe"),
    ]
    .into_iter()
    .find(|p| p.exists())
}

fn conda_roots(home: &Path) -> Vec<PathBuf> {
    vec![
        home.join("miniforge3"),
        home.join("mambaforge"),
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
        PathBuf::from("/opt/miniforge3"),
    ]
}

fn find_baitbench_env(home: &Path) -> Option<PathBuf> {
    for root in conda_roots(home) {
        let candidate = root.join("envs").join("baitbench");
        if candidate.is_dir() {
            return Some(candidate);
        }
    }
    None
}

fn find_sidecar_binary() -> Result<PathBuf, String> {
    let exe = std::env::current_exe().map_err(|e| e.to_string())?;
    let dir = exe.parent().ok_or("Cannot find exe directory")?;
    std::fs::read_dir(dir)
        .map_err(|e| e.to_string())?
        .filter_map(|e| e.ok())
        .map(|e| e.path())
        .find(|p| {
            p.file_name()
                .unwrap_or_default()
                .to_string_lossy()
                .starts_with("baitbench-")
        })
        .ok_or_else(|| "baitbench sidecar binary not found in app bundle".to_string())
}

fn miniforge_filename() -> Option<&'static str> {
    #[cfg(all(target_os = "macos", target_arch = "aarch64"))]
    return Some("Miniforge3-MacOSX-arm64.sh");
    #[cfg(all(target_os = "macos", target_arch = "x86_64"))]
    return Some("Miniforge3-MacOSX-x86_64.sh");
    #[cfg(all(target_os = "windows", target_arch = "x86_64"))]
    return Some("Miniforge3-Windows-x86_64.exe");
    #[cfg(target_os = "linux")]
    return Some("Miniforge3-Linux-x86_64.sh");
    #[allow(unreachable_code)]
    None
}

async fn download_file(app: &AppHandle, url: &str, dest: &Path) -> Result<(), String> {
    use tokio::io::AsyncWriteExt;

    let response = reqwest::get(url)
        .await
        .map_err(|e| format!("Download failed: {}", e))?;

    if !response.status().is_success() {
        return Err(format!("Download failed: HTTP {}", response.status()));
    }

    let total = response.content_length();
    let mut downloaded: u64 = 0;
    let mut file = tokio::fs::File::create(dest)
        .await
        .map_err(|e| format!("Cannot create temp file: {}", e))?;

    let mut stream = response.bytes_stream();
    while let Some(chunk) = stream.next().await {
        let chunk = chunk.map_err(|e| format!("Download error: {}", e))?;
        file.write_all(&chunk)
            .await
            .map_err(|e| format!("Write error: {}", e))?;
        downloaded += chunk.len() as u64;
        if let Some(total) = total {
            let pct = downloaded as f32 / total as f32;
            let mb = downloaded as f32 / 1_048_576.0;
            emit_progress(app, "download", &format!("Downloading… {:.1} MB", mb), Some(pct));
        }
    }

    Ok(())
}

async fn run_installer(installer: &Path, install_dir: &Path) -> Result<(), String> {
    #[cfg(target_os = "windows")]
    {
        let status = TokioCommand::new(installer)
            .args([
                "/S",
                "/AddToPath=0",
                "/RegisterPython=0",
                &format!("/D={}", install_dir.display()),
            ])
            .status()
            .await
            .map_err(|e| format!("Installer failed: {}", e))?;
        if !status.success() {
            return Err(format!(
                "Miniforge3 installer exited {}",
                status.code().unwrap_or(-1)
            ));
        }
        return Ok(());
    }

    #[cfg(not(target_os = "windows"))]
    {
        use std::os::unix::fs::PermissionsExt;
        let mut perms = std::fs::metadata(installer)
            .map_err(|e| e.to_string())?
            .permissions();
        perms.set_mode(0o755);
        std::fs::set_permissions(installer, perms).map_err(|e| e.to_string())?;

        let status = TokioCommand::new("bash")
            .arg(installer)
            .args(["-b", "-p"])
            .arg(install_dir)
            .status()
            .await
            .map_err(|e| format!("Installer failed: {}", e))?;
        if !status.success() {
            return Err(format!(
                "Miniforge3 installer exited {}",
                status.code().unwrap_or(-1)
            ));
        }
        Ok(())
    }
}

// ── Existing commands (unchanged API) ─────────────────────────────────────────

#[tauri::command]
pub async fn detect_conda_envs() -> Result<Vec<CondaEnv>, String> {
    let home = dirs_next::home_dir().ok_or("Cannot determine home directory")?;
    let mut envs: Vec<CondaEnv> = Vec::new();

    if let Ok(prefix) = std::env::var("CONDA_PREFIX") {
        let p = PathBuf::from(&prefix);
        if let Some(env) = env_from_path(p) {
            envs.push(CondaEnv {
                name: format!("{} (active)", env.name),
                path: env.path,
            });
        }
    }

    for root in conda_roots(&home) {
        if let Some(env) = env_from_path(root.clone()) {
            if !envs.iter().any(|e| e.path == env.path) {
                envs.push(CondaEnv {
                    name: format!("{} (base)", env.name),
                    path: env.path,
                });
            }
        }
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
        if !bin_exists_in_env(base, bin) {
            missing.push(bin.to_string());
        }
    }
    for bin in OPTIONAL_BINS {
        if !bin_exists_in_env(base, bin) {
            warnings.push(format!("{} not found (optional — needed for reports)", bin));
        }
    }
    Ok(ValidationResult {
        valid: missing.is_empty(),
        missing,
        warnings,
    })
}

// ── New wizard commands ───────────────────────────────────────────────────────

/// Step 0: scan for an existing valid baitbench environment.
#[tauri::command]
pub async fn check_setup() -> Result<SetupCheck, String> {
    let home = dirs_next::home_dir().ok_or("Cannot determine home directory")?;

    let mut conda_executable: Option<PathBuf> = None;
    if let Ok(exe) = std::env::var("CONDA_EXE") {
        let p = PathBuf::from(&exe);
        if p.exists() {
            conda_executable = Some(p);
        }
    }
    if conda_executable.is_none() {
        for root in conda_roots(&home) {
            if let Some(p) = find_conda_in_root(&root) {
                conda_executable = Some(p);
                break;
            }
        }
    }

    let baitbench_env = find_baitbench_env(&home);
    let baitbench_env_valid = baitbench_env
        .as_ref()
        .map(|p| REQUIRED_BINS.iter().all(|b| bin_exists_in_env(p, b)))
        .unwrap_or(false);

    Ok(SetupCheck {
        conda_found: conda_executable.is_some(),
        conda_executable: conda_executable.map(|p| p.to_string_lossy().to_string()),
        baitbench_env: baitbench_env.map(|p| p.to_string_lossy().to_string()),
        baitbench_env_valid,
    })
}

/// Step 2: download and silently install Miniforge3.
/// Returns the path to the conda executable.
#[tauri::command]
pub async fn install_conda(app: AppHandle) -> Result<String, String> {
    let home = dirs_next::home_dir().ok_or("Cannot determine home directory")?;
    let install_dir = home.join("miniforge3");

    if install_dir.exists() {
        if let Some(p) = find_conda_in_root(&install_dir) {
            return Ok(p.to_string_lossy().to_string());
        }
    }

    let filename = miniforge_filename()
        .ok_or("Automatic conda install is not supported on this platform")?;
    let url = format!(
        "https://github.com/conda-forge/miniforge/releases/latest/download/{}",
        filename
    );

    emit_progress(&app, "download", &format!("Downloading {}…", filename), Some(0.0));

    let installer_path = std::env::temp_dir().join(filename);
    download_file(&app, &url, &installer_path).await?;

    emit_progress(&app, "install", "Installing Miniforge3 (this takes a moment)…", None);
    run_installer(&installer_path, &install_dir).await?;

    let conda_exe = find_conda_in_root(&install_dir)
        .ok_or("Miniforge3 installed but conda executable not found")?;

    emit_progress(&app, "install", "Conda installed successfully.", Some(1.0));
    Ok(conda_exe.to_string_lossy().to_string())
}

/// Step 3: create the baitbench conda environment from the bundled environment.yml.
/// Returns the env path.
#[tauri::command]
pub async fn create_baitbench_env(
    app: AppHandle,
    conda_executable: String,
) -> Result<String, String> {
    let resource_dir = app.path().resource_dir().map_err(|e| e.to_string())?;
    #[cfg(target_os = "windows")]
    let env_yml = resource_dir.join("environment-windows.yml");
    #[cfg(not(target_os = "windows"))]
    let env_yml = resource_dir.join("environment.yml");

    if !env_yml.exists() {
        return Err(format!(
            "Bundled environment file not found at {}",
            env_yml.display()
        ));
    }

    emit_progress(
        &app,
        "create_env",
        "Creating baitbench environment (this can take 5–10 minutes)…",
        None,
    );

    let mut child = TokioCommand::new(&conda_executable)
        .args(["env", "create", "-f"])
        .arg(&env_yml)
        .args(["-n", "baitbench"])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("Failed to start conda: {}", e))?;

    if let Some(stdout) = child.stdout.take() {
        let app_c = app.clone();
        let mut lines = BufReader::new(stdout).lines();
        tokio::spawn(async move {
            while let Ok(Some(line)) = lines.next_line().await {
                emit_progress(&app_c, "create_env", &line, None);
            }
        });
    }
    if let Some(stderr) = child.stderr.take() {
        let app_c = app.clone();
        let mut lines = BufReader::new(stderr).lines();
        tokio::spawn(async move {
            while let Ok(Some(line)) = lines.next_line().await {
                emit_progress(&app_c, "create_env", &line, None);
            }
        });
    }

    let status = child.wait().await.map_err(|e| e.to_string())?;
    if !status.success() {
        return Err(format!(
            "conda env create failed (exit {})",
            status.code().unwrap_or(-1)
        ));
    }

    let home = dirs_next::home_dir().ok_or("Cannot determine home directory")?;
    let env_path = find_baitbench_env(&home)
        .ok_or("baitbench env was created but its path could not be determined")?;

    #[cfg(target_os = "windows")]
    {
        emit_progress(&app, "win_tools", "Downloading native Windows tools…", None);
        if let Err(e) = fetch_windows_binaries(&app, &env_path).await {
            emit_progress(&app, "win_tools", &format!("Warning: {}", e), None);
        }
    }

    emit_progress(&app, "create_env", "Environment created successfully.", Some(1.0));
    Ok(env_path.to_string_lossy().to_string())
}

/// Done screen: install the bundled sidecar binary to a system PATH location.
/// Returns the install path as a string.
#[tauri::command]
pub async fn install_cli_to_path(_app: AppHandle) -> Result<String, String> {
    let sidecar = find_sidecar_binary()?;

    #[cfg(target_os = "windows")]
    return install_cli_windows(&sidecar);

    #[cfg(not(target_os = "windows"))]
    return install_cli_unix(&sidecar);
}

// ── Platform-specific CLI install ─────────────────────────────────────────────

#[cfg(not(target_os = "windows"))]
fn install_cli_unix(sidecar: &Path) -> Result<String, String> {
    use std::os::unix::fs::PermissionsExt;

    fn copy_and_chmod(src: &Path, dest: &Path) -> Result<(), String> {
        std::fs::copy(src, dest).map_err(|e| e.to_string())?;
        let mut perms = std::fs::metadata(dest)
            .map_err(|e| e.to_string())?
            .permissions();
        perms.set_mode(0o755);
        std::fs::set_permissions(dest, perms).map_err(|e| e.to_string())?;
        Ok(())
    }

    // Prefer /usr/local/bin (typically on PATH, often user-writable on macOS)
    let usr_local = PathBuf::from("/usr/local/bin/baitbench");
    if copy_and_chmod(sidecar, &usr_local).is_ok() {
        return Ok(usr_local.to_string_lossy().to_string());
    }

    // Fall back to ~/.local/bin
    let home = dirs_next::home_dir().ok_or("Cannot determine home directory")?;
    let local_bin = home.join(".local/bin");
    std::fs::create_dir_all(&local_bin).map_err(|e| e.to_string())?;
    let dest = local_bin.join("baitbench");
    copy_and_chmod(sidecar, &dest)?;
    Ok(dest.to_string_lossy().to_string())
}

#[cfg(target_os = "windows")]
fn install_cli_windows(sidecar: &Path) -> Result<String, String> {
    let home = dirs_next::home_dir().ok_or("Cannot determine home directory")?;
    let user_bin = home.join("bin");
    std::fs::create_dir_all(&user_bin).map_err(|e| e.to_string())?;
    let dest = user_bin.join("baitbench.exe");
    std::fs::copy(sidecar, &dest).map_err(|e| e.to_string())?;

    // Append directory to user PATH via PowerShell (idempotent)
    let bin_str = user_bin.to_string_lossy();
    let ps = format!(
        r#"$p=[Environment]::GetEnvironmentVariable('PATH','User'); if ($p -notlike '*{0}*') {{ [Environment]::SetEnvironmentVariable('PATH',"$p;{0}",'User') }}"#,
        bin_str
    );
    let _ = std::process::Command::new("powershell")
        .args(["-NoProfile", "-Command", &ps])
        .output();

    Ok(dest.to_string_lossy().to_string())
}

// ── Windows native binary download ───────────────────────────────────────────

#[cfg(target_os = "windows")]
async fn fetch_windows_binaries(app: &AppHandle, env_path: &Path) -> Result<(), String> {
    use std::io::Read;

    let lib_bin = env_path.join("Library").join("bin");
    std::fs::create_dir_all(&lib_bin).map_err(|e| e.to_string())?;

    let mm2_version = "2.28";
    let mm2_url = format!(
        "https://github.com/lh3/minimap2/releases/download/v{}/minimap2-{}_win64.zip",
        mm2_version, mm2_version
    );
    let mm2_zip = std::env::temp_dir().join("minimap2_win64.zip");

    emit_progress(app, "win_tools", "Downloading minimap2 for Windows…", None);
    download_file(app, &mm2_url, &mm2_zip).await?;

    let zip_data = std::fs::read(&mm2_zip).map_err(|e| e.to_string())?;
    let cursor = std::io::Cursor::new(zip_data);
    let mut archive = zip::ZipArchive::new(cursor).map_err(|e| e.to_string())?;
    for i in 0..archive.len() {
        let mut zf = archive.by_index(i).map_err(|e| e.to_string())?;
        if zf.name().ends_with("minimap2.exe") {
            let dest = lib_bin.join("minimap2.exe");
            let mut buf = Vec::new();
            zf.read_to_end(&mut buf).map_err(|e| e.to_string())?;
            std::fs::write(&dest, &buf).map_err(|e| e.to_string())?;
            break;
        }
    }

    emit_progress(app, "win_tools", "minimap2 installed.", None);
    Ok(())
}
