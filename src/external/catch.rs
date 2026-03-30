use anyhow::{bail, Context, Result};
use std::fs::File;
use std::path::Path;
use std::path::PathBuf;
use std::process::Command;

/// Check that CATCH's design.py and python are available on PATH.
pub fn check_available() -> Result<()> {
    if which_design_py().is_none() {
        bail!(
            "CATCH design.py not found on PATH. \
             Install via: conda install -c bioconda catch"
        );
    }

    let python_ok = Command::new("python")
        .arg("--version")
        .stdout(std::process::Stdio::null())
        .stderr(std::process::Stdio::null())
        .status()
        .map(|s| s.success())
        .unwrap_or(false);

    if !python_ok {
        bail!("python not found on PATH (required to run CATCH design.py)");
    }

    Ok(())
}

/// Locate design.py on PATH by running `which design.py`.
fn which_design_py() -> Option<PathBuf> {
    Command::new("which")
        .arg("design.py")
        .output()
        .ok()
        .and_then(|out| {
            if out.status.success() {
                let path = String::from_utf8_lossy(&out.stdout).trim().to_string();
                if path.is_empty() {
                    None
                } else {
                    Some(PathBuf::from(path))
                }
            } else {
                None
            }
        })
}

/// Run CATCH's design.py to design probes from target sequences.
///
/// ```text
/// design.py <input> -pl <probe_length> -o <output> [extra_args...]
/// ```
///
/// On macOS with Python >= 3.13, the default multiprocessing start method
/// changed from `fork` to `forkserver`/`spawn`. CATCH relies on fork semantics
/// to share global variables with worker processes, so we force `fork` mode
/// via the `OBJC_DISABLE_INITIALIZE_FORK_SAFETY` env var and a Python wrapper.
pub fn design(
    input: &Path,
    output: &Path,
    probe_length: usize,
    extra_args: &str,
    log_file: &Path,
) -> Result<()> {
    let log = File::create(log_file)
        .with_context(|| format!("Cannot create log: {}", log_file.display()))?;

    // Find design.py on PATH so we can invoke it via python
    let design_py = which_design_py()
        .context("Cannot locate design.py on PATH")?;

    // Use python to run design.py with fork mode forced before any imports
    let wrapper = format!(
        "import multiprocessing; multiprocessing.set_start_method('fork'); exec(open({:?}).read())",
        design_py.to_string_lossy()
    );

    let mut cmd = Command::new("python");
    cmd.arg("-c")
        .arg(&wrapper)
        // design.py positional/flag args follow
        .arg(input)
        .arg("-pl")
        .arg(probe_length.to_string())
        .arg("-o")
        .arg(output);

    // Parse extra_args string into individual arguments
    for arg in shell_words(extra_args) {
        if !arg.is_empty() {
            cmd.arg(&arg);
        }
    }

    // macOS safety: allow fork in Objective-C runtime
    cmd.env("OBJC_DISABLE_INITIALIZE_FORK_SAFETY", "YES");

    cmd.stdout(log.try_clone()?)
        .stderr(log);

    let status = cmd
        .status()
        .context("Failed to execute CATCH design.py")?;

    if !status.success() {
        bail!(
            "CATCH design.py failed (exit code {:?}). Check log: {}",
            status.code(),
            log_file.display()
        );
    }

    // Verify output was created
    if !output.exists() {
        bail!(
            "CATCH design.py completed but output file not found: {}",
            output.display()
        );
    }

    Ok(())
}

/// Split a string into shell-like words, respecting quoted strings.
/// Simple implementation: splits on whitespace, handles double and single quotes.
fn shell_words(s: &str) -> Vec<String> {
    let mut words = Vec::new();
    let mut current = String::new();
    let mut in_single_quote = false;
    let mut in_double_quote = false;

    for ch in s.chars() {
        match ch {
            '\'' if !in_double_quote => {
                in_single_quote = !in_single_quote;
            }
            '"' if !in_single_quote => {
                in_double_quote = !in_double_quote;
            }
            ' ' | '\t' if !in_single_quote && !in_double_quote => {
                if !current.is_empty() {
                    words.push(std::mem::take(&mut current));
                }
            }
            _ => {
                current.push(ch);
            }
        }
    }
    if !current.is_empty() {
        words.push(current);
    }

    words
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_shell_words_basic() {
        assert_eq!(
            shell_words("-ps 60 -m 5 -e 0 --filter-with-lsh-minhash 0.6"),
            vec!["-ps", "60", "-m", "5", "-e", "0", "--filter-with-lsh-minhash", "0.6"]
        );
    }

    #[test]
    fn test_shell_words_empty() {
        assert_eq!(shell_words(""), Vec::<String>::new());
    }

    #[test]
    fn test_shell_words_quoted() {
        assert_eq!(
            shell_words("--arg \"value with spaces\""),
            vec!["--arg", "value with spaces"]
        );
    }
}
