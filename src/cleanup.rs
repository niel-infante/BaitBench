use std::fs;
use std::path::Path;

/// Delete a list of files from the given directory. Missing files are silently skipped.
pub fn cleanup_files(outdir: &Path, filenames: &[&str]) {
    for name in filenames {
        let path = outdir.join(name);
        if path.exists() {
            match fs::remove_file(&path) {
                Ok(()) => log::info!("Cleanup: removed {}", path.display()),
                Err(e) => log::warn!("Cleanup: failed to remove {}: {}", path.display(), e),
            }
        }
    }
}

/// Remove all subdirectories within the given directory.
pub fn cleanup_all_subdirs(outdir: &Path) {
    let entries = match fs::read_dir(outdir) {
        Ok(e) => e,
        Err(e) => {
            log::warn!("Cleanup: cannot read directory {}: {}", outdir.display(), e);
            return;
        }
    };

    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            match fs::remove_dir_all(&path) {
                Ok(()) => log::info!("Cleanup: removed directory {}", path.display()),
                Err(e) => log::warn!("Cleanup: failed to remove directory {}: {}", path.display(), e),
            }
        }
    }
}
