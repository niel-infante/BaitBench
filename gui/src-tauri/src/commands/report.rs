use tauri::{AppHandle, WebviewWindowBuilder, WebviewUrl};

#[tauri::command]
pub async fn open_report(path: String, app: AppHandle) -> Result<(), String> {
    // Convert to a file:// URL
    let url = format!("file://{}", path);
    let webview_url = WebviewUrl::External(url.parse().map_err(|e: url::ParseError| e.to_string())?);

    WebviewWindowBuilder::new(&app, "report", webview_url)
        .title("BaitBench Report")
        .inner_size(1200.0, 900.0)
        .min_inner_size(800.0, 600.0)
        .resizable(true)
        .build()
        .map_err(|e| e.to_string())?;

    Ok(())
}
