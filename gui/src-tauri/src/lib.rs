mod commands;

use commands::pipeline::AppState;

pub fn run() {
    tauri::Builder::default()
        .plugin(tauri_plugin_shell::init())
        .plugin(tauri_plugin_dialog::init())
        .plugin(tauri_plugin_fs::init())
        .plugin(tauri_plugin_store::Builder::default().build())
        .plugin(tauri_plugin_notification::init())
        .manage(AppState::default())
        .invoke_handler(tauri::generate_handler![
            commands::setup::detect_conda_envs,
            commands::setup::validate_conda_env,
            commands::pipeline::run_pipeline,
            commands::pipeline::get_pipeline_status,
            commands::pipeline::cancel_pipeline,
            commands::report::open_report,
        ])
        .run(tauri::generate_context!())
        .expect("error while running tauri application");
}
