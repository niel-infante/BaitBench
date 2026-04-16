cd gui
make release
That single command:

Builds baitbench CLI (cargo build --release)
Copies it as the Tauri sidecar
Runs npm run tauri:build (produces .app + .dmg)
Runs package-release.mjs which:
Ad-hoc signs the .app with codesign --sign -
Copies both the .app/.dmg and the standalone baitbench CLI into gui/dist-release/
What users get:

File	Purpose
BaitBench.app	GUI — drag to Applications, double-click to launch
BaitBench_x.x.x_aarch64.dmg	Installer for the GUI
baitbench	Standalone CLI — copy to /usr/local/bin
Gatekeeper (unsigned distribution): Since the app isn't notarized, users need to right-click → Open on first launch, or run xattr -d com.apple.quarantine BaitBench.app. The script prints this reminder. If you get an Apple Developer account later, xcrun notarytool submit is the path to fully transparent installs.




Yes. Tauri uses a set of pre-sized PNG/ICNS files generated from a single source image.

**Source image requirements:**

* Square (1:1 aspect ratio)
* At least **1024×1024 pixels**
* **RGBA PNG** (transparency supported but not required)

**How to set it:**

Tauri has a built-in command that takes your source image and generates all the required icon sizes automatically:

```bash
cd gui
npm run tauri icon /path/to/your-icon.png
```

This overwrites `gui/src-tauri/icons/` with the full set of correctly-sized files (`32x32.png`, `128x128.png`, `128x128@2x.png`, `icon.icns` for macOS, `icon.ico` for Windows, etc.) and updates `tauri.conf.json` to reference them.

Then rebuild:

```bash
make release
```



