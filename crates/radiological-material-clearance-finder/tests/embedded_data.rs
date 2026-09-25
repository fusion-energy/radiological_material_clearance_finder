//! The tables under `data/` are compiled into the library with `include_str!`,
//! so every byte of them ships in every binary that links this crate. A
//! provenance field holding an absolute path would put the machine the table
//! was generated on into all of them: yamc-core's wheel check rejected 0.2.0
//! for exactly that, three `/home/<user>/openmc/...` strings.

use std::path::Path;

/// Prefixes that only ever name a particular machine's filesystem.
const MACHINE_PATHS: [&str; 5] = ["/home/", "/Users/", "/root/", "/tmp/", ":\\"];

fn json_files(dir: &Path, found: &mut Vec<std::path::PathBuf>) {
    for entry in std::fs::read_dir(dir).expect("data directory is readable") {
        let path = entry.expect("directory entry").path();
        if path.is_dir() {
            json_files(&path, found);
        } else if path.extension().is_some_and(|e| e == "json") {
            found.push(path);
        }
    }
}

#[test]
fn no_embedded_table_names_a_local_path() {
    let mut files = Vec::new();
    json_files(
        &Path::new(env!("CARGO_MANIFEST_DIR")).join("data"),
        &mut files,
    );
    assert!(
        files.len() >= 9,
        "expected the nine shipped tables, found {files:?}"
    );

    let mut leaks = Vec::new();
    for file in &files {
        let text = std::fs::read_to_string(file).expect("table is UTF-8");
        for (number, line) in text.lines().enumerate() {
            if MACHINE_PATHS.iter().any(|p| line.contains(p)) {
                leaks.push(format!(
                    "{}:{}: {}",
                    file.display(),
                    number + 1,
                    line.trim()
                ));
            }
        }
    }
    assert!(
        leaks.is_empty(),
        "local paths would ship in the binary:\n{}",
        leaks.join("\n")
    );
}
