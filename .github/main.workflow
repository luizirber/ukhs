workflow "Quickstart" {
  on = "push"
  resolves = ["quickstart"]
}

action "quickstart" {
  uses = "icepuma/rust-action@1.0.4"
  args = "cargo fmt -- --check && cargo clippy -- -Dwarnings && cargo test && cargo bench"
}
