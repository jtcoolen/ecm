## How to use

If using the Nix package manager, simply enter `nix-shell` in your shell to get the required dependencies.

Otherwise, the following dependencies are needed:
- Rust
- GMP (as well as libmpc and mpfr, although not used in the project)
- gnuplot
- hwloc

Compile with `cargo build`. The executable is found in `./target/debug/main`. Type `./target/debug/main -h` for usage instructions.
Run test suites with `cargo test` and benchmarks with `cargo bench`.

Program multi-threaded by default.
