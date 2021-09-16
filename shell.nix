{ pkgs ? import <nixpkgs> {} }:
with pkgs;
mkShell {
  name = "ecmv";
  buildInputs = [
    gcc gnumake m4 binutils gmp libmpc mpfr gnuplot
    cargo hwloc rustc rustfmt
  ];
}
