# termocline-grid

The single definition of "what a grid cell is" for this project: [`Field2D`], a
flat row-major 2D array, and [`Grid`], the Arakawa C-grid geometry that says
where `h`, `u` and `v` live relative to a cell.

It holds data structures and indexing math only — no physics, no time
stepping, no I/O. The staggering conventions it encodes come from
[ADR-0003](../docs/planning/adr/0003-numerical-scheme.md): `h` at cell centers,
`u` at east/west faces, `v` at north/south faces. Solver code uses the named
offsets from this crate instead of raw `+1`/`-1` index arithmetic.

[`Field2D`]: src/lib.rs
[`Grid`]: src/lib.rs
