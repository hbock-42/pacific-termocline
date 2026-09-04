//! Acceptance tests for T-00.4.
//!
//! Two acceptance criteria, both about indexing only — this crate carries no
//! physics, so there is no numerical scheme here and nothing whose order of
//! accuracy could be measured. Expected values come from the row-major
//! definition written out by hand and from the staggering stated in
//! [ADR-0003](../../docs/planning/adr/0003-numerical-scheme.md), never from
//! running the code.

use termocline_grid::{Field2D, Grid, Staggering, H_STAGGERING, U_STAGGERING, V_STAGGERING};

/// A deliberately non-square shape, so a transposed index formula cannot pass
/// by symmetry.
const NX: usize = 7;
const NY: usize = 3;

// --- Acceptance criterion 1: (i, j) -> flat -> (i, j) round-trips. ---

#[test]
fn index_round_trips_for_every_cell() {
    let field = Field2D::filled(NX, NY, 0.0_f64).expect("7x3 is a valid shape");

    for j in 0..NY {
        for i in 0..NX {
            let flat = field.flat_index(i, j).expect("in-bounds cell");
            assert_eq!(
                field.cell_of(flat),
                Some((i, j)),
                "flat index {flat} should map back to ({i}, {j})"
            );
        }
    }
}

#[test]
fn flat_index_is_row_major() {
    let field = Field2D::filled(NX, NY, 0u8).expect("7x3 is a valid shape");

    // Row-major: consecutive `i` are adjacent in memory, a step in `j` skips a
    // whole row. Written out independently of the implementation.
    for j in 0..NY {
        for i in 0..NX {
            assert_eq!(field.flat_index(i, j), Some(j * NX + i));
        }
    }
}

#[test]
fn every_flat_offset_round_trips_to_a_distinct_cell() {
    let field = Field2D::filled(NX, NY, 0u8).expect("7x3 is a valid shape");
    assert_eq!(field.len(), NX * NY);

    let cells: Vec<(usize, usize)> = (0..field.len())
        .map(|flat| field.cell_of(flat).expect("in-bounds offset"))
        .collect();

    for (flat, &(i, j)) in cells.iter().enumerate() {
        assert_eq!(field.flat_index(i, j), Some(flat));
    }

    let mut distinct = cells.clone();
    distinct.sort_unstable();
    distinct.dedup();
    assert_eq!(
        distinct.len(),
        field.len(),
        "each of the {} offsets must address a distinct cell",
        field.len()
    );
}

#[test]
fn out_of_bounds_indices_have_no_flat_offset() {
    let field = Field2D::filled(NX, NY, 0u8).expect("7x3 is a valid shape");

    assert_eq!(field.flat_index(NX, 0), None);
    assert_eq!(field.flat_index(0, NY), None);
    assert_eq!(field.cell_of(NX * NY), None);
    assert_eq!(field.get(NX, NY - 1), None);
}

#[test]
fn values_are_stored_and_read_back_at_the_cell_they_were_written_to() {
    let mut field = Field2D::filled(NX, NY, 0.0_f64).expect("7x3 is a valid shape");
    *field.get_mut(5, 2).expect("in-bounds cell") = 1.5;

    assert_eq!(field.get(5, 2), Some(&1.5));
    // Neighbours must be untouched: a row-major/column-major mix-up would
    // write into (2, 5)'s slot instead.
    assert_eq!(field.get(2, 1), Some(&0.0));
    assert_eq!(field.as_slice()[2 * NX + 5], 1.5);
}

#[test]
fn degenerate_and_mismatched_shapes_are_rejected() {
    assert!(Field2D::filled(0, NY, 0.0_f64).is_err());
    assert!(Field2D::filled(NX, 0, 0.0_f64).is_err());
    assert!(Field2D::from_vec(NX, NY, vec![0.0_f64; NX * NY - 1]).is_err());
    assert!(Field2D::from_vec(NX, NY, vec![0.0_f64; NX * NY]).is_ok());
}

// --- Acceptance criterion 2: C-grid staggering is named, not magic. ---

#[test]
fn staggering_of_each_prognostic_variable_matches_adr_0003() {
    // ADR-0003: "`h` at cell centers, `u` at cell east/west faces, `v` at cell
    // north/south faces."
    assert_eq!(H_STAGGERING, Staggering::CellCenter);
    assert_eq!(U_STAGGERING, Staggering::EastWestFace);
    assert_eq!(V_STAGGERING, Staggering::NorthSouthFace);
}

#[test]
fn staggering_offsets_are_half_a_cell_on_the_staggered_axis() {
    // Offsets are in cell widths, measured from the cell's southwest corner.
    // On an Arakawa C-grid the center sits half a cell in on both axes, and a
    // face value sits on the corner in the direction it is staggered.
    assert_eq!(Staggering::CellCenter.offset_in_cells(), (0.5, 0.5));
    assert_eq!(Staggering::EastWestFace.offset_in_cells(), (0.0, 0.5));
    assert_eq!(Staggering::NorthSouthFace.offset_in_cells(), (0.5, 0.0));
}

#[test]
fn face_fields_carry_one_extra_line_of_points_on_their_staggered_axis() {
    let grid = Grid::new(NX, NY).expect("valid grid");

    // A closed basin of NX by NY cells has NX + 1 east/west faces per row and
    // NY + 1 north/south faces per column.
    assert_eq!(grid.field_shape(H_STAGGERING), (NX, NY));
    assert_eq!(grid.field_shape(U_STAGGERING), (NX + 1, NY));
    assert_eq!(grid.field_shape(V_STAGGERING), (NX, NY + 1));
}

#[test]
fn allocated_fields_match_the_shape_their_staggering_asks_for() {
    let grid = Grid::new(NX, NY).expect("valid grid");

    let h: Field2D<f64> = grid.allocate(H_STAGGERING, 0.0);
    let u: Field2D<f64> = grid.allocate(U_STAGGERING, 0.0);
    let v: Field2D<f64> = grid.allocate(V_STAGGERING, 0.0);

    assert_eq!((h.nx(), h.ny()), (NX, NY));
    assert_eq!((u.nx(), u.ny()), (NX + 1, NY));
    assert_eq!((v.nx(), v.ny()), (NX, NY + 1));
    assert_eq!(h.len(), NX * NY);
}

#[test]
fn offsets_stay_dimensionless_so_the_face_line_is_addressable() {
    let grid = Grid::new(NX, NY).expect("valid grid");

    // Offsets are fractions of a cell, not lengths: turning one into metres
    // needs a cell spacing, which is a scenario parameter and lives outside
    // this crate. Cell (3, 1)'s u point is on its western edge, three whole
    // cells east of the basin's western boundary and half a cell north of
    // cell row 1's southern edge.
    let (offset_x, offset_y) = U_STAGGERING.offset_in_cells();
    assert_eq!((3.0 + offset_x, 1.0 + offset_y), (3.0, 1.5));
    let (offset_x, offset_y) = V_STAGGERING.offset_in_cells();
    assert_eq!((3.0 + offset_x, 1.0 + offset_y), (3.5, 1.0));

    // The extra face line is addressable: u at i = NX is the basin's eastern
    // boundary, v at j = NY its northern one; h has no such line.
    let u: Field2D<f64> = grid.allocate(U_STAGGERING, 0.0);
    let v: Field2D<f64> = grid.allocate(V_STAGGERING, 0.0);
    let h: Field2D<f64> = grid.allocate(H_STAGGERING, 0.0);
    assert!(u.flat_index(NX, 0).is_some());
    assert!(v.flat_index(0, NY).is_some());
    assert_eq!(u.flat_index(NX + 1, 0), None);
    assert_eq!(h.flat_index(NX, 0), None);
}

#[test]
fn degenerate_grid_geometry_is_rejected() {
    assert!(Grid::new(0, NY).is_err());
    assert!(Grid::new(NX, 0).is_err());
    assert!(Grid::new(NX, NY).is_ok());
}
