//! Serialized constraint-model round-trip (serde feature only).
#![cfg(feature = "serde")]

use rcad_gcs::prelude::*;

/// The serialized constraint model round-trips through JSON (serde feature).
#[test]
fn constraint_model_json_roundtrip() {
    let mut sk = Sketch::new();
    let l = sk.add_line(0.0, 0.0, 10.0, 0.0);
    sk.add_constraint(Constraint::Horizontal(l));
    sk.add_constraint(Constraint::LineLength {
        line: l,
        length: 10.0,
    });

    let json = serde_json::to_string(&sk).expect("serialize");
    let back: Sketch = serde_json::from_str(&json).expect("deserialize");
    assert_eq!(back.entities.len(), sk.entities.len());
    assert_eq!(back.constraints().len(), sk.constraints().len());

    // The round-tripped sketch still solves to the same geometry.
    let mut sk2 = back;
    let result = sk2.solve();
    assert!(result.converged, "residual {}", result.residual);
}
