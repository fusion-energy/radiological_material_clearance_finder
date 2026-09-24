//! The crate's behaviour as yamc and yani see it, without the Python layer.
//!
//! The pytest suite exercises the same logic through the bindings. These cover
//! the regulatory rules that matter most, so a Rust-only change is checked by
//! `cargo test` alone.

use std::collections::BTreeMap;

use radiological_material_clearance_finder::{
    clearable_routes, clearance_index, clearance_indices, get_limit_set, limit_sets,
    nrc_waste_class, time_to_clear, uk_waste_category, ActivityUnit, ClearanceOptions, DecayData,
    Error, LimitSet, LimitSetData, Material, NrcWasteClass, UkCategory, BECQUEREL_PER_CURIE,
};

fn close(a: f64, b: f64, rel: f64) -> bool {
    (a - b).abs() <= rel * a.abs().max(b.abs())
}

fn limits(pairs: &[(&str, f64)]) -> BTreeMap<String, f64> {
    pairs.iter().map(|&(k, v)| (k.to_string(), v)).collect()
}

fn daughters(pairs: &[(&str, &[&str])]) -> BTreeMap<String, Vec<String>> {
    pairs
        .iter()
        .map(|&(p, ds)| (p.to_string(), ds.iter().map(|d| d.to_string()).collect()))
        .collect()
}

fn assay(pairs: &[(&str, f64)]) -> Material {
    Material::from_specific_activities(pairs.iter().copied()).unwrap()
}

fn assess(
    material: &Material,
    set: &LimitSet,
) -> radiological_material_clearance_finder::ClearanceResult {
    clearance_index(material, set, ClearanceOptions::default()).unwrap()
}

fn toy(name: &str, pairs: &[(&str, f64)]) -> LimitSetData {
    LimitSetData::new(name, "toy", "Bq/g", limits(pairs))
}

#[test]
fn the_readme_example_is_what_the_python_package_gives() {
    let steel =
        Material::from_atom_counts([("Fe56", 8.4e22), ("Co60", 2.1e9), ("Cs137", 3.1e8)]).unwrap();
    let result = assess(&steel, &get_limit_set("UK_EPR16_out_of_scope").unwrap());
    assert_eq!(result.index, 11.244691988476553);
    assert!(!result.clearable());
    assert!(result.uncovered.is_empty());
    assert_eq!(result.dominant(1)[0].0, "Co60");
}

#[test]
fn the_shipped_sets_load_with_the_sizes_the_sources_hold() {
    let names = limit_sets(None);
    assert_eq!(names.len(), 22);
    for (name, count) in [
        ("UK_EPR16_out_of_scope", 282),
        ("UK_IRR17_notification", 300),
        ("StrlSchV_unrestricted", 763),
        ("EU_BSS_clearance", 291),
        ("IAEA_GSR3_clearance", 257),
        ("Fetter", 81),
    ] {
        assert_eq!(get_limit_set(name).unwrap().limits().len(), count, "{name}");
    }
    assert_eq!(limit_sets(Some("uk")).len(), 6);
    match get_limit_set("NOT_A_SET") {
        Err(Error::UnknownLimitSet(m)) => assert!(m.contains("Available")),
        other => panic!("expected UnknownLimitSet, got {other:?}"),
    }
}

#[test]
fn pure_cobalt_60_has_the_textbook_specific_activity() {
    let material = Material::from_masses([("Co60", 1.0)]).unwrap();
    let curies_per_gram = material.specific_activity().unwrap() / BECQUEREL_PER_CURIE;
    assert!(close(curies_per_gram, 1131.0, 1e-3));
}

#[test]
fn atom_densities_fix_the_density_and_unit_conversions_agree() {
    let material = Material::from_atom_densities([("Fe56", 0.0849), ("Co60", 1e-10)])
        .unwrap()
        .with_volume(1000.0)
        .unwrap();
    let density = material.mass_density().unwrap();
    assert!(close(density, 7.89, 1e-2));
    let per_gram = material.activity(ActivityUnit::BqPerG).unwrap();
    assert!(close(
        material.activity(ActivityUnit::Bq).unwrap(),
        per_gram * material.mass().unwrap(),
        1e-12
    ));
    assert!(close(
        material.activity(ActivityUnit::CiPerM3).unwrap(),
        per_gram * density * 1e6 / BECQUEREL_PER_CURIE,
        1e-12
    ));
    assert!(
        material.clone().with_density(7.87).is_err(),
        "a density could only contradict it"
    );
}

#[test]
fn volumetric_activity_without_a_density_is_insufficient_data() {
    let material = Material::from_atom_counts([("Fe56", 1e22), ("Co60", 1e12)]).unwrap();
    assert!(matches!(
        material.activity(ActivityUnit::CiPerM3),
        Err(Error::InsufficientData(_))
    ));
    // And clearance_indices skips such sets rather than failing.
    let names: Vec<String> = clearance_indices(&material, None, ClearanceOptions::default())
        .unwrap()
        .into_iter()
        .map(|r| r.limit_set)
        .collect();
    assert!(names.iter().any(|n| n == "StrlSchV_unrestricted"));
    assert!(!names.iter().any(|n| n == "Fetter"));
}

#[test]
fn bad_inputs_are_rejected() {
    assert!(Material::from_atom_counts([("Co60", -1.0)]).is_err());
    assert!(Material::from_atom_counts([("Co60", f64::NAN)]).is_err());
    assert!(Material::from_atom_counts([("Fe", 1.0)]).is_err());
    assert!(assay(&[("Co60", 1.0)]).with_density(0.0).is_err());
    assert!(LimitSet::new(toy("BAD", &[("Co60", 0.0)])).is_err());
    assert!(LimitSet::new(toy("BAD", &[("Co60", f64::INFINITY)])).is_err());
    assert!(LimitSet::new(toy("BAD", &[("U240", 1000.0), ("U-240", 0.1)])).is_err());
    let mut units = toy("BAD", &[("Co60", 1.0)]);
    units.units = "sieverts".into();
    assert!(LimitSet::new(units).is_err());
}

#[test]
fn names_are_normalised_and_duplicates_summed() {
    let material = Material::from_atom_counts([("Co-60", 1e12), ("co60", 1e12)]).unwrap();
    assert_eq!(material.nuclides().collect::<Vec<_>>(), ["Co60"]);
    assert_eq!(material.atoms().unwrap()["Co60"], 2e12);
    assert!(
        Material::from_atom_counts([("Co60", 1e12), ("Cs137", 1e12)])
            .unwrap()
            .looks_truncated()
    );
    assert!(
        !Material::from_atom_counts([("Fe56", 1e22), ("Co60", 1e12)])
            .unwrap()
            .looks_truncated()
    );
}

#[test]
fn a_catch_all_limit_covers_unlisted_nuclides_unless_switched_off() {
    let mut data = toy("TOY_DEFAULT", &[("Co60", 1.0)]);
    data.default_limit = Some(0.01);
    let set = LimitSet::new(data).unwrap();
    let material = assay(&[("Co60", 1.0), ("Fe55", 0.02)]);
    let result = assess(&material, &set);
    assert!(close(result.index, 3.0, 1e-12));
    assert_eq!(result.defaulted, ["Fe55"]);

    let off = ClearanceOptions {
        apply_default_limit: false,
        ..Default::default()
    };
    let result = clearance_index(&material, &set, off).unwrap();
    assert!(close(result.index, 1.0, 1e-12));
    assert_eq!(result.uncovered, [("Fe55".to_string(), 0.02)]);
}

#[test]
fn a_parent_only_accounts_for_the_daughter_activity_it_supports() {
    let mut data = toy("TOY_CREDIT", &[("Sr90", 1.0), ("Y90", 1.0)]);
    data.secular_equilibrium = daughters(&[("Sr90", &["Y90"])]);
    let set = LimitSet::new(data).unwrap();

    let equilibrium = assess(&assay(&[("Sr90", 100.0), ("Y90", 100.0)]), &set);
    assert_eq!(equilibrium.excluded.len(), 1);
    assert!(close(equilibrium.index, 100.0, 1e-12));

    let lopsided = assess(&assay(&[("Sr90", 1.0), ("Y90", 1000.0)]), &set);
    assert!(lopsided.excluded.is_empty());
    assert_eq!(lopsided.credited, [("Y90".to_string(), 1.0)]);
    assert!(close(lopsided.index, 1000.0, 1e-12));
    assert!(close(lopsided.assessed_activity("Y90"), 999.0, 1e-12));

    let off = ClearanceOptions {
        exclude_daughters: false,
        ..Default::default()
    };
    let doubled = clearance_index(&assay(&[("Sr90", 1.0), ("Y90", 1.0)]), &set, off).unwrap();
    assert!(close(doubled.index, 2.0, 1e-12));
}

#[test]
fn a_parent_with_no_limit_cannot_account_for_its_daughters() {
    let material = assay(&[("Np237", 1e4), ("Pa233", 1e4)]);
    let result = assess(&material, &get_limit_set("StrlSchV_soil").unwrap());
    assert!(result.excluded.is_empty());
    assert!(!result.clearable());
    assert!(result.uncovered.iter().any(|(n, _)| n == "Np237"));
}

#[test]
fn a_lenient_equilibrium_limit_needs_real_equilibrium() {
    let set = get_limit_set("UK_IRR17_notification").unwrap();
    let trace = assess(&assay(&[("U240", 1.0), ("Np240_m1", 1e-12)]), &set);
    assert!(close(trace.index, 100.0, 1e-9));
    let equilibrium = assess(
        &assay(&[("U240", 1.0), ("Np240_m1", 1.0), ("Np240", 1.0)]),
        &set,
    );
    assert_eq!(equilibrium.limit_used("U240"), Some(100.0));
    assert!(equilibrium.clearable());
}

#[test]
fn a_sec_only_limit_is_applied_rather_than_left_unreachable() {
    let result = assess(
        &assay(&[("U238", 1e6)]),
        &get_limit_set("UK_IRR17_natural").unwrap(),
    );
    assert_eq!(result.limit_used("U238"), Some(1.0));
    assert!(result.uncovered.is_empty());
}

#[test]
fn the_short_half_life_scope_rule_is_a_whole_material_test() {
    let mut data = toy("TOY_SCOPE", &[("N16", 1.0), ("Co60", 1.0), ("Fe56", 1.0)]);
    data.min_half_life_scope = Some(100.0);
    let set = LimitSet::new(data).unwrap();
    assert!(assess(&assay(&[("N16", 1e6)]), &set).clearable());
    let mixed = assess(&assay(&[("N16", 1e6), ("Co60", 1e6)]), &set);
    assert!(!mixed.out_of_scope);
    assert!(close(mixed.index, 2e6, 1e-12));
    assert!(
        !assess(&assay(&[("Fe56", 1e6)]), &set).out_of_scope,
        "no radionuclides at all"
    );
}

#[test]
fn metal_and_dynamic_rows() {
    let mut data = toy("TOY_METAL", &[("C14", 8.0)]);
    data.metal_overrides = limits(&[("C14", 80.0), ("Ni59", 220.0)]);
    let set = LimitSet::new(data).unwrap();
    let material = assay(&[("C14", 8.0), ("Ni59", 220.0)]);
    assert!(close(assess(&material, &set).index, 1.0, 1e-12));
    let metal = ClearanceOptions {
        metal: true,
        ..Default::default()
    };
    assert!(close(
        clearance_index(&material, &set, metal).unwrap().index,
        1.1,
        1e-12
    ));

    let mut data = toy("TOY_DYNAMIC", &[]);
    data.dynamic_rule = Some("nrc_short_lived_class_a".into());
    let set = LimitSet::new(data).unwrap();
    let result = assess(&assay(&[("Mn54", 700.0), ("Cs137", 1e6)]), &set);
    assert_eq!(result.limits_used, [("Mn54".to_string(), 700.0)]);
    assert!(close(result.index, 1.0, 1e-12));
}

#[test]
fn nrc_class_rises_with_activity() {
    let bq_per_g = |ci_per_m3: f64| ci_per_m3 * 3.7e10 / (7.87 * 1e6);
    let classes: Vec<NrcWasteClass> = [0.5, 10.0, 1000.0, 1e5]
        .into_iter()
        .map(|target| {
            let material = assay(&[("Cs137", bq_per_g(target))])
                .with_density(7.87)
                .unwrap();
            nrc_waste_class(&material, false).unwrap()
        })
        .collect();
    assert_eq!(
        classes,
        [
            NrcWasteClass::A,
            NrcWasteClass::B,
            NrcWasteClass::C,
            NrcWasteClass::Gtcc
        ]
    );
}

#[test]
fn uk_categories() {
    assert_eq!(
        uk_waste_category(&assay(&[("H3", 30.0)])).unwrap().category,
        UkCategory::Vllw
    );
    assert_eq!(
        uk_waste_category(&assay(&[("H3", 50.0)])).unwrap().category,
        UkCategory::Llw
    );
    assert_eq!(
        uk_waste_category(&assay(&[("Am241", 5e3)]))
            .unwrap()
            .category,
        UkCategory::Ilw
    );
}

#[test]
fn time_to_clear_interpolates_and_refuses_ingrowth() {
    let set = LimitSet::new(toy("TOY_COOL", &[("Co60", 1.0)])).unwrap();
    let at = |index: f64| assay(&[("Co60", index)]);
    let materials = [at(10.0), at(0.1)];
    let series = [(0.0, &materials[0]), (100.0, &materials[1])];
    let found = time_to_clear(series, &set, ClearanceOptions::default(), false).unwrap();
    // Log interpolation puts an index of 1 exactly half way between 10 and 0.1.
    assert!(close(found.unwrap(), 50.0, 1e-12));

    let rebound = [at(10.0), at(0.1), at(5.0)];
    let series = [
        (0.0, &rebound[0]),
        (100.0, &rebound[1]),
        (200.0, &rebound[2]),
    ];
    assert!(matches!(
        time_to_clear(series, &set, ClearanceOptions::default(), false),
        Err(Error::Ingrowth(_))
    ));
    assert!(
        time_to_clear(series, &set, ClearanceOptions::default(), true)
            .unwrap()
            .is_some()
    );
    assert!(time_to_clear(
        [(0.0, &rebound[0])],
        &set,
        ClearanceOptions::default(),
        false
    )
    .is_err());
}

#[test]
fn clearable_routes_are_ordered_by_margin() {
    let material = Material::from_atom_counts([("Fe56", 1e22), ("Co60", 1e6)]).unwrap();
    let routes = clearable_routes(&material, None, ClearanceOptions::default()).unwrap();
    assert!(!routes.is_empty());
    let margins: Vec<f64> = routes
        .iter()
        .map(|n| {
            let set = get_limit_set(n).unwrap();
            assess(&material, &set).index / set.threshold()
        })
        .collect();
    assert!(margins.windows(2).all(|w| w[0] <= w[1]));
}

#[test]
fn decay_data_from_a_chain_and_overrides() {
    let data = DecayData::from_chain_xml_str(
        r#"<depletion_chain><nuclide name="Co60" half_life="1.0"/><nuclide name="Cs137"/></depletion_chain>"#,
        "chain",
    )
    .unwrap();
    assert_eq!(data.half_life("Co60").unwrap(), Some(1.0));
    assert_eq!(data.half_life("Cs137").unwrap(), None);
    for (body, message) in [
        (
            r#"<nuclide name="Co60" half_life="-1"/>"#,
            "positive and finite",
        ),
        (r#"<nuclide name="Co60" half_life="abc"/>"#, "not a number"),
        ("", "no readable nuclide entries"),
    ] {
        let err = DecayData::from_chain_xml_str(&format!("<c>{body}</c>"), "chain").unwrap_err();
        assert!(err.to_string().contains(message), "{err}");
    }
    assert!(DecayData::from_chain_xml_str("<c><unclosed>", "chain").is_err());

    let stable = DecayData::shared_default()
        .with_overrides([("Co-60", None)])
        .unwrap();
    assert_eq!(stable.half_life("Co60").unwrap(), None);
    assert!(matches!(
        DecayData::shared_default().half_life("Og296"),
        Err(Error::UnknownNuclide { .. })
    ));
}
