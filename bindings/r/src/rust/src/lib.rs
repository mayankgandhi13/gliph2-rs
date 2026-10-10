//! R bindings for gliph2-rs (extendr).
//!
//! The FFI surface is deliberately thin: character vectors in, named lists of
//! equal-length columns out, with every option passed in one `params` list so
//! the R side owns the defaults and the data-frame assembly.

use std::collections::HashMap;

use extendr_api::prelude::*;
use gliph2::convergence::{gliph2 as run_gliph2, Association, Gliph2Params};
use gliph2::io::{Repertoire, TcrRecord};

/// Runs `f` on a rayon pool of `n` threads (0 = rayon's default). A local
/// pool, rather than the global one, so repeated calls from one R session
/// keep working.
fn with_threads<T: Send>(n: usize, f: impl FnOnce() -> T + Send) -> T {
    if n == 0 {
        return f();
    }
    match rayon::ThreadPoolBuilder::new().num_threads(n).build() {
        Ok(pool) => pool.install(f),
        Err(_) => f(),
    }
}

struct Params<'a>(HashMap<&'a str, Robj>);

impl<'a> Params<'a> {
    fn new(list: &'a List) -> Self {
        Params(list.iter().filter(|(_, v)| !v.is_null()).collect())
    }

    fn f64(&self, name: &str, default: f64) -> f64 {
        self.0.get(name).and_then(|r| r.as_real()).unwrap_or(default)
    }

    fn opt_f64(&self, name: &str) -> Option<f64> {
        self.0.get(name).and_then(|r| r.as_real())
    }

    fn usize(&self, name: &str, default: usize) -> usize {
        self.0.get(name).and_then(|r| r.as_real()).map_or(default, |v| v.max(0.0) as usize)
    }

    fn opt_u64(&self, name: &str) -> Option<u64> {
        self.0.get(name).and_then(|r| r.as_real()).map(|v| v.max(0.0) as u64)
    }

    fn bool(&self, name: &str, default: bool) -> bool {
        self.0.get(name).and_then(|r| r.as_bool()).unwrap_or(default)
    }

    fn str(&self, name: &str) -> Option<String> {
        self.0.get(name).and_then(|r| r.as_str()).map(str::to_owned).filter(|s| !s.is_empty())
    }

    fn f64_vec(&self, name: &str, default: Vec<f64>) -> Vec<f64> {
        self.0.get(name).and_then(|r| r.as_real_vector()).filter(|v| !v.is_empty()).unwrap_or(default)
    }
}

fn build_params(p: &Params) -> Gliph2Params {
    let mut g = Gliph2Params {
        lcminp: p.f64("lcminp", 0.01),
        lcminove: p.f64_vec("lcminove", vec![1000.0, 100.0, 10.0]),
        kmer_mindepth: p.opt_u64("kmer_mindepth").unwrap_or(3),
        motif_distance_cutoff: p.usize("motif_distance_cutoff", 3),
        boundary_size: p.usize("boundary_size", 3),
        motif_length: p.f64_vec("motif_length", vec![2.0, 3.0, 4.0]).iter().map(|v| *v as usize).collect(),
        discontinuous: p.bool("discontinuous_motifs", false),
        min_seq_length: p.usize("min_seq_length", 0),
        require_cf: p.bool("accept_sequences_with_C_F_start_end", true),
        all_aa_interchangeable: p.bool("all_aa_interchangeable", false),
        global_vgene: p.bool("global_vgene", false),
        cluster_min_size: p.usize("cluster_min_size", 2),
        local: p.bool("local_similarities", true),
        global: p.bool("global_similarities", true),
        global_max_p: p.opt_f64("global_max_p"),
        global_max_q: p.opt_f64("global_max_q"),
        min_subjects: p.usize("min_subjects", 0),
        association: None,
    };
    if let (Some(case), Some(control)) = (p.str("case"), p.str("control")) {
        g.association = Some(Association {
            case,
            control,
            min_donors: p.usize("assoc_min_donors", 10),
            max_p: p.opt_f64("assoc_max_p"),
            max_q: p.opt_f64("assoc_max_q"),
            permute_seed: p.opt_u64("assoc_permute_seed"),
        });
    }
    g
}

/// Builds a named list. Every call site passes equal-length slices, so the
/// error case is unreachable; extendr 0.9 turns a returned `Err` into a Rust
/// panic across the FFI boundary, so reachable validation lives in R.
fn columns(names: &[&str], values: Vec<Robj>) -> List {
    debug_assert_eq!(names.len(), values.len());
    List::from_names_and_values(names, values).expect("column names and values have equal length")
}

/// Builds a repertoire from parallel character vectors; `""` means missing.
fn repertoire(cdr3: Vec<String>, v_gene: Vec<String>, subject: Vec<String>, count: Vec<f64>) -> Repertoire {
    let pick = |v: &[String], i: usize| v.get(i).filter(|s| !s.is_empty()).cloned();
    Repertoire::from_records(cdr3.iter().enumerate().filter_map(|(i, c)| {
        gliph2::io::normalize_cdr3(c.as_bytes()).map(|cdr3| TcrRecord {
            cdr3: cdr3.into_owned(),
            v_gene: pick(&v_gene, i),
            j_gene: None,
            subject: pick(&subject, i),
            count: count.get(i).copied().filter(|v| *v >= 1.0).map_or(1, |v| v as u32),
        })
    }))
}

/// Internal entry point behind `gliph2()`. Returns a list of column lists.
/// Input is validated in R, so this never fails: empty input yields empty
/// columns.
/// @export
#[extendr]
fn gliph2_impl(
    cdr3: Vec<String>,
    v_gene: Vec<String>,
    subject: Vec<String>,
    count: Vec<f64>,
    reference: Vec<String>,
    params: List,
) -> List {
    let p = Params::new(&params);
    let gp = build_params(&p);
    let n_threads = p.usize("n_cores", 0);

    let rep = repertoire(cdr3, v_gene, subject, count);
    let reference: Vec<Vec<u8>> =
        reference.iter().filter_map(|s| gliph2::io::normalize_cdr3(s.as_bytes()).map(|c| c.into_owned())).collect();

    let res = with_threads(n_threads, || run_gliph2(&rep, &reference, &gp));
    let g = &res.groups;
    let text = |f: &dyn Fn(&gliph2::convergence::ConvergenceGroup) -> String| -> Robj {
        g.iter().map(f).collect::<Vec<String>>().into()
    };
    let int = |f: &dyn Fn(&gliph2::convergence::ConvergenceGroup) -> i32| -> Robj {
        g.iter().map(f).collect::<Vec<i32>>().into()
    };
    let real = |f: &dyn Fn(&gliph2::convergence::ConvergenceGroup) -> f64| -> Robj {
        g.iter().map(f).collect::<Vec<f64>>().into()
    };

    let mut names = vec![
        "type",
        "tag",
        "cluster_size",
        "unique_cdr3_sample",
        "unique_cdr3_ref",
        "OvE",
        "fisher.score",
        "fdr.q",
        "n_subjects",
    ];
    let mut values = vec![
        text(&|x| x.kind.as_str().to_owned()),
        text(&|x| x.tag(&rep)),
        int(&|x| x.cluster_size as i32),
        int(&|x| x.unique_cdr3_sample as i32),
        real(&|x| x.unique_cdr3_ref as f64),
        real(&|x| x.ove),
        real(&|x| x.fisher_score),
        real(&|x| x.q_value),
        int(&|x| x.n_subjects as i32),
    ];
    if res.association.is_some() {
        names.extend(["case_donors", "control_donors", "assoc.p", "assoc.q"]);
        values.push(int(&|x| x.case_donors as i32));
        values.push(int(&|x| x.control_donors as i32));
        values.push(real(&|x| x.assoc_p));
        values.push(real(&|x| x.assoc_q));
    }
    // Members as space-separated CDR3s, as turboGliph writes them.
    names.push("members");
    values.push(text(&|x| {
        x.members
            .iter()
            .map(|&m| String::from_utf8_lossy(&rep.sequences[m as usize]).into_owned())
            .collect::<Vec<String>>()
            .join(" ")
    }));
    let groups = columns(&names, values);

    let m = &res.selected_motifs;
    let motifs = columns(
        &["motif", "num_in_sample", "num_in_ref", "fisher.score", "num_fold"],
        vec![
            m.iter().map(|x| x.motif.to_string().replace('%', ".")).collect::<Vec<String>>().into(),
            m.iter().map(|x| x.num_in_sample as f64).collect::<Vec<f64>>().into(),
            m.iter().map(|x| x.num_in_ref as f64).collect::<Vec<f64>>().into(),
            m.iter().map(|x| x.fisher_score).collect::<Vec<f64>>().into(),
            m.iter().map(|x| x.num_fold).collect::<Vec<f64>>().into(),
        ],
    );

    let assoc: Robj = match res.association {
        Some(s) => columns(
            &["n_case", "n_control", "tested", "significant_q05"],
            vec![
                (s.n_case as i32).into(),
                (s.n_control as i32).into(),
                (s.tested as f64).into(),
                (s.significant_q05 as f64).into(),
            ],
        )
        .into(),
        None => NULL.into(),
    };

    columns(
        &["groups", "motifs", "n_sample", "n_reference", "n_input_cdr3", "association"],
        vec![
            groups.into(),
            motifs.into(),
            (res.n_sample as f64).into(),
            (res.n_reference as f64).into(),
            (rep.len() as f64).into(),
            assoc,
        ],
    )
}

extendr_module! {
    mod gliph2rs;
    fn gliph2_impl;
}
