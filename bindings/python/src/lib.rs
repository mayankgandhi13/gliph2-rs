//! Python bindings for gliph2-rs (PyO3).
//!
//! The FFI surface is deliberately thin: lists of strings and one params dict
//! in, a dict of equal-length columns out. Defaults, validation and the
//! friendly signature live in the Python wrapper.

use std::collections::HashMap;

use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

use gliph2::convergence::{gliph2 as run_gliph2, Association, ConvergenceGroup, Gliph2Params, HlaTable};
use gliph2::io::{Repertoire, TcrRecord};

struct Params<'a, 'py>(&'a Bound<'py, PyDict>);

impl<'py> Params<'_, 'py> {
    fn get<T>(&self, key: &str) -> Option<T>
    where
        T: for<'a> FromPyObject<'a, 'py>,
    {
        self.0.get_item(key).ok().flatten().filter(|v| !v.is_none()).and_then(|v| v.extract().ok())
    }

    fn f64(&self, key: &str, default: f64) -> f64 {
        self.get(key).unwrap_or(default)
    }

    fn usize(&self, key: &str, default: usize) -> usize {
        self.get(key).unwrap_or(default)
    }

    fn bool(&self, key: &str, default: bool) -> bool {
        self.get(key).unwrap_or(default)
    }

    fn string(&self, key: &str) -> Option<String> {
        self.get::<String>(key).filter(|s| !s.is_empty())
    }
}

fn build_params(p: &Params) -> Gliph2Params {
    let mut g = Gliph2Params {
        lcminp: p.f64("lcminp", 0.01),
        lcminove: p.get("lcminove").unwrap_or_else(|| vec![1000.0, 100.0, 10.0]),
        kmer_mindepth: p.get("kmer_mindepth").unwrap_or(3),
        motif_distance_cutoff: p.usize("motif_distance_cutoff", 3),
        boundary_size: p.usize("boundary_size", 3),
        motif_length: p.get("motif_length").unwrap_or_else(|| vec![2, 3, 4]),
        discontinuous: p.bool("discontinuous_motifs", false),
        min_seq_length: p.usize("min_seq_length", 0),
        require_cf: p.bool("accept_sequences_with_c_f_start_end", true),
        all_aa_interchangeable: p.bool("all_aa_interchangeable", false),
        global_vgene: p.bool("global_vgene", false),
        cluster_min_size: p.usize("cluster_min_size", 2),
        local: p.bool("local_similarities", true),
        global: p.bool("global_similarities", true),
        global_max_p: p.get("global_max_p"),
        global_max_q: p.get("global_max_q"),
        min_subjects: p.usize("min_subjects", 0),
        association: None,
        hla: None,
    };
    if let (Some(case), Some(control)) = (p.string("case"), p.string("control")) {
        g.association = Some(Association {
            case,
            control,
            min_donors: p.usize("assoc_min_donors", 10),
            max_p: p.get("assoc_max_p"),
            max_q: p.get("assoc_max_q"),
            permute_seed: p.get("assoc_permute_seed"),
        });
    }
    // HLA typing arrives as two parallel lists.
    let donors: Vec<String> = p.get("hla_donor").unwrap_or_default();
    let alleles: Vec<String> = p.get("hla_allele").unwrap_or_default();
    if !donors.is_empty() {
        let mut by_donor: HashMap<String, Vec<String>> = HashMap::new();
        for (d, a) in donors.iter().zip(&alleles) {
            if !d.is_empty() && !a.is_empty() {
                by_donor.entry(d.clone()).or_default().push(a.clone());
            }
        }
        g.hla = Some(HlaTable { by_donor });
    }
    g
}

fn repertoire(cdr3: Vec<String>, v_gene: Vec<String>, subject: Vec<String>, count: Vec<u32>) -> Repertoire {
    let pick = |v: &[String], i: usize| v.get(i).filter(|s| !s.is_empty()).cloned();
    Repertoire::from_records(cdr3.iter().enumerate().filter_map(|(i, c)| {
        gliph2::io::normalize_cdr3(c.as_bytes()).map(|cdr3| TcrRecord {
            cdr3: cdr3.into_owned(),
            v_gene: pick(&v_gene, i),
            j_gene: None,
            subject: pick(&subject, i),
            count: count.get(i).copied().filter(|v| *v >= 1).unwrap_or(1),
        })
    }))
}

/// Runs `f` on a rayon pool of `n` threads (0 = rayon's default). A local
/// pool, rather than the global one, so repeated calls keep working.
fn with_threads<T: Send>(n: usize, f: impl FnOnce() -> T + Send) -> T {
    if n == 0 {
        return f();
    }
    match rayon::ThreadPoolBuilder::new().num_threads(n).build() {
        Ok(pool) => pool.install(f),
        Err(_) => f(),
    }
}

fn column<'py, T, F>(py: Python<'py>, groups: &[ConvergenceGroup], f: F) -> PyResult<Bound<'py, PyList>>
where
    T: IntoPyObject<'py>,
    F: Fn(&ConvergenceGroup) -> T,
{
    PyList::new(py, groups.iter().map(f))
}

/// Internal entry point behind `gliph2()`. Input is validated in Python, so
/// empty input simply yields empty columns.
#[pyfunction]
#[pyo3(signature = (cdr3, v_gene, subject, count, reference, params))]
fn run<'py>(
    py: Python<'py>,
    cdr3: Vec<String>,
    v_gene: Vec<String>,
    subject: Vec<String>,
    count: Vec<u32>,
    reference: Vec<String>,
    params: &Bound<'py, PyDict>,
) -> PyResult<Bound<'py, PyDict>> {
    let p = Params(params);
    let gp = build_params(&p);
    let threads = p.usize("threads", 0);

    let rep = repertoire(cdr3, v_gene, subject, count);
    let reference: Vec<Vec<u8>> =
        reference.iter().filter_map(|s| gliph2::io::normalize_cdr3(s.as_bytes()).map(|c| c.into_owned())).collect();

    let res = py.detach(|| with_threads(threads, || run_gliph2(&rep, &reference, &gp)));
    let g = &res.groups;

    let groups = PyDict::new(py);
    groups.set_item("type", column(py, g, |x| x.kind.as_str())?)?;
    groups.set_item("tag", column(py, g, |x| x.tag(&rep))?)?;
    groups.set_item("cluster_size", column(py, g, |x| x.cluster_size)?)?;
    groups.set_item("unique_cdr3_sample", column(py, g, |x| x.unique_cdr3_sample)?)?;
    groups.set_item("unique_cdr3_ref", column(py, g, |x| x.unique_cdr3_ref)?)?;
    groups.set_item("OvE", column(py, g, |x| x.ove)?)?;
    groups.set_item("fisher_score", column(py, g, |x| x.fisher_score)?)?;
    groups.set_item("fdr_q", column(py, g, |x| x.q_value)?)?;
    groups.set_item("n_subjects", column(py, g, |x| x.n_subjects)?)?;
    if res.association.is_some() {
        groups.set_item("case_donors", column(py, g, |x| x.case_donors)?)?;
        groups.set_item("control_donors", column(py, g, |x| x.control_donors)?)?;
        groups.set_item("assoc_p", column(py, g, |x| x.assoc_p)?)?;
        groups.set_item("assoc_q", column(py, g, |x| x.assoc_q)?)?;
    }
    if res.hla.is_some() {
        let names = &res.hla_alleles;
        groups.set_item("hla_allele", column(py, g, |x| names.get(x.hla_allele as usize).map(String::as_str))?)?;
        groups.set_item("hla_donors", column(py, g, |x| x.hla_donors)?)?;
        groups.set_item("hla_with_allele", column(py, g, |x| x.hla_with_allele)?)?;
        groups.set_item("hla_p", column(py, g, |x| x.hla_p)?)?;
        groups.set_item("hla_q", column(py, g, |x| x.hla_q)?)?;
    }
    // Members as a list of CDR3 lists, sorted as the engine emits them.
    let members = PyList::new(
        py,
        g.iter()
            .map(|x| {
                PyList::new(py, x.members.iter().map(|&m| String::from_utf8_lossy(&rep.sequences[m as usize]).into_owned()))
            })
            .collect::<PyResult<Vec<_>>>()?,
    )?;
    groups.set_item("members", members)?;

    let m = &res.selected_motifs;
    let motifs = PyDict::new(py);
    motifs.set_item("motif", PyList::new(py, m.iter().map(|x| x.motif.to_string().replace('%', ".")))?)?;
    motifs.set_item("num_in_sample", PyList::new(py, m.iter().map(|x| x.num_in_sample))?)?;
    motifs.set_item("num_in_ref", PyList::new(py, m.iter().map(|x| x.num_in_ref))?)?;
    motifs.set_item("fisher_score", PyList::new(py, m.iter().map(|x| x.fisher_score))?)?;
    motifs.set_item("num_fold", PyList::new(py, m.iter().map(|x| x.num_fold))?)?;

    let out = PyDict::new(py);
    out.set_item("groups", groups)?;
    out.set_item("motifs", motifs)?;
    out.set_item("n_sample", res.n_sample)?;
    out.set_item("n_reference", res.n_reference)?;
    out.set_item("n_input_cdr3", rep.len())?;
    match res.association {
        Some(s) => {
            let d = PyDict::new(py);
            d.set_item("n_case", s.n_case)?;
            d.set_item("n_control", s.n_control)?;
            d.set_item("tested", s.tested)?;
            d.set_item("significant_q05", s.significant_q05)?;
            out.set_item("association", d)?;
        }
        None => out.set_item("association", py.None())?,
    }
    match res.hla {
        Some(s) => {
            let d = PyDict::new(py);
            d.set_item("typed_donors", s.typed_donors)?;
            d.set_item("scored", s.scored)?;
            d.set_item("tests", s.tests)?;
            d.set_item("significant_q05", s.significant_q05)?;
            out.set_item("hla", d)?;
        }
        None => out.set_item("hla", py.None())?,
    }
    Ok(out)
}

#[pymodule]
fn _gliph2rs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(run, m)?)?;
    Ok(())
}
