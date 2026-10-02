use std::fs::File;
use std::io::{BufWriter, Write};
use std::process::ExitCode;
use std::time::Instant;

use gliph2::{io, synthetic, Gliph2Params, Params};

const USAGE: &str = "\
gliph2-rs — deterministic GLIPH2 core

USAGE:
  gliph2-rs cluster --input FILE [--reference FILE] [--out PREFIX]
                    [--threads N] [--mismatch-trim N] [--min-cluster N]
                    [--max-p P] [--fdr] [--gapped] [--min-count N] [--no-motifs]
  gliph2-rs gliph2  --input FILE --reference FILE [--out PREFIX] [--threads N]
                    [--paper-params] [--lcminp P] [--lcminove A[,B,C]]
                    [--kmer-mindepth N] [--min-seq-length N] [--motif-distance N]
                    [--all-aa-interchangeable] [--global-vgene] [--discontinuous]
                    [--min-cluster N] [--no-local] [--no-global]
  gliph2-rs synth   --n N --out FILE [--seed S] [--reference]
  gliph2-rs expand  --input FILE --factor F --out FILE [--seed S]

`cluster` writes PREFIX_clusters.tsv and PREFIX_motifs.tsv and prints a
one-line JSON timing summary to stdout. Inputs are TSV with a header
(GLIPH2, VDJtools or AIRR column names) or one CDR3 per line.

  --fdr     threshold Benjamini-Hochberg q-values instead of raw p-values
  --gapped  also test gapped motifs such as S%G and SL%Q

`gliph2` reports GLIPH2 convergence groups (local motif groups and global
struct groups, not merged) with turboGliph's gliph2() semantics and defaults,
writing PREFIX_groups.tsv. --paper-params uses the GLIPH2 parameter-file
values instead (lcminp 0.001, lcminove 10, min length 8, all aa
interchangeable).";

struct Args(Vec<String>);

impl Args {
    fn value(&self, flag: &str) -> Option<&str> {
        self.0.iter().position(|a| a == flag).and_then(|i| self.0.get(i + 1)).map(String::as_str)
    }
    fn flag(&self, flag: &str) -> bool {
        self.0.iter().any(|a| a == flag)
    }
    fn parse<T: std::str::FromStr>(&self, flag: &str) -> Result<Option<T>, String> {
        self.value(flag)
            .map(|v| v.parse().map_err(|_| format!("invalid value for {flag}: {v}")))
            .transpose()
    }
    fn require(&self, flag: &str) -> Result<&str, String> {
        self.value(flag).ok_or_else(|| format!("missing {flag}"))
    }
}

fn write_seqs(path: &str, seqs: &[Vec<u8>]) -> std::io::Result<()> {
    let mut w = BufWriter::new(File::create(path)?);
    writeln!(w, "CDR3b")?;
    for s in seqs {
        w.write_all(s)?;
        w.write_all(b"\n")?;
    }
    w.flush()
}

fn cluster(a: &Args) -> Result<(), String> {
    let wall = Instant::now();
    if let Some(n) = a.parse::<usize>("--threads")? {
        rayon::ThreadPoolBuilder::new().num_threads(n).build_global().map_err(|e| e.to_string())?;
    }
    let mut params = Params::default();
    if let Some(v) = a.parse("--mismatch-trim")? {
        params.local.mismatch_trim = v;
    }
    if let Some(v) = a.parse("--min-cluster")? {
        params.assembly.min_size = v;
    }
    if let Some(v) = a.parse("--max-p")? {
        params.enrichment.max_p = v;
    }
    if let Some(v) = a.parse("--min-count")? {
        params.enrichment.min_count = v;
    }
    params.skip_motifs = a.flag("--no-motifs");
    params.enrichment.fdr = a.flag("--fdr");
    params.motif.gapped = a.flag("--gapped");

    let input = a.require("--input")?;
    let now = Instant::now();
    let rep = io::read_repertoire(input).map_err(|e| format!("{input}: {e}"))?;
    let n_rows = rep.rows.len();
    let reference = match a.value("--reference") {
        Some(p) => Some(io::read_sequences(p).map_err(|e| format!("{p}: {e}"))?),
        None => None,
    };
    let t_io = now.elapsed();
    if rep.skipped > 0 {
        eprintln!("skipped {} rows with invalid CDR3s", rep.skipped);
    }
    if reference.is_none() && !params.skip_motifs {
        eprintln!("no --reference given; motif enrichment skipped");
    }

    let res = gliph2::run(&rep.sequences, reference.as_deref(), &params);

    let now = Instant::now();
    let prefix = a.value("--out").unwrap_or("gliph2rs");
    let mut w = BufWriter::new(File::create(format!("{prefix}_clusters.tsv")).map_err(|e| e.to_string())?);
    gliph2::write_clusters(&mut w, &rep, &res.clusters).map_err(|e| e.to_string())?;
    let mut w = BufWriter::new(File::create(format!("{prefix}_motifs.tsv")).map_err(|e| e.to_string())?);
    gliph2::write_motifs(&mut w, &res.motifs).map_err(|e| e.to_string())?;
    let t_out = now.elapsed();

    let t = &res.timings;
    println!(
        "{{\"rows\":{},\"unique\":{},\"reference\":{},\"threads\":{},\"local_edges\":{},\"motifs\":{},\"clusters\":{},\
\"t_io\":{:.4},\"t_local\":{:.4},\"t_ref_table\":{:.4},\"t_enrich\":{:.4},\"t_assembly\":{:.4},\"t_output\":{:.4},\"t_total\":{:.4}}}",
        n_rows,
        rep.len(),
        reference.as_ref().map_or(0, |r| r.len()),
        rayon::current_num_threads(),
        res.local_edges.len(),
        res.motifs.len(),
        res.clusters.len(),
        t_io.as_secs_f64(),
        t.local.as_secs_f64(),
        t.reference_table.as_secs_f64(),
        t.enrichment.as_secs_f64(),
        t.assembly.as_secs_f64(),
        t_out.as_secs_f64(),
        wall.elapsed().as_secs_f64()
    );
    Ok(())
}

fn gliph2_cmd(a: &Args) -> Result<(), String> {
    let wall = Instant::now();
    if let Some(n) = a.parse::<usize>("--threads")? {
        rayon::ThreadPoolBuilder::new().num_threads(n).build_global().map_err(|e| e.to_string())?;
    }
    let mut p = if a.flag("--paper-params") { Gliph2Params::gliph2_paper() } else { Gliph2Params::default() };
    if let Some(v) = a.parse("--lcminp")? {
        p.lcminp = v;
    }
    if let Some(v) = a.value("--lcminove") {
        p.lcminove = v.split(',').map(|x| x.parse().map_err(|_| format!("invalid --lcminove: {v}"))).collect::<Result<_, _>>()?;
    }
    if let Some(v) = a.parse("--kmer-mindepth")? {
        p.kmer_mindepth = v;
    }
    if let Some(v) = a.parse("--min-seq-length")? {
        p.min_seq_length = v;
    }
    if let Some(v) = a.parse("--motif-distance")? {
        p.motif_distance_cutoff = v;
    }
    if let Some(v) = a.parse("--min-cluster")? {
        p.cluster_min_size = v;
    }
    p.all_aa_interchangeable |= a.flag("--all-aa-interchangeable");
    p.global_vgene = a.flag("--global-vgene");
    p.discontinuous = a.flag("--discontinuous");
    p.local = !a.flag("--no-local");
    p.global = !a.flag("--no-global");

    let input = a.require("--input")?;
    let refp = a.require("--reference")?;
    let now = Instant::now();
    let rep = io::read_repertoire(input).map_err(|e| format!("{input}: {e}"))?;
    let reference = io::read_sequences(refp).map_err(|e| format!("{refp}: {e}"))?;
    let t_io = now.elapsed();

    let now = Instant::now();
    let res = gliph2::gliph2(&rep, &reference, &p);
    let t_run = now.elapsed();

    let prefix = a.value("--out").unwrap_or("gliph2rs");
    let mut w = BufWriter::new(File::create(format!("{prefix}_groups.tsv")).map_err(|e| e.to_string())?);
    gliph2::convergence::write_groups(&mut w, &rep, &res.groups).map_err(|e| e.to_string())?;
    w.flush().map_err(|e| e.to_string())?;
    let n_local = res.groups.iter().filter(|g| g.kind == gliph2::convergence::GroupType::Local).count();
    println!(
        "{{\"rows\":{},\"unique\":{},\"sample\":{},\"reference\":{},\"threads\":{},\"motifs\":{},\"local_groups\":{},\"global_groups\":{},\"t_io\":{:.4},\"t_gliph2\":{:.4},\"t_total\":{:.4}}}",
        rep.rows.len(),
        rep.len(),
        res.n_sample,
        res.n_reference,
        rayon::current_num_threads(),
        res.selected_motifs.len(),
        n_local,
        res.groups.len() - n_local,
        t_io.as_secs_f64(),
        t_run.as_secs_f64(),
        wall.elapsed().as_secs_f64()
    );
    Ok(())
}

fn main() -> ExitCode {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let Some(cmd) = argv.first().cloned() else {
        eprintln!("{USAGE}");
        return ExitCode::FAILURE;
    };
    let a = Args(argv);
    let result = match cmd.as_str() {
        "cluster" => cluster(&a),
        "gliph2" => gliph2_cmd(&a),
        "synth" => (|| {
            let n = a.parse("--n")?.ok_or("missing --n")?;
            let seed = a.parse("--seed")?.unwrap_or(1);
            let seqs = if a.flag("--reference") { synthetic::reference(n, seed) } else { synthetic::repertoire(n, seed) };
            write_seqs(a.require("--out")?, &seqs).map_err(|e| e.to_string())
        })(),
        "expand" => (|| {
            let input = a.require("--input")?;
            let base = io::read_sequences(input).map_err(|e| format!("{input}: {e}"))?;
            let factor = a.parse("--factor")?.ok_or("missing --factor")?;
            let seed = a.parse("--seed")?.unwrap_or(1);
            write_seqs(a.require("--out")?, &synthetic::expand(&base, factor, seed)).map_err(|e| e.to_string())
        })(),
        "-h" | "--help" | "help" => {
            println!("{USAGE}");
            Ok(())
        }
        other => Err(format!("unknown command: {other}\n\n{USAGE}")),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("error: {e}");
            ExitCode::FAILURE
        }
    }
}
