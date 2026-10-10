//! Deterministic synthetic repertoires for scaling benchmarks.
//!
//! Sequences look roughly like CDR3β: a germline-like N-terminal stem, a
//! random interior, and a J-like C-terminal stem. A fraction are planted as
//! near-copies of a seed sequence (to create local edges) or carry a planted
//! motif (to create enrichment), so benchmarks exercise every stage.

const AA_INTERIOR: &[u8] = b"AAGGGSSSSTTLLQQEEDDNNPPRRYYVVKKIIHHFFMWC";
const STEMS_N: &[&[u8]] = &[b"CASS", b"CASR", b"CAST", b"CSAR", b"CASG", b"CAWS"];
const STEMS_C: &[&[u8]] = &[b"YEQYF", b"TQYF", b"NEQFF", b"GYTF", b"TEAFF", b"QPQHF"];
const PLANTED: &[&[u8]] = &[b"RGWG", b"WDQY", b"KHMP", b"FWMK", b"HWNC"];

/// xorshift64* — tiny, fast and stable across platforms and versions.
#[derive(Debug, Clone)]
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Rng(seed.wrapping_mul(0x9E3779B97F4A7C15) | 1)
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545F4914F6CDD1D)
    }

    pub fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }

    pub fn chance(&mut self, p: f64) -> bool {
        ((self.next_u64() >> 11) as f64 / (1u64 << 53) as f64) < p
    }
}

fn random_cdr3(rng: &mut Rng) -> Vec<u8> {
    let mut s = STEMS_N[rng.below(STEMS_N.len())].to_vec();
    let interior = 3 + rng.below(7);
    for _ in 0..interior {
        s.push(AA_INTERIOR[rng.below(AA_INTERIOR.len())]);
    }
    s.extend_from_slice(STEMS_C[rng.below(STEMS_C.len())]);
    s
}

fn mutate_one(rng: &mut Rng, s: &[u8]) -> Vec<u8> {
    let mut v = s.to_vec();
    let p = 3 + rng.below(v.len().saturating_sub(6).max(1));
    let p = p.min(v.len() - 1);
    v[p] = AA_INTERIOR[rng.below(AA_INTERIOR.len())];
    v
}

/// A background-like reference repertoire (no planted structure).
pub fn reference(n: usize, seed: u64) -> Vec<Vec<u8>> {
    let mut rng = Rng::new(seed);
    (0..n).map(|_| random_cdr3(&mut rng)).collect()
}

/// An input repertoire with planted local families and enriched motifs.
pub fn repertoire(n: usize, seed: u64) -> Vec<Vec<u8>> {
    let mut rng = Rng::new(seed ^ 0xA5A5_5A5A);
    let mut out: Vec<Vec<u8>> = Vec::with_capacity(n);
    while out.len() < n {
        if !out.is_empty() && rng.chance(0.15) {
            let src = out[rng.below(out.len())].clone();
            out.push(mutate_one(&mut rng, &src));
        } else if rng.chance(0.02) {
            let mut s = random_cdr3(&mut rng);
            let m = PLANTED[rng.below(PLANTED.len())];
            let at = 4.min(s.len());
            s.splice(at..at, m.iter().copied());
            out.push(s);
        } else {
            out.push(random_cdr3(&mut rng));
        }
    }
    out
}

/// Expands a real repertoire `factor`-fold: originals first, then copies with
/// one random interior substitution each.
pub fn expand(base: &[Vec<u8>], factor: usize, seed: u64) -> Vec<Vec<u8>> {
    let mut rng = Rng::new(seed);
    let mut out = base.to_vec();
    if base.is_empty() {
        return out;
    }
    for _ in 1..factor.max(1) {
        for s in base {
            out.push(if s.len() > 6 { mutate_one(&mut rng, s) } else { s.clone() });
        }
    }
    out
}
