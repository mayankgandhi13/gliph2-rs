//! Union-find cluster assembly over local and global (motif) edges.
//!
//! Edges are applied in a fixed order and unions always keep the smaller id
//! as root, so the resulting partition and cluster numbering depend only on
//! the input — not on thread scheduling.

use crate::enrichment::MotifResult;
use crate::io::SequenceId;
use crate::local_cluster::LocalEdge;
use crate::motif::Kmer;

/// Disjoint-set forest with path compression; the root is always the
/// smallest member id.
#[derive(Debug, Clone)]
pub struct UnionFind {
    parent: Vec<SequenceId>,
}

impl UnionFind {
    pub fn new(n: usize) -> Self {
        UnionFind { parent: (0..n as SequenceId).collect() }
    }

    pub fn find(&mut self, mut x: SequenceId) -> SequenceId {
        let mut root = x;
        while self.parent[root as usize] != root {
            root = self.parent[root as usize];
        }
        while self.parent[x as usize] != root {
            let next = self.parent[x as usize];
            self.parent[x as usize] = root;
            x = next;
        }
        root
    }

    pub fn union(&mut self, a: SequenceId, b: SequenceId) -> bool {
        let (ra, rb) = (self.find(a), self.find(b));
        if ra == rb {
            return false;
        }
        let (lo, hi) = if ra < rb { (ra, rb) } else { (rb, ra) };
        self.parent[hi as usize] = lo;
        true
    }
}

#[derive(Debug, Clone, Default)]
pub struct Cluster {
    /// 0-based, assigned in order of each cluster's smallest member.
    pub id: usize,
    /// Sorted sequence ids.
    pub members: Vec<SequenceId>,
    pub n_local_edges: usize,
    /// Enriched motifs that contributed to this cluster, sorted.
    pub motifs: Vec<Kmer>,
}

#[derive(Debug, Clone, Copy)]
pub struct AssemblyParams {
    /// Drop clusters smaller than this (singletons are dropped by default).
    pub min_size: usize,
}

impl Default for AssemblyParams {
    fn default() -> Self {
        AssemblyParams { min_size: 2 }
    }
}

/// Merges sequences connected by a local edge or by sharing an enriched motif.
pub fn assemble_clusters(
    n_sequences: usize,
    local: &[LocalEdge],
    motifs: &[MotifResult],
    params: &AssemblyParams,
) -> Vec<Cluster> {
    let mut uf = UnionFind::new(n_sequences);
    for e in local {
        uf.union(e.a, e.b);
    }
    for m in motifs {
        if let Some((&first, rest)) = m.members.split_first() {
            for &x in rest {
                uf.union(first, x);
            }
        }
    }

    // Roots are the smallest member, so scanning ids in order yields clusters
    // already sorted by smallest member, each with sorted members.
    let mut slot: Vec<usize> = vec![usize::MAX; n_sequences];
    let mut clusters: Vec<Cluster> = Vec::new();
    for i in 0..n_sequences as SequenceId {
        let r = uf.find(i) as usize;
        if slot[r] == usize::MAX {
            slot[r] = clusters.len();
            clusters.push(Cluster::default());
        }
        clusters[slot[r]].members.push(i);
    }
    for e in local {
        clusters[slot[uf.find(e.a) as usize]].n_local_edges += 1;
    }
    for m in motifs {
        if let Some(&first) = m.members.first() {
            let c = &mut clusters[slot[uf.find(first) as usize]];
            if c.motifs.last() != Some(&m.motif) {
                c.motifs.push(m.motif);
            }
        }
    }

    let mut out: Vec<Cluster> = clusters
        .into_iter()
        .filter(|c| c.members.len() >= params.min_size.max(1))
        .collect();
    for (i, c) in out.iter_mut().enumerate() {
        c.id = i;
        c.motifs.sort_unstable();
        c.motifs.dedup();
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn edge(a: u32, b: u32) -> LocalEdge {
        LocalEdge { a, b, position: 0 }
    }

    #[test]
    fn merges_local_and_motif_edges() {
        let motif = MotifResult {
            motif: Kmer::new(b"SLG"),
            input_count: 2,
            reference_count: 0,
            fold: f64::INFINITY,
            p_value: 0.0,
            members: vec![2, 5],
        };
        let c = assemble_clusters(7, &[edge(0, 2), edge(3, 4)], &[motif], &AssemblyParams::default());
        let members: Vec<_> = c.iter().map(|c| c.members.clone()).collect();
        assert_eq!(members, vec![vec![0, 2, 5], vec![3, 4]]);
        assert_eq!(c[0].n_local_edges, 1);
        assert_eq!(c[0].motifs.len(), 1);
        assert_eq!(c[1].id, 1);
    }

    #[test]
    fn edge_order_does_not_matter() {
        let edges = [edge(5, 9), edge(1, 5), edge(0, 3), edge(3, 8), edge(2, 7)];
        let base = assemble_clusters(10, &edges, &[], &AssemblyParams::default());
        let mut rev = edges.to_vec();
        rev.reverse();
        let other = assemble_clusters(10, &rev, &[], &AssemblyParams::default());
        let m = |c: &[Cluster]| c.iter().map(|c| c.members.clone()).collect::<Vec<_>>();
        assert_eq!(m(&base), m(&other));
    }
}
