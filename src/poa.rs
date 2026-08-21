//! Stripped-down and refactored implemenation of the POA engine from crate 
//! `poa_consensus`, modified for our use case with a reusable Poa graph object.
//! 
//! https://github.com/Psy-Fer/poa-consensus
//! 
//! The open source license from the poa_consensus crate from August 2026 folows.
//! All code in this create is subject to the same MIT license.
//! 
//! MIT License
//! 
//! Copyright (c) 2026 James Ferguson
//! 
//! Permission is hereby granted, free of charge, to any person obtaining a copy
//! of this software and associated documentation files (the "Software"), to deal
//! in the Software without restriction, including without limitation the rights
//! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//! copies of the Software, and to permit persons to whom the Software is
//! furnished to do so, subject to the following conditions:
//! 
//! The above copyright notice and this permission notice shall be included in all
//! copies or substantial portions of the Software.
//! 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

/* -----------------------------------------------------------------------------
PoaConfig
----------------------------------------------------------------------------- */
/// POA alignment modes. Global always extends to the boundaries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AlignmentMode {
    Global,
    SemiGlobal,
}

/// Configuration parameters for a POA graph engine. These parameters are used
/// as provided for all iterative calls to an instantiated graph engine.
#[derive(Debug, Clone)]
pub struct PoaConfig {
    /// 0 = unbanded (full NW over DAG). Default 50.
    pub band_width: usize,

    /// Enable adaptive band width using the abPOA formula: w = b + f*L.
    pub adaptive_band: bool,
    /// Base component of the adaptive band formula (abPOA default: 10).
    pub adaptive_band_b: usize,
    /// Length-proportional component of the adaptive band formula (abPOA default: 0.01).
    pub adaptive_band_f: f32,

    /// Score for a base match. Positive. Default 2.
    pub match_score: i32,
    /// Penalty for a base mismatch. Negative. Default -4.
    pub mismatch_score: i32,
    /// One-time penalty when a gap opens. Negative. Default -4.
    pub gap_open: i32,
    /// Per-base penalty inside a gap. Negative. Default -3 (less than gap_open).
    pub gap_extend: i32,

    /// Number of sequences that must cover a boundary node for it to appear in 
    /// the reported consensus. Boundary nodes below this coverage are trimmed
    /// but internal low-coverage spans remain.
    pub min_boundary_coverage: u32,

    /// Global vs. semi-global alignment. Default AlignmentMode::SemiGlobal.
    pub alignment_mode: AlignmentMode,
}
impl Default for PoaConfig {
    fn default() -> Self {
        PoaConfig {
            band_width: 50,
            adaptive_band: true,
            adaptive_band_b: 10,
            adaptive_band_f: 0.01,

            // // abPOA-style scoring: gaps and mismatches are harsh relative to
            // // match (+2), so the aligner does not open cheap spurious gaps in
            // // homopolymer/periodic runs. The old +1/-1/-2/-1 made gaps too cheap
            // // (a gap-open cost only 2 matches), scattering homopolymer alignments
            // // at high error and over-calling repeats; abPOA-like scoring roughly
            // // halves that error with no regression on clean inputs (validated on
            // // the robustness matrix + 3-way comparison, 2026-07-28).
            match_score: 2,
            mismatch_score: -4,
            gap_open: -4,
            gap_extend: -3,
            
            min_boundary_coverage: 0,
            alignment_mode: AlignmentMode::SemiGlobal,
        }
    }
}

/* -----------------------------------------------------------------------------
Poa graph engine support
----------------------------------------------------------------------------- */
const NEG_INF: i32 = i32::MIN / 4; // /4 so additions can't overflow

// data type aliases
type NodeIndex = usize;
type SeqCount  = u32;
type Weight    = i32;
type BaseByte  = u8;

/// An out-edge and how many sequences traversed it (Match/Insert founding 
/// traffic).
#[derive(Clone, Copy, Debug)]
struct Edge {
    to: NodeIndex,
    weight: Weight,
}

/// A single-base node in the graph and metadata of out-edges, precedessor
/// nodes, and coverage.
struct Node {
    base:    u8,
    out:     Vec<Edge>,
    inc:     Vec<NodeIndex>, // predecessor node indices
    aligned: Vec<NodeIndex>, // sibling nodes in the same MSA column (alternative bases)
    cov: SeqCount, // sequences that Matched/founded this node (per-base depth)
    del: SeqCount, // sequences that Deleted (skipped) this node
}

/// Alignment op from traceback (forward sequence order).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Op {
    Match(NodeIndex), // sequence base aligned to this node (match or mismatch)
    Ins(BaseByte),    // sequence base with no graph node
    Del(NodeIndex),   // graph node skipped by the sequence
}

/// Read a banded matrix cell (M or I). Returns NEG_INF for j=0 or out-of-band.
#[inline]
fn bget(flat: &[i32], off: &[usize], lo: &[usize], hi: &[usize], r: usize, j: usize) -> i32 {
    if j == 0 || j < lo[r] || j > hi[r] {
        NEG_INF
    } else {
        flat[off[r] + (j - lo[r])]
    }
}

/// Read a banded D cell. j=0 lives in the separate `d0` column.
#[inline]
fn bgetd(
    flat: &[i32],
    d0: &[i32],
    off: &[usize],
    lo: &[usize],
    hi: &[usize],
    r: usize,
    j: usize,
) -> i32 {
    if j == 0 {
        d0[r]
    } else if j < lo[r] || j > hi[r] {
        NEG_INF
    } else {
        flat[off[r] + (j - lo[r])]
    }
}

/* -------------------------------------------------------------------------
Poa graph engine
------------------------------------------------------------------------- */
/// A Poa consensus graph engine:
///  - instantiated with PoaConfig parameters and the anticipated required 
///    capacity using `let poa = Poa::with_capacity()`,
///  - seeded with a scaffold sequence using `poa.seed_new_graph()`,
///  - populated with additional sequences using `poa.add_seq()`,
///  - executed to consensus using `poa.get_heaviest_path()`, and 
///  - used iteratively beginning again with `poa.seed_new_graph()`.
pub struct Poa {
    cfg: PoaConfig,
    base_capacity: usize, // used for instantiating Vecs, provided via with_capacity()
    pub(crate) n_seqs: usize,
    nodes:    Vec<Node>,
    seq_lens: Vec<usize>,
    /// Per-sequence provenance: the ordered node indices each sequence occupies
    /// (Match/Insert nodes, in sequence order; Deletes skip). Seq 0 is the seed.
    /// Node indices are stable (nodes are only appended), so a path recorded
    /// when a sequence was added stays valid for the life of the graph. This 
    /// drives per-partition rebuild membership and the coverage invariant.
    seq_paths: Vec<Vec<NodeIndex>>,
    /// Per-sequence *matched* adjacencies: edge (a, b) is recorded only when the
    /// sequence reached `b` directly after `a` with NO intervening Delete — i.e. 
    /// a genuine match/insert adjacency, excluding delete-bypass "resume" edges.
    /// This is the poa2 equivalent of legacy `edge_reads` (matched-only axis):
    /// bubble detection and arm-membership phasing use these, NOT the unified
    /// `Edge.weight` (which folds in delete-bypass traffic and would make a
    /// length variant's short-allele bypass look like a competing arm).
    seq_matched_edges: Vec<Vec<(NodeIndex, NodeIndex)>>,
    mrow:  Vec<i32>, /// single pre-allocated Vecs that live as long as the graph engine
    irow:  Vec<i32>,
    drow:  Vec<i32>,
    mbrow: Vec<(u8, u32)>,
    dbrow: Vec<(u8, u32)>,
}
impl Poa {
    /* -------------------------------------------------------------------------
    public caller functions needed to build a consensus, in usage order 
    ------------------------------------------------------------------------- */
    /// Create a new reusable Poa graph engine with the indicated sequence and 
    /// sequence length capacity. The returned engine can be used for all consensus 
    /// building that uses the same configuration parameters. Capacity will
    /// grow and persist stably as needed if the number of sequences or their 
    /// lengths exceed the values provided here.
    pub fn with_capacity(
        cfg: PoaConfig,
        n_seqs:  usize,
        n_bases: usize,
    ) -> Self {
        Poa {
            cfg,
            base_capacity: n_bases,
            n_seqs: 0,
            nodes:     Vec::with_capacity(n_bases),
            seq_lens:  Vec::with_capacity(n_seqs),
            seq_paths: Vec::with_capacity(n_seqs),
            seq_matched_edges: Vec::with_capacity(n_seqs),
            mrow:  Vec::with_capacity(n_bases),
            irow:  Vec::with_capacity(n_bases),
            drow:  Vec::with_capacity(n_bases),
            mbrow: Vec::with_capacity(n_bases),
            dbrow: Vec::with_capacity(n_bases),
        }
    }

    /// (Re)seed a Poa graph with a new scaffold sequence. This action clears 
    /// any previously built graph to begin the next graph build iteration.
    /// 
    /// Note that if the seed is a reference scaffold rather than a read, the 
    /// scaffold will get an equal vote during consensus building, and 
    /// unresolvable ties will prefer to report the scaffold sequence. This is 
    /// desirable when performing three-strand error correction. Use one chosen 
    /// read as seed if reference scaffolding is not needed.
    /// 
    /// If not None, `band_width` will update that POA setting for the set of 
    /// sequences anticipated as consensus inputs for this resolver iteration.
    pub fn seed_new_graph(
        &mut self, 
        seed: &[BaseByte], 
        band_width: Option<usize>
    ){
        self.n_seqs = 1;
        self.nodes.clear();
        self.seq_lens.clear();
        self.seq_paths.clear();
        self.seq_matched_edges.clear();
        for (i, &b) in seed.iter().enumerate() {
            let inc = if i == 0 { vec![] } else { vec![i - 1] };
            let out = if i + 1 < seed.len() {
                vec![Edge {
                    to: i + 1,
                    weight: 1,
                }]
            } else {
                vec![]
            };
            self.nodes.push(Node {
                base: b,
                out,
                inc,
                aligned: vec![],
                cov: 1,
                del: 0,
            });
        }
        let seed_path: Vec<usize> = (0..seed.len()).collect();
        let seed_medges: Vec<(usize, usize)> = (0..seed.len().saturating_sub(1))
            .map(|i| (i, i + 1))
            .collect();
        self.seq_lens.push(seed.len());
        self.seq_paths.push(seed_path);
        self.seq_matched_edges.push(seed_medges);
        if let Some(band_width) = band_width{
            self.cfg.band_width = band_width;
        }
    }

    /// Reverse complement a sequence prior to adding it to the graph. Callers 
    /// are responsible for calling `Poa::reverse_complement()` as needed prior 
    /// to calling `poa.add_seq()` to ensure that all sequences are in the same 
    /// orientation as the graph seed/scaffold.
    pub fn reverse_complement(seq: &[BaseByte]) -> Vec<BaseByte> {
        seq.iter().rev().map(|&b| match b {
            b'A' | b'a' => b'T', // all bases are returned as uppercase
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'N' | b'n' => b'N', // S, W, N are self-complementary
            b'S' | b's' => b'S',
            b'W' | b'w' => b'W',
            b'R' | b'r' => b'Y', // complemented IUPAC codes (R↔Y, K↔M, B↔V, D↔H)
            b'Y' | b'y' => b'R',
            b'K' | b'k' => b'M',
            b'M' | b'm' => b'K',
            b'B' | b'b' => b'V',
            b'V' | b'v' => b'B',
            b'D' | b'd' => b'H',
            b'H' | b'h' => b'D',
            other => other,
        }).collect()
    }

    /// Align another sequence into the building graph. Its base content will be
    /// used to build the final consensus. All sequences must be in the same
    /// strand orientation as the graph seed/scaffold for results to be 
    /// meaningful (see `Poa::reverse_complement()`).
    pub fn add_seq(&mut self, seq: &[BaseByte]) {
        if seq.is_empty() { return; }
        let (topo, rank_of) = self.topo_order();
        let ops = self.align(seq, &topo, &rank_of);
        let (path, medges) = self.integrate(seq, &ops);
        self.seq_paths.push(path);
        self.seq_matched_edges.push(medges);
        self.seq_lens.push(seq.len());
        self.n_seqs += 1;
    }

    /// Return the base sequence of the heaviest-path graph consensus as UTF-8 
    /// base bytes after all sequences have been added.
    pub fn get_heaviest_path(&self) -> Vec<BaseByte> {
        self.heaviest_path_nodes()
            .iter()
            .map(|&nd| self.nodes[nd].base)
            .collect()
    }
    /* -------------------------------------------------------------------------
    internal graph functions 
    ------------------------------------------------------------------------- */
    /// Kahn's topological order over the current DAG. Returns (topo, rank_of).
    fn topo_order(&self) -> (Vec<usize>, Vec<usize>) {
        let n = self.nodes.len();
        let mut indeg = vec![0usize; n];
        for nd in &self.nodes {
            for e in &nd.out {
                indeg[e.to] += 1;
            }
        }
        let mut stack: Vec<usize> = (0..n).filter(|&i| indeg[i] == 0).collect();
        stack.sort_unstable(); // deterministic
        let mut topo = Vec::with_capacity(n);
        let mut rank_of = vec![0usize; n];
        while let Some(u) = stack.pop() {
            rank_of[u] = topo.len();
            topo.push(u);
            for e in &self.nodes[u].out {
                indeg[e.to] -= 1;
                if indeg[e.to] == 0 {
                    stack.push(e.to);
                }
            }
        }
        (topo, rank_of)
    }

    /// Full partial-order affine-gap GLOBAL alignment of `seq` to the graph.
    /// Standard Gotoh with three states, generalized to a DAG: node `t`'s M/D
    /// relax over all predecessors. Returns ops in forward order.
    fn align(&mut self, seq: &[BaseByte], topo: &[usize], rank_of: &[usize]) -> Vec<Op> {
        const VSTART: u32 = u32::MAX; // pred-rank sentinel: virtual start
        let n_nodes = self.nodes.len();
        let n_bases = seq.len();
        let (mm, mis, go, ge) = (
            self.cfg.match_score,
            self.cfg.mismatch_score,
            self.cfg.gap_open,
            self.cfg.gap_extend,
        );
        let semi = self.cfg.alignment_mode == AlignmentMode::SemiGlobal;

        // ---- Adaptive static-diagonal-union band (abPOA-style anti-fold) ----
        // Per row we score only query columns [lo, hi] = union of the
        // predecessor-following adaptive center with the static graph-geometry
        // diagonal, padded by `half`. The static diagonal is fixed by graph
        // geometry so it is ALWAYS in-band no matter how far the adaptive center
        // drifts -- this is what stops a repeat fold (abPOA GET_AD_DP_*).
        // Storage is BANDED: each row keeps only its [lo, hi] window in flat
        // arrays with per-row offsets, so memory is O(nodes × band) not
        // O(nodes × seq_len).
        let banded0 = self.cfg.band_width > 0 || self.cfg.adaptive_band;
        const MIN_BANDWIDTH: usize = 5;
        let half = if self.cfg.adaptive_band{
            let adaptive =
                self.cfg.adaptive_band_b + (self.cfg.adaptive_band_f * n_bases as f32).ceil() as usize;
            adaptive.max(self.cfg.band_width).max(MIN_BANDWIDTH)
        } else {
            self.cfg.band_width.max(MIN_BANDWIDTH) 
        };

        let mut remain = vec![0usize; n_nodes]; // heaviest-out-edge node-distance to sink
        for &t in topo.iter().rev() {
            let mut best_w = -1i32;
            let mut best_c: Option<usize> = None;
            for e in &self.nodes[t].out {
                if e.weight > best_w {
                    best_w = e.weight;
                    best_c = Some(e.to);
                }
            }
            remain[t] = best_c.map_or(0, |c| remain[c] + 1);
        }
        // Band only spanning seqs: `n_bases - remain` assumes the seq covers the
        // graph (false for a short partial seq). Partials are short -> unbanded.
        let graph_span = remain.iter().max().copied().unwrap_or(0) + 1;
        let banded = banded0 && (n_bases as f64) >= 0.8 * (graph_span as f64);

        // Banded per-row storage (indexed by rank). off[r]..off[r]+width is the
        // row's window; lo_[r]/hi_[r] its query-column bounds. j=0 column (D only,
        // leading deletions) is stored separately in d0/d0bk.
        let mut lo_ = vec![1usize; n_nodes];
        let mut hi_ = vec![0usize; n_nodes];
        let mut off = vec![0usize; n_nodes];
        let cap = if banded {
            n_nodes * (2 * half + 8).min(n_bases + 1)
        } else {
            n_nodes * (n_bases + 1)
        };
        let mut mflat: Vec<i32> = Vec::with_capacity(cap);
        let mut iflat: Vec<i32> = Vec::with_capacity(cap);
        let mut dflat: Vec<i32> = Vec::with_capacity(cap);
        let mut mbk: Vec<(u8, u32)> = Vec::with_capacity(cap);
        let mut dbk: Vec<(u8, u32)> = Vec::with_capacity(cap);
        let mut d0 = vec![NEG_INF; n_nodes];
        let mut d0bk = vec![(0u8, VSTART); n_nodes];
        let mut best_j = vec![1usize; n_nodes]; // per-rank winning column (adaptive center)

        for &t in topo.iter() {
            let r = rank_of[t];
            let base = self.nodes[t].base;
            let is_source = self.nodes[t].inc.is_empty();
            let (lo, hi) = if banded {
                let adaptive = if is_source {
                    1
                } else {
                    self.nodes[t]
                        .inc
                        .iter()
                        .map(|&p| best_j[rank_of[p]] + 1)
                        .max()
                        .unwrap_or(1)
                };
                let stat = n_bases.saturating_sub(remain[t]).clamp(1, n_bases.max(1));
                (
                    adaptive.min(stat).saturating_sub(half).max(1),
                    (adaptive.max(stat) + half).min(n_bases).max(1),
                )
            } else {
                (1, n_bases.max(1))
            };
            let hi = hi.max(lo);
            let width = hi + 1 - lo;
            lo_[r] = lo;
            hi_[r] = hi;
            off[r] = mflat.len();

            // ---- j = 0 column: D0 only (leading deletion chain) ----
            {
                let mut best = NEG_INF;
                let mut bk = (0u8, VSTART);
                if is_source {
                    // leading deletion from the virtual start
                    if go + ge > best {
                        best = go + ge;
                        bk = (0, VSTART);
                    }
                }
                for &p in &self.nodes[t].inc {
                    let pr = rank_of[p] as u32;
                    let de = d0[rank_of[p]];
                    if de != NEG_INF && de + ge > best {
                        best = de + ge;
                        bk = (2, pr);
                    }
                }
                d0[r] = best;
                d0bk[r] = bk;
            }
            self.mrow.clear();
            self.mrow.resize(width, NEG_INF);
            self.irow.clear();
            self.irow.resize(width, NEG_INF);
            self.drow.clear();
            self.drow.resize(width, NEG_INF);
            self.mbrow.clear();
            self.mbrow.resize(width, (0u8, VSTART));
            self.dbrow.clear();
            self.dbrow.resize(width, (0u8, VSTART));
            let mut row_best = NEG_INF;
            let mut row_best_j = best_j[r];

            for j in lo..=hi {
                let k = j - lo;
                // ---- M[t][j] ----
                let sc = if base == seq[j - 1] { mm } else { mis };
                let mut best = NEG_INF;
                let mut bk = (0u8, VSTART);
                if (is_source || semi) && j == 1 && sc > best {
                    // free start: seq[0] begins the alignment here
                    best = sc;
                    bk = (0, VSTART);
                }
                for &p in &self.nodes[t].inc {
                    let pr = rank_of[p];
                    let vm = bget(&mflat, &off, &lo_, &hi_, pr, j - 1);
                    let vi = bget(&iflat, &off, &lo_, &hi_, pr, j - 1);
                    let vd = bgetd(&dflat, &d0, &off, &lo_, &hi_, pr, j - 1);
                    for (st, v) in [(0u8, vm), (1u8, vi), (2u8, vd)] {
                        if v != NEG_INF && v + sc > best {
                            best = v + sc;
                            bk = (st, pr as u32);
                        }
                    }
                }
                self.mrow[k] = best;
                self.mbrow[k] = bk;
                // ---- I[t][j] ---- (seq base inserted; same row, j-1) ----
                {
                    let (mo, io) = if j > lo {
                        (self.mrow[k - 1], self.irow[k - 1])
                    } else {
                        (NEG_INF, NEG_INF)
                    };
                    let mut best = NEG_INF;
                    if mo != NEG_INF {
                        best = mo + go + ge;
                    }
                    if io != NEG_INF && io + ge > best {
                        best = io + ge;
                    }
                    self.irow[k] = best;
                }
                // ---- D[t][j] ---- (node t skipped; predecessors at same j) ----
                {
                    let mut best = NEG_INF;
                    let mut bk = (0u8, VSTART);
                    for &p in &self.nodes[t].inc {
                        let pr = rank_of[p];
                        let mo = bget(&mflat, &off, &lo_, &hi_, pr, j);
                        let io = bget(&iflat, &off, &lo_, &hi_, pr, j);
                        let de = bgetd(&dflat, &d0, &off, &lo_, &hi_, pr, j);
                        if mo != NEG_INF && mo + go + ge > best {
                            best = mo + go + ge;
                            bk = (0, pr as u32);
                        }
                        if io != NEG_INF && io + go + ge > best {
                            best = io + go + ge;
                            bk = (1, pr as u32);
                        }
                        if de != NEG_INF && de + ge > best {
                            best = de + ge;
                            bk = (2, pr as u32);
                        }
                    }
                    self.drow[k] = best;
                    self.dbrow[k] = bk;
                }
                let cell = self.mrow[k].max(self.irow[k]).max(self.drow[k]);
                if cell != NEG_INF && cell > row_best {
                    row_best = cell;
                    row_best_j = j;
                }
            }
            best_j[r] = row_best_j;
            mflat.extend_from_slice(&self.mrow);
            iflat.extend_from_slice(&self.irow);
            dflat.extend_from_slice(&self.drow);
            mbk.extend_from_slice(&self.mbrow);
            dbk.extend_from_slice(&self.dbrow);
        }

        // Terminal (seq fully consumed at j=l). Global: end at a sink; semi:
        // end at ANY node (free trailing graph gap).
        let mut best = NEG_INF;
        let mut cur: (u8, u32) = (0, VSTART);
        for &t in topo.iter() {
            if !semi && !self.nodes[t].out.is_empty() {
                continue;
            }
            let r = rank_of[t];
            let vm = bget(&mflat, &off, &lo_, &hi_, r, n_bases);
            let vi = bget(&iflat, &off, &lo_, &hi_, r, n_bases);
            let vd = bgetd(&dflat, &d0, &off, &lo_, &hi_, r, n_bases);
            for (st, v) in [(0u8, vm), (1u8, vi), (2u8, vd)] {
                if v != NEG_INF && v > best {
                    best = v;
                    cur = (st, r as u32);
                }
            }
        }

        // Traceback.
        let mut ops: Vec<Op> = Vec::with_capacity(self.base_capacity);
        let (mut state, mut rr) = cur;
        let mut j = n_bases;
        while rr != VSTART {
            let r = rr as usize;
            let t = topo[r];
            match state {
                0 => {
                    ops.push(Op::Match(t));
                    let (pst, pr) = mbk[off[r] + (j - lo_[r])];
                    j -= 1;
                    state = pst;
                    rr = pr;
                }
                1 => {
                    ops.push(Op::Ins(seq[j - 1]));
                    j -= 1;
                    if j == 0 {
                        rr = VSTART;
                    } else {
                        let o = bget(&mflat, &off, &lo_, &hi_, r, j);
                        let e = bget(&iflat, &off, &lo_, &hi_, r, j);
                        // The forward recurrence is I[j] = max(M[j-1]+go+ge,
                        // I[j-1]+ge), so arriving from M pays the gap-open. The
                        // predecessor is M only when M[j] + go >= I[j]; comparing
                        // raw `o >= e` (go omitted) picks M too eagerly and cuts a
                        // multi-base insert run one base short, mis-emitting an
                        // inserted base as a Match. Guard the NEG_INF M cell so an
                        // unreachable M can never win the comparison.
                        state = if o != NEG_INF && o + go >= e { 0 } else { 1 };
                    }
                }
                _ => {
                    ops.push(Op::Del(t));
                    let (pst, pr) = if j == 0 {
                        d0bk[r]
                    } else {
                        dbk[off[r] + (j - lo_[r])]
                    };
                    state = pst;
                    rr = pr;
                }
            }
            if rr == VSTART {
                while j > 0 {
                    ops.push(Op::Ins(seq[j - 1]));
                    j -= 1;
                }
                break;
            }
        }
        ops.reverse();
        ops
    }

    /// Add an edge to the graph.
    fn add_edge(&mut self, from: NodeIndex, to: NodeIndex) {
        for edge in &mut self.nodes[from].out {
            if edge.to == to {
                edge.weight += 1;
                return;
            }
        }
        self.nodes[from].out.push(Edge { to, weight: 1 });
        self.nodes[to].inc.push(from);
    }

    /// Add a node to the graph.
    fn new_node(&mut self, base: BaseByte) -> NodeIndex {
        let node_index = self.nodes.len();
        self.nodes.push(Node {
            base,
            out: vec![],
            inc: vec![],
            aligned: vec![],
            cov: 0,
            del: 0,
        });
        node_index
    }

    /// SPOA/abPOA-style node fusion: the seq base aligned (by the DP) to
    /// column-node `t`. If `t` (or one of its aligned siblings) already carries
    /// `base`, reuse that node; otherwise create a new substitution node and
    /// mutually cross-link it into `t`'s alignment column. This keeps all
    /// alternative bases at a position as siblings, so later sequences reuse 
    /// them (correct allele counting) instead of fragmenting into fresh nodes.
    fn aligned_or_new(&mut self, t: NodeIndex, base: BaseByte) -> NodeIndex {
        if self.nodes[t].base == base { return t; }
        for &a in &self.nodes[t].aligned {
            if self.nodes[a].base == base { return a; }
        }
        let base_node_index = self.new_node(base);
        let mut column = self.nodes[t].aligned.clone();
        column.push(t);
        for &column_node_index in &column {
            self.nodes[column_node_index].aligned.push(base_node_index);
            self.nodes[base_node_index].aligned.push(column_node_index);
        }
        base_node_index
    }

    /// Integrate a seq's alignment into the graph and return `(path, medges)`:
    /// the ordered node indices it occupies (Match/Insert nodes, Deletes
    /// skipped) and its *matched* adjacencies (edges reached with no intervening
    /// Delete — the delete-bypass "resume" edge is excluded). The caller stores
    /// these in `seq_paths` / `seq_matched_edges`.
    fn integrate(&mut self, seq: &[BaseByte], ops: &[Op]) -> (Vec<NodeIndex>, Vec<(NodeIndex, NodeIndex)>) {
        // `prev` = last graph node the seq is currently attached to.
        // Walk seq positions: Match and Ins each consume one seq base; Del none.
        let mut prev: Option<NodeIndex> = None;
        let mut j = 0_usize;
        let mut path = Vec::with_capacity(seq.len());
        let mut medges = Vec::with_capacity(seq.len());
        // Whether a Delete occurred since the last occupied node: if so the next
        // Match/Ins reconnects via a bypass edge, which is NOT a matched
        // adjacency (mirrors legacy pure-bypass resume).
        let mut bypassed = false;
        for op in ops {
            match *op {
                Op::Match(t) => {
                    let base = seq[j];
                    j += 1;
                    let node_index = self.aligned_or_new(t, base); // reuse on match, sub-node on mismatch
                    self.nodes[node_index].cov += 1;
                    if let Some(from) = prev {
                        self.add_edge(from, node_index);
                        if !bypassed {
                            medges.push((from, node_index));
                        }
                    }
                    prev = Some(node_index);
                    path.push(node_index);
                    bypassed = false;
                }
                Op::Ins(_) => {
                    let base = seq[j];
                    j += 1;
                    let node_index = self.new_node(base);
                    self.nodes[node_index].cov += 1;
                    if let Some(from) = prev {
                        self.add_edge(from, node_index);
                        if !bypassed {
                            medges.push((from, node_index));
                        }
                    }
                    prev = Some(node_index);
                    path.push(node_index);
                    bypassed = false;
                }
                Op::Del(t) => {
                    // node skipped by this seq: record the delete (used by the
                    // analysis layer's Match-vs-Delete column entropy) and mark
                    // that the next reconnection is a bypass, not a match.
                    self.nodes[t].del += 1;
                    bypassed = true;
                }
            }
        }
        (path, medges)
    }

    /// abPOA-style heaviest bundling: reverse pass computes, per node, the
    /// heaviest out-edge (tie-broken by giving preference to first recorded 
    /// best out-edge from a node, i.e., to the scaffold sequence); walk the  
    /// chosen chain from a source. Return the consensus node path as indices.
    fn heaviest_path_nodes(&self) -> Vec<NodeIndex> {
        let (topo, _) = self.topo_order();
        let n_nodes = self.nodes.len();
        if n_nodes == 0 { return vec![]; }
        let mut scores = vec![0_i64; n_nodes]; 
        let mut next_node_indices = vec![NodeIndex::MAX; n_nodes];
        for &t in topo.iter().rev() {
            let mut best_weight: Weight = -1;
            let mut best_score = i64::MIN;
            let mut best_column = NodeIndex::MAX;
            // the first out-edge for a scaffold node will be to the next scaffold node
            for edge in &self.nodes[t].out { 
                let score = scores[edge.to];
                // eprintln!("node {}, {} to {}, score {},  {:?}", 
                //     t,
                //     std::str::from_utf8(&[self.nodes[t].base]).unwrap(), 
                //     std::str::from_utf8(&[self.nodes[edge.to].base]).unwrap(), 
                //     score, 
                //     edge
                // );
                // for a scaffold node, this loop triggers once for the scaffold-to-scaffold edge
                // and won't trigger again unless scaffold-to-alt has a heavier weight
                // TODO: consider bringing back score tie-breaker only at non-scaffold t nodes
                if edge.weight > best_weight {
                    best_weight = edge.weight;
                    best_score  = score;
                    best_column = edge.to;
                }
            }
            if best_column != NodeIndex::MAX {
                scores[t] = best_weight as i64 + best_score;
                next_node_indices[t] = best_column;
            }
            // eprintln!("{} to {}, score {}", t, next_node_indices[t], scores[t]);
        }
        let mut start = topo[0];
        let mut best_score = i64::MIN;
        for &t in &topo {
            if self.nodes[t].inc.is_empty() && scores[t] > best_score {
                best_score = scores[t];
                start = t;
            }
        }
        let mut path = Vec::with_capacity(self.base_capacity);
        let mut current = start;
        // `visited` is cycle-safety insurance: the graph is a DAG by construction
        // (topo_order debug-asserts it), so `next_node_indices` cannot cycle — 
        // but node fusion has no explicit back-edge guard, so if a cycle ever 
        // slipped through in a release build this bounds the walk instead of 
        // looping.
        let mut visited = vec![false; n_nodes];
        while current != NodeIndex::MAX && !visited[current] {
            visited[current] = true;
            path.push(current);
            current = next_node_indices[current];
        }
        // Boundary trim: drop leading/trailing consensus nodes below
        // min_boundary_coverage. This removes low-coverage flanks and 
        // trailing/leading repeat-unit extensions supported by only a minority 
        // of seqs (e.g. a few longer seqs over-extending a homogeneous 
        // repeat). abPOA/SPOA do the equivalent. Interior low-coverage (a 
        // genuine spanning-seqs gap) is preserved.
        let floor = self.cfg.min_boundary_coverage;
        let s = path.iter().position( |&nd| self.nodes[nd].cov >= floor);
        let e = path.iter().rposition(|&nd| self.nodes[nd].cov >= floor);
        match (s, e) {
            (Some(s), Some(e)) if s <= e => path[s..=e].to_vec(),
            _ => path,
        }
    }
}
