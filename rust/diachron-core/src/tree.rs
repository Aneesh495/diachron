//! Rooted phylogenetic tree representation.
//!
//! Nodes are stored in a flat arena. Each node knows its parent branch length
//! (the length of the edge above it, in time or substitution units depending on
//! the caller's convention) and its children. Tips (leaves) carry a taxon index
//! used to look up observed character states. Internal nodes may be
//! multifurcating.
//!
//! Ownership: the `Tree` owns all nodes. Indices are stable for the lifetime of
//! the tree. The root has no parent edge (its `branch_len` is ignored).
//!
//! Node times: `height[i]` is the age of node i, in years before present, with
//! tips usually at 0 (contemporary) or at a positive age (dated/extinct tips).
//! Branch length above node i equals `height[parent] - height[i]` when times are
//! used; callers that work directly in substitution units set `branch_len`.

/// A single node in the arena.
#[derive(Clone, Debug)]
pub struct Node {
    pub parent: Option<usize>,
    pub children: Vec<usize>,
    /// Edge length above this node. For tips and internal nodes in a
    /// substitution-unit tree this is set directly; in a timed tree it is
    /// derived from node heights.
    pub branch_len: f64,
    /// Taxon index if this is a tip, else None.
    pub taxon: Option<usize>,
    /// Node age in years before present (0 = contemporary tip).
    pub height: f64,
}

#[derive(Clone, Debug)]
pub struct Tree {
    pub nodes: Vec<Node>,
    pub root: usize,
}

impl Tree {
    pub fn new() -> Self {
        Tree { nodes: Vec::new(), root: 0 }
    }

    pub fn add_node(&mut self, parent: Option<usize>, branch_len: f64, taxon: Option<usize>) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(Node {
            parent,
            children: Vec::new(),
            branch_len,
            taxon,
            height: 0.0,
        });
        if let Some(p) = parent {
            self.nodes[p].children.push(idx);
        }
        idx
    }

    pub fn is_tip(&self, idx: usize) -> bool {
        self.nodes[idx].children.is_empty()
    }

    pub fn n_tips(&self) -> usize {
        self.nodes.iter().filter(|n| n.children.is_empty()).count()
    }

    /// Post-order traversal (children before parents), returned as node indices.
    /// This is the order the pruning recursion consumes.
    pub fn postorder(&self) -> Vec<usize> {
        let mut order = Vec::with_capacity(self.nodes.len());
        let mut stack = vec![(self.root, false)];
        while let Some((node, processed)) = stack.pop() {
            if processed {
                order.push(node);
            } else {
                stack.push((node, true));
                for &c in &self.nodes[node].children {
                    stack.push((c, false));
                }
            }
        }
        order
    }

    /// Recompute `branch_len` for every node from node heights. A branch length
    /// is parent height minus child height; it must be nonnegative for a valid
    /// time tree, which the caller is expected to enforce via constraints.
    pub fn set_branch_lengths_from_heights(&mut self) {
        let n = self.nodes.len();
        for i in 0..n {
            if let Some(p) = self.nodes[i].parent {
                self.nodes[i].branch_len = self.nodes[p].height - self.nodes[i].height;
            }
        }
    }

    /// True if every child is strictly younger than its parent. Zero-length
    /// branches are allowed (sampled ancestors / hard polytomies) but negatives
    /// are not.
    pub fn times_valid(&self) -> bool {
        for (i, node) in self.nodes.iter().enumerate() {
            if let Some(p) = node.parent {
                if self.nodes[p].height < node.height - 1e-9 {
                    return false;
                }
                let _ = i;
            }
        }
        true
    }
}

impl Default for Tree {
    fn default() -> Self {
        Self::new()
    }
}
