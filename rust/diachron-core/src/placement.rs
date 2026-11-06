//! Two-state parsimony scores for every insertion edge of a detached subtree.

use crate::tree::Tree;

type StateCosts = [u32; 2];

const IMPOSSIBLE: u32 = 1_000_000;

fn project(cost: StateCosts) -> StateCosts {
    [cost[0].min(cost[1] + 1), cost[1].min(cost[0] + 1)]
}

fn walk(root: usize, children: &[Vec<usize>]) -> Vec<usize> {
    let mut result = vec![root];
    let mut next = 0;
    while next < result.len() {
        result.extend(&children[result[next]]);
        next += 1;
    }
    result
}

fn subtree_costs(
    tree: &Tree,
    children: &[Vec<usize>],
    order: &[usize],
    column: &[Option<usize>],
    costs: &mut [StateCosts],
) {
    for &node in order.iter().rev() {
        let mut state = [0, 0];
        if let Some(taxon) = tree.nodes[node].taxon {
            if let Some(observed) = column[taxon] {
                state[1 - observed] = IMPOSSIBLE;
            }
        } else {
            for &child in &children[node] {
                let contribution = project(costs[child]);
                state[0] += contribution[0];
                state[1] += contribution[1];
            }
        }
        costs[node] = state;
    }
}

/// Return each reachable backbone node's total placement cost, with the
/// backbone root representing the stem. Detached nodes have no score.
pub fn placement_costs(
    tree: &Tree,
    detached: usize,
    columns: &[Vec<Option<usize>>],
) -> Result<Vec<Option<u64>>, &'static str> {
    let n = tree.nodes.len();
    let junction = tree
        .nodes
        .get(detached)
        .and_then(|node| node.parent)
        .ok_or("cannot detach the root")?;
    let siblings = &tree.nodes[junction].children;
    if siblings.len() != 2 {
        return Err("placement requires a binary attachment junction");
    }
    let sibling = *siblings
        .iter()
        .find(|&&node| node != detached)
        .ok_or("attachment junction has no sibling")?;
    let mut children: Vec<_> = tree
        .nodes
        .iter()
        .map(|node| node.children.clone())
        .collect();
    let root = if let Some(parent) = tree.nodes[junction].parent {
        let entry = children[parent]
            .iter_mut()
            .find(|node| **node == junction)
            .ok_or("attachment junction is missing from its parent")?;
        *entry = sibling;
        tree.root
    } else {
        sibling
    };
    let detached_order = walk(detached, &children);
    let backbone = walk(root, &children);
    let mut costs = vec![[0, 0]; n];
    let mut outside = vec![[0, 0]; n];
    let mut totals = vec![None; n];
    for &node in &backbone {
        totals[node] = Some(0);
    }
    for column in columns {
        subtree_costs(tree, &children, &detached_order, column, &mut costs);
        subtree_costs(tree, &children, &backbone, column, &mut costs);
        let detached = project(costs[detached]);
        let stem = project(costs[root]);
        *totals[root].as_mut().unwrap() +=
            u64::from((stem[0] + detached[0]).min(stem[1] + detached[1]));
        outside[root] = [0, 0];
        for &parent in &backbone {
            let mut total = outside[parent];
            for &child in &children[parent] {
                let contribution = project(costs[child]);
                total[0] += contribution[0];
                total[1] += contribution[1];
            }
            for &child in &children[parent] {
                let contribution = project(costs[child]);
                let side = project([total[0] - contribution[0], total[1] - contribution[1]]);
                let absent = side[0] + contribution[0] + detached[0];
                let present = side[1] + contribution[1] + detached[1];
                *totals[child].as_mut().unwrap() += u64::from(absent.min(present));
                outside[child] = side;
            }
        }
    }
    Ok(totals)
}
