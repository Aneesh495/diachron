//! Two-state parsimony scores for every insertion edge of a detached subtree.

use crate::subst::{BinaryGainLoss, SubstModel};
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

struct Backbone {
    root: usize,
    children: Vec<Vec<usize>>,
    order: Vec<usize>,
    detached_order: Vec<usize>,
}

fn detach(tree: &Tree, detached: usize) -> Result<Backbone, &'static str> {
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
    let order = walk(root, &children);
    Ok(Backbone {
        root,
        children,
        order,
        detached_order,
    })
}

/// Return each reachable backbone node's total placement cost, with the
/// backbone root representing the stem. Detached nodes have no score.
pub fn placement_costs(
    tree: &Tree,
    detached: usize,
    columns: &[Vec<Option<usize>>],
) -> Result<Vec<Option<u64>>, &'static str> {
    let n = tree.nodes.len();
    let Backbone {
        root,
        children,
        order: backbone,
        detached_order,
    } = detach(tree, detached)?;
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

type Pair = [f64; 2];
type Transition = [Pair; 2];

fn transition(model: &BinaryGainLoss, duration: f64) -> Transition {
    let p = model.transition(duration);
    [[p.at(0, 0), p.at(0, 1)], [p.at(1, 0), p.at(1, 1)]]
}

fn dot(left: Pair, right: Pair) -> f64 {
    left[0] * right[0] + left[1] * right[1]
}

fn multiply(left: Pair, right: Pair) -> Pair {
    [left[0] * right[0], left[1] * right[1]]
}

fn normalize(value: Pair) -> Pair {
    let scale = value[0].max(value[1]);
    if scale > 0.0 {
        [value[0] / scale, value[1] / scale]
    } else {
        value
    }
}

fn up(p: &Transition, value: Pair) -> Pair {
    [dot(p[0], value), dot(p[1], value)]
}

fn down(p: &Transition, value: Pair) -> Pair {
    [
        p[0][0] * value[0] + p[1][0] * value[1],
        p[0][1] * value[0] + p[1][1] * value[1],
    ]
}

fn likelihood_messages(
    tree: &Tree,
    children: &[Vec<usize>],
    order: &[usize],
    column: &[Option<usize>],
    transitions: &[Transition],
    inside: &mut [Pair],
) {
    for &node in order.iter().rev() {
        inside[node] = if let Some(taxon) = tree.nodes[node].taxon {
            match column[taxon] {
                Some(0) => [1., 0.],
                Some(1) => [0., 1.],
                _ => [1., 1.],
            }
        } else {
            children[node].iter().fold([1., 1.], |value, &child| {
                normalize(multiply(value, up(&transitions[child], inside[child])))
            })
        };
    }
}

struct Split {
    above: Transition,
    below: Transition,
    detached: Transition,
}

/// Compare vertical likelihood at the midpoint of three age bands per edge.
/// The reference mass of the backbone and detached message scale cancel
/// across placements. The score is used only as a proposal guide.
pub fn placement_likelihoods(
    tree: &Tree,
    detached: usize,
    model: &BinaryGainLoss,
    columns: &[Vec<Option<usize>>],
) -> Result<Vec<Option<[f64; 3]>>, &'static str> {
    let Backbone {
        root,
        children,
        order,
        detached_order,
    } = detach(tree, detached)?;
    let n = tree.nodes.len();
    let height = |node: usize| tree.nodes[node].height;
    let mut transitions = vec![[[0., 0.]; 2]; n];
    for &parent in order.iter().chain(&detached_order) {
        for &child in &children[parent] {
            let duration = height(parent) - height(child);
            if duration < 0.0 {
                return Err("invalid ages on placement backbone");
            }
            transitions[child] = transition(model, duration);
        }
    }
    let mut splits: Vec<Option<[Split; 3]>> = (0..n).map(|_| None).collect();
    let mut scores = vec![None; n];
    for &parent in &order {
        for &child in &children[parent] {
            let lower = height(child).max(height(detached));
            let upper = height(parent);
            if lower < upper {
                splits[child] = Some(std::array::from_fn(|band| {
                    let age = lower + (upper - lower) * (2 * band + 1) as f64 / 6.;
                    Split {
                        above: transition(model, upper - age),
                        below: transition(model, age - height(child)),
                        detached: transition(model, age - height(detached)),
                    }
                }));
                scores[child] = Some([0.; 3]);
            }
        }
    }
    let mut inside = vec![[0., 0.]; n];
    let mut outside = vec![[0., 0.]; n];
    let mut suffix = Vec::new();
    let stationary = model.stationary();
    for column in columns {
        likelihood_messages(
            tree,
            &children,
            &detached_order,
            column,
            &transitions,
            &mut inside,
        );
        likelihood_messages(tree, &children, &order, column, &transitions, &mut inside);
        outside[root] = [stationary[0], stationary[1]];
        for &parent in &order {
            let siblings = &children[parent];
            suffix.resize(siblings.len() + 1, [1., 1.]);
            suffix[siblings.len()] = [1., 1.];
            for (i, &child) in siblings.iter().enumerate().rev() {
                suffix[i] = normalize(multiply(
                    suffix[i + 1],
                    up(&transitions[child], inside[child]),
                ));
            }
            let mut prefix = outside[parent];
            for (i, &child) in siblings.iter().enumerate() {
                let context = normalize(multiply(prefix, suffix[i + 1]));
                let child_message = inside[child];
                let parent_to_child = down(&transitions[child], context);
                if let (Some(splits), Some(scores)) = (&splits[child], &mut scores[child]) {
                    let reference = dot(parent_to_child, child_message).ln();
                    for (split, score) in splits.iter().zip(scores) {
                        let from_parent = down(&split.above, context);
                        let from_child = up(&split.below, child_message);
                        let from_detached = up(&split.detached, inside[detached]);
                        *score +=
                            dot(multiply(from_parent, from_child), from_detached).ln() - reference;
                    }
                }
                outside[child] = normalize(parent_to_child);
                prefix = normalize(multiply(prefix, up(&transitions[child], child_message)));
            }
        }
    }
    if scores
        .iter()
        .flatten()
        .flatten()
        .any(|value| !value.is_finite())
    {
        return Err("nonfinite likelihood placement score");
    }
    Ok(scores)
}
