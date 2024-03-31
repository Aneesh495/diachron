//! Felsenstein pruning likelihood with log-space scaling.
//!
//! Given a rooted tree, a substitution model, and observed states at the tips,
//! `log_likelihood` computes the probability of the observed character under the
//! model, marginalizing over all internal node states. The recursion stores, for
//! each node and each state s, the conditional likelihood
//!
//!     L_node(s) = P(data below node | node is in state s),
//!
//! computed bottom up. For a tip the vector is an indicator of the observed
//! state; a missing observation is the all-ones vector (every state is equally
//! consistent with "we did not observe this cell"), which is what makes missing
//! and absent produce different likelihoods. At an internal node with children
//! c, the per-child contribution is sum_j P_ij(t_c) L_c(j), and the node vector
//! is the product of child contributions.
//!
//! Numerical stability: products of many small probabilities underflow. At each
//! internal node we factor out the maximum entry of the conditional vector into
//! an accumulated log-scale, so the stored vector stays O(1). The final log
//! likelihood adds the accumulated scales back.
//!
//! Multifurcations are handled directly: an internal node simply has more than
//! two children and the product runs over all of them.

use crate::linalg::Mat;
use crate::subst::SubstModel;
use crate::tree::Tree;

mod gradient;
pub use gradient::binary_mixture_branch_gradient;

/// Observed tip states. `states[taxon]` is `Some(s)` for an observed state `s`
/// in `0..n_states`, or `None` for a missing observation.
pub struct TipStates<'a> {
    pub states: &'a [Option<usize>],
}

/// Compute the log-likelihood of a single character on a tree under a model.
///
/// `root_prior` is the probability distribution over the root state (usually the
/// model's stationary distribution). The returned value is a natural log.
pub fn log_likelihood(
    tree: &Tree,
    model: &dyn SubstModel,
    tips: &TipStates,
    root_prior: &[f64],
) -> f64 {
    PruningPlan::new(tree, model, 1.0).score(tips.states, root_prior)
}

/// Tree traversal and transitions shared across sites of one rate category.
pub(crate) struct PruningPlan<'a> {
    pub(crate) tree: &'a Tree,
    pub(crate) order: Vec<usize>,
    pub(crate) transitions: Vec<Mat>,
    n_states: usize,
}

impl<'a> PruningPlan<'a> {
    pub(crate) fn new(tree: &'a Tree, model: &dyn SubstModel, rate: f64) -> Self {
        Self {
            tree,
            order: tree.postorder(),
            transitions: tree
                .nodes
                .iter()
                .map(|node| model.transition(node.branch_len * rate))
                .collect(),
            n_states: model.n_states(),
        }
    }

    pub(crate) fn score(&self, states: &[Option<usize>], root_prior: &[f64]) -> f64 {
        let tree = self.tree;
        let n_states = self.n_states;
        assert_eq!(root_prior.len(), n_states);
        let n_nodes = tree.nodes.len();
        // Conditional likelihood vectors, one length-n_states block per node.
        let mut cond = vec![0.0f64; n_nodes * n_states];
        // Accumulated log scale per node.
        let mut logscale = vec![0.0f64; n_nodes];

        for &node in &self.order {
            let base = node * n_states;
            if tree.is_tip(node) {
                let taxon = tree.nodes[node].taxon.expect("tip without taxon");
                match states[taxon] {
                    Some(s) => {
                        for k in 0..n_states {
                            cond[base + k] = if k == s { 1.0 } else { 0.0 };
                        }
                    }
                    None => {
                        // Missing: consistent with every state.
                        for k in 0..n_states {
                            cond[base + k] = 1.0;
                        }
                    }
                }
                logscale[node] = 0.0;
            } else {
                // Start from all-ones and multiply in each child's contribution.
                for k in 0..n_states {
                    cond[base + k] = 1.0;
                }
                let mut child_scale = 0.0;
                for &child in &tree.nodes[node].children {
                    let p = &self.transitions[child];
                    let cbase = child * n_states;
                    // partial[i] = sum_j P_ij(t) * L_child(j)
                    for i in 0..n_states {
                        let mut acc = 0.0;
                        for j in 0..n_states {
                            acc += p.at(i, j) * cond[cbase + j];
                        }
                        cond[base + i] *= acc;
                    }
                    child_scale += logscale[child];
                    // Scale after each child so a wide polytomy stays finite too.
                    let block = &mut cond[base..base + n_states];
                    let maxv = block.iter().copied().fold(0.0, f64::max);
                    if maxv <= 0.0 {
                        return f64::NEG_INFINITY;
                    }
                    for value in block {
                        *value /= maxv;
                    }
                    child_scale += maxv.ln();
                }
                logscale[node] = child_scale;
            }
        }

        // Combine at the root with the prior.
        let rbase = tree.root * n_states;
        let mut site = 0.0;
        for k in 0..n_states {
            site += root_prior[k] * cond[rbase + k];
        }
        if site <= 0.0 {
            return f64::NEG_INFINITY;
        }
        site.ln() + logscale[tree.root]
    }

    /// Binary workspaces are reused across columns and contact components.
    fn score_binary(
        &self,
        states: &[Option<usize>],
        prior: &[f64],
        partials: &mut [[f64; 2]],
    ) -> f64 {
        let mut log_scale = 0.0;
        for &node in &self.order {
            if let Some(taxon) = self.tree.nodes[node].taxon {
                partials[node] = match states[taxon] {
                    Some(0) => [1.0, 0.0],
                    Some(1) => [0.0, 1.0],
                    None => [1.0, 1.0],
                    _ => return f64::NEG_INFINITY,
                };
                continue;
            }

            let mut pair = [1.0, 1.0];
            for &child in &self.tree.nodes[node].children {
                let [absent, present] = partials[child];
                let transition = &self.transitions[child];
                pair[0] *= transition.at(0, 0) * absent + transition.at(0, 1) * present;
                pair[1] *= transition.at(1, 0) * absent + transition.at(1, 1) * present;
                let scale = pair[0].max(pair[1]);
                if scale <= 0.0 {
                    return f64::NEG_INFINITY;
                }
                pair = [pair[0] / scale, pair[1] / scale];
                log_scale += scale.ln();
            }
            partials[node] = pair;
        }

        // Every child subtree is used once, so all rescaling factors can be
        // accumulated globally rather than carried in a per-node buffer.
        let [absent, present] = partials[self.tree.root];
        (prior[0] * absent + prior[1] * present).ln() + log_scale
    }
}

/// Log-likelihood of a block of characters that share the same tree and model,
/// summed over sites. `columns[c][taxon]` is the state of taxon at character c.
pub fn log_likelihood_block(
    tree: &Tree,
    model: &dyn SubstModel,
    columns: &[Vec<Option<usize>>],
    root_prior: &[f64],
) -> f64 {
    let plan = PruningPlan::new(tree, model, 1.0);
    if model.n_states() == 2 {
        let mut partials = vec![[0.0; 2]; tree.nodes.len()];
        return columns
            .iter()
            .map(|col| plan.score_binary(col, root_prior, &mut partials))
            .sum();
    }
    columns.iter().map(|col| plan.score(col, root_prior)).sum()
}

const COUNT_BUCKETS: usize = 6;
const ZERO_FLAG: usize = 3;
type Counts = [f64; COUNT_BUCKETS];
type StateCounts = [Counts; 2];

fn combine_bucket(a: usize, b: usize) -> usize {
    let has_zero = (a >= ZERO_FLAG || b >= ZERO_FLAG) as usize;
    let positives = (a % ZERO_FLAG + b % ZERO_FLAG).min(2);
    has_zero * ZERO_FLAG + positives
}

fn combine_counts(left: &Counts, right: &Counts) -> Counts {
    let mut combined = [0.0; COUNT_BUCKETS];
    for (a, &pa) in left.iter().enumerate() {
        if pa == 0.0 {
            continue;
        }
        for (b, &pb) in right.iter().enumerate() {
            combined[combine_bucket(a, b)] += pa * pb;
        }
    }
    combined
}

fn child_counts(plan: &PruningPlan<'_>, child: usize, state: usize, message: &StateCounts) -> Counts {
    let mut result = [0.0; COUNT_BUCKETS];
    let transition = &plan.transitions[child];
    for (bucket, probability) in result.iter_mut().enumerate() {
        *probability = transition.at(state, 0) * message[0][bucket]
            + transition.at(state, 1) * message[1][bucket];
    }
    result
}

fn count_node<'a>(
    plan: &PruningPlan<'_>, node: usize, multiplicity: &[usize],
    message_at: impl Fn(usize) -> &'a StateCounts,
) -> StateCounts {
    let mut mass = [[0.0; COUNT_BUCKETS]; 2];
    if let Some(taxon) = plan.tree.nodes[node].taxon {
        let copies = multiplicity[taxon];
        mass[0][if copies > 0 { ZERO_FLAG } else { 0 }] = 1.0;
        mass[1][copies.min(2)] = 1.0;
    } else {
        mass[0][0] = 1.0;
        mass[1][0] = 1.0;
        for &child in &plan.tree.nodes[node].children {
            for (state, left) in mass.iter_mut().enumerate() {
                *left = combine_counts(left, &child_counts(plan, child, state, message_at(child)));
            }
        }
    }
    mass
}

fn fill_count_messages(plan: &PruningPlan<'_>, multiplicity: &[usize], messages: &mut [StateCounts]) {
    for &node in &plan.order {
        messages[node] = count_node(plan, node, multiplicity, |child| &messages[child]);
    }
}

fn count_messages(plan: &PruningPlan<'_>, multiplicity: &[usize]) -> Vec<StateCounts> {
    let mut messages = vec![[[0.0; COUNT_BUCKETS]; 2]; plan.tree.nodes.len()];
    fill_count_messages(plan, multiplicity, &mut messages);
    messages
}

fn root_retention(root: &StateCounts, min_present: usize, prior: &[f64]) -> f64 {
    root.iter().zip(prior)
        .map(|(mass, pi)| pi * mass[ZERO_FLAG + min_present..].iter().sum::<f64>())
        .sum()
}

/// Probability of the observed-mask retention event, computed directly in the
/// selected buckets to avoid subtraction on nearly constant trees.
fn retained_probability(
    plan: &PruningPlan<'_>, multiplicity: &[usize], min_present: usize, prior: &[f64],
) -> f64 {
    let messages = count_messages(plan, multiplicity);
    root_retention(&messages[plan.tree.root], min_present, prior)
}

struct ContactSelectionPath {
    affected: Vec<bool>,
    order: Vec<usize>,
}

struct ContactSelectionWorkspace {
    paths: Vec<ContactSelectionPath>,
    baseline: Vec<StateCounts>,
    changed: Vec<StateCounts>,
    retained: Vec<f64>,
}

impl ContactSelectionWorkspace {
    fn new(tree: &Tree, order: &[usize], contacts: &[(usize, usize, f64)]) -> Self {
        let mut tips = vec![0; tree.n_tips()];
        for (index, node) in tree.nodes.iter().enumerate() {
            if let Some(taxon) = node.taxon {
                tips[taxon] = index;
            }
        }
        let paths = contacts.iter().map(|&(donor, recipient, _)| {
            let mut affected = vec![false; tree.nodes.len()];
            for tip in [tips[donor], tips[recipient]] {
                let mut current = Some(tip);
                while let Some(node) = current {
                    if affected[node] {
                        break;
                    }
                    affected[node] = true;
                    current = tree.nodes[node].parent;
                }
            }
            let order = order.iter().copied().filter(|&node| affected[node]).collect();
            ContactSelectionPath { affected, order }
        }).collect();
        let buffer = vec![[[0.0; COUNT_BUCKETS]; 2]; tree.nodes.len()];
        Self { paths, baseline: buffer.clone(), changed: buffer, retained: vec![0.0; contacts.len() + 1] }
    }

    fn mixture_mass(
        &mut self, plans: &[PruningPlan<'_>], contacts: &[(usize, usize, f64)],
        observed: &[usize], minimum: usize, prior: &[f64],
    ) -> f64 {
        let mut multiplicity = observed.to_vec();
        self.retained.fill(0.0);
        for plan in plans {
            let root = plan.tree.root;
            fill_count_messages(plan, observed, &mut self.baseline);
            self.retained[0] += root_retention(&self.baseline[root], minimum, prior);
            for (index, &(donor, recipient, _)) in contacts.iter().enumerate() {
                if observed[recipient] == 0 {
                    continue;
                }
                multiplicity[donor] += 1;
                multiplicity[recipient] = 0;
                let path = &self.paths[index];
                // Only ancestors of the replaced observations can see a new
                // count message. Other subtrees keep the original mask mass.
                for &node in &path.order {
                    self.changed[node] = count_node(plan, node, &multiplicity, |child| {
                        if path.affected[child] { &self.changed[child] } else { &self.baseline[child] }
                    });
                }
                self.retained[index + 1] += root_retention(&self.changed[root], minimum, prior);
                multiplicity[donor] -= 1;
                multiplicity[recipient] = 1;
            }
        }
        let mut sum = self.retained[0];
        for (index, &(_, recipient, intensity)) in contacts.iter().enumerate() {
            sum += intensity * self.retained[if observed[recipient] == 0 { 0 } else { index + 1 }];
        }
        sum
    }
}

fn choose_weight(weights: &[f64], rng: &mut crate::rng::Rng) -> Option<usize> {
    let total: f64 = weights.iter().sum();
    if !total.is_finite() || total <= 0.0 {
        return None;
    }
    let mut threshold = total * rng.next_f64();
    let mut last = None;
    for (index, &weight) in weights.iter().enumerate() {
        if weight > 0.0 {
            last = Some(index);
            if threshold < weight {
                return last;
            }
            threshold -= weight;
        }
    }
    last
}

fn retained_assignment(
    plan: &PruningPlan<'_>, messages: &[StateCounts], prior: &[f64],
    minimum: usize, rng: &mut crate::rng::Rng,
) -> Option<Vec<i64>> {
    let root = plan.tree.root;
    let mut weights = [0.0; 2 * COUNT_BUCKETS];
    for state in 0..2 {
        for bucket in 0..COUNT_BUCKETS {
            if minimum == 0 || bucket >= ZERO_FLAG + minimum {
                weights[state * COUNT_BUCKETS + bucket] = prior[state] * messages[root][state][bucket];
            }
        }
    }
    let first = choose_weight(&weights, rng)?;
    let mut pending = vec![(root, first / COUNT_BUCKETS, first % COUNT_BUCKETS)];
    let mut leaves = vec![0; plan.tree.n_tips()];
    while let Some((node, state, mut bucket)) = pending.pop() {
        if let Some(taxon) = plan.tree.nodes[node].taxon {
            leaves[taxon] = state as i64;
            continue;
        }
        let children = &plan.tree.nodes[node].children;
        let mut prefix = Vec::with_capacity(children.len() + 1);
        let mut first = [0.0; COUNT_BUCKETS];
        first[0] = 1.0;
        prefix.push(first);
        for (index, &child) in children.iter().enumerate() {
            prefix.push(combine_counts(&prefix[index], &child_counts(plan, child, state, &messages[child])));
        }
        // Draw the last child's state and counts conditional on the aggregate.
        // Its remaining prefix count becomes the constraint for the preceding child.
        for (index, &child) in children.iter().enumerate().rev() {
            let mut choices = [0.0; 2 * COUNT_BUCKETS * COUNT_BUCKETS];
            let transition = &plan.transitions[child];
            for child_state in 0..2 {
                for a in 0..COUNT_BUCKETS {
                    for b in 0..COUNT_BUCKETS {
                        if combine_bucket(a, b) == bucket {
                            let slot = (child_state * COUNT_BUCKETS + a) * COUNT_BUCKETS + b;
                            choices[slot] = prefix[index][a] * transition.at(state, child_state)
                                * messages[child][child_state][b];
                        }
                    }
                }
            }
            let chosen = choose_weight(&choices, rng)?;
            let child_state = chosen / (COUNT_BUCKETS * COUNT_BUCKETS);
            pending.push((child, child_state, chosen % COUNT_BUCKETS));
            bucket = (chosen / COUNT_BUCKETS) % COUNT_BUCKETS;
        }
    }
    Some(leaves)
}

fn observation_multiplicity(mask: &[bool], contact: Option<&(usize, usize, f64)>) -> Vec<usize> {
    let mut counts: Vec<_> = mask.iter().map(|&observed| usize::from(observed)).collect();
    if let Some(&(donor, recipient, _)) = contact {
        counts[donor] += counts[recipient];
        counts[recipient] = 0;
    }
    counts
}

/// Exact conditional simulation. Select a rate and copy component by its
/// retained mass, then sample the inherited states conditional on retention.
/// Requests sharing a mask and component reuse the same upward messages.
pub fn binary_conditioned_draws(
    tree: &Tree, model: &dyn SubstModel, masks: &[Vec<bool>], prior: &[f64],
    rates: &[f64], minimum: usize, contacts: &[(usize, usize, f64)],
    rng: &mut crate::rng::Rng,
) -> Result<Vec<Vec<i64>>, &'static str> {
    use std::collections::{BTreeMap, HashMap};

    let plans: Vec<_> = rates.iter().map(|rate| PruningPlan::new(tree, model, *rate)).collect();
    let components = contacts.len() + 1;
    let mut mixtures = HashMap::<Vec<bool>, Vec<f64>>::new();
    let mut requests = BTreeMap::<(Vec<bool>, usize), Vec<usize>>::new();
    for (index, mask) in masks.iter().enumerate() {
        if !mixtures.contains_key(mask) {
            let work = (mixtures.len() + 1).saturating_mul(tree.nodes.len())
                .saturating_mul(plans.len()).saturating_mul(components);
            if work > 20_000_000 {
                return Err("conditional simulation exceeds the mask and component budget");
            }
            let mut weights = Vec::with_capacity(plans.len() * components);
            for plan in &plans {
                for component in 0..components {
                    let contact = component.checked_sub(1).map(|edge| &contacts[edge]);
                    let intensity = contact.map_or(1.0, |edge| edge.2);
                    let mass = if minimum == 0 || intensity == 0.0 {
                        1.0
                    } else {
                        retained_probability(plan, &observation_multiplicity(mask, contact), minimum, prior)
                    };
                    weights.push(mass * intensity);
                }
            }
            mixtures.insert(mask.clone(), weights);
        }
        let component = choose_weight(&mixtures[mask], rng)
            .ok_or("selection event has zero finite mass")?;
        requests.entry((mask.clone(), component)).or_default().push(index);
    }
    let mut output = vec![Vec::new(); masks.len()];
    for ((mask, component), indices) in requests {
        let plan = &plans[component / components];
        let contact = (component % components).checked_sub(1).map(|edge| &contacts[edge]);
        let messages = count_messages(plan, &observation_multiplicity(&mask, contact));
        for index in indices {
            let mut leaves = retained_assignment(plan, &messages, prior, minimum, rng)
                .ok_or("conditional assignment has zero finite mass")?;
            if let Some(&(donor, recipient, _)) = contact {
                leaves[recipient] = leaves[donor];
            }
            for (value, observed) in leaves.iter_mut().zip(&mask) {
                if !observed {
                    *value = -1;
                }
            }
            output[index] = leaves;
        }
    }
    Ok(output)
}

fn log_mean_exp(log_probs: &[f64]) -> f64 {
    let max = log_probs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if max == f64::NEG_INFINITY {
        return max;
    }
    let total: f64 = log_probs.iter().map(|p| (p - max).exp()).sum();
    max + (total / log_probs.len() as f64).ln()
}

/// Equal-weight site rates and a finite contact mixture with baseline weight 1.
/// At most one edge copies per character. Ascertainment conditions the mixture
/// after summing both the pattern mass and selection mass over components.
pub fn binary_mixture_scores(
    tree: &Tree,
    model: &dyn SubstModel,
    columns: &[Vec<Option<usize>>],
    root_prior: &[f64],
    rate_categories: &[f64],
    min_present: usize,
    contacts: &[(usize, usize, f64)],
) -> Vec<f64> {
    use std::collections::HashMap;

    let plans: Vec<_> = rate_categories
        .iter()
        .map(|rate| PruningPlan::new(tree, model, *rate))
        .collect();
    let mut selection = if min_present > 0 && !contacts.is_empty() {
        plans.first().map(|plan| ContactSelectionWorkspace::new(tree, &plan.order, contacts))
    } else {
        None
    };
    let mut normalizers = HashMap::new();
    let mut mixtures = vec![0.0; plans.len()];
    let weight_sum = 1.0 + contacts.iter().map(|(_, _, w)| w).sum::<f64>();
    let mut working_column = vec![None; tree.n_tips()];
    let mut partials = vec![[0.0; 2]; tree.nodes.len()];
    let mut components = vec![0.0; contacts.len() + 1];
    let mut score = |column: &[Option<usize>], mixtures: &mut [f64]| {
        for (result, plan) in mixtures.iter_mut().zip(&plans) {
            *result = plan.score_binary(column, root_prior, &mut partials);
        }
        log_mean_exp(mixtures)
    };
    columns
        .iter()
        .map(|column| {
            let log_normalizer = if min_present == 0 {
                weight_sum.ln()
            } else {
                let present = column.iter().filter(|s| **s == Some(1)).count();
                if present < min_present || !column.contains(&Some(0)) {
                    return f64::NEG_INFINITY;
                }
                let observed: Vec<_> = column
                    .iter()
                    .map(|state| usize::from(state.is_some()))
                    .collect();
                *normalizers.entry(observed.clone()).or_insert_with(|| {
                    let sum = if let Some(workspace) = selection.as_mut() {
                        workspace.mixture_mass(&plans, contacts, &observed, min_present, root_prior)
                    } else {
                        plans.iter().map(|plan| {
                            retained_probability(plan, &observed, min_present, root_prior)
                        }).sum::<f64>()
                    };
                    (sum / plans.len() as f64).ln()
                })
            };
            if !log_normalizer.is_finite() {
                return f64::NEG_INFINITY;
            }
            let base = score(column, &mut mixtures);
            components[0] = base;
            for (&(donor, recipient, intensity), component) in
                contacts.iter().zip(&mut components[1..])
            {
                let log_mass = match (column[donor], column[recipient]) {
                    (_, None) => base,
                    (Some(d), Some(r)) if d != r => f64::NEG_INFINITY,
                    (_, observed) => {
                        working_column.copy_from_slice(column);
                        working_column[donor] = observed;
                        working_column[recipient] = None;
                        score(&working_column, &mut mixtures)
                    }
                };
                *component = log_mass + intensity.ln();
            }
            let log_weighted_sum = log_mean_exp(&components) + (components.len() as f64).ln();
            log_weighted_sum - log_normalizer
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::subst::BinaryGainLoss;

    #[test]
    fn cached_contact_selection_preserves_exact_multirate_mask_mass() {
        let mut tree = Tree::new();
        let root = tree.add_node(None, 0.0, None);
        let left = tree.add_node(Some(root), 0.6, None);
        let right = tree.add_node(Some(root), 0.2, None);
        for (taxon, ancestor) in [left, left, right, right].into_iter().enumerate() {
            tree.add_node(Some(ancestor), 0.7 + taxon as f64 / 3.0, Some(taxon));
        }
        tree.root = root;
        let model = BinaryGainLoss::new(0.4, 0.8);
        let prior = model.stationary();
        let plans: Vec<_> = [0.3, 1.0, 2.0].into_iter()
            .map(|rate| PruningPlan::new(&tree, &model, rate)).collect();
        let contacts = [(0, 1, 0.2), (2, 0, 0.4), (1, 3, 0.03)];
        let mut workspace = ContactSelectionWorkspace::new(&tree, &plans[0].order, &contacts);
        for bits in 0..16 {
            let mask: Vec<_> = (0..4).map(|index| bits & (1 << index) != 0).collect();
            for minimum in [1, 2] {
                let mass = |contact| plans.iter().map(|plan| {
                    retained_probability(plan, &observation_multiplicity(&mask, contact), minimum, &prior)
                }).sum::<f64>();
                let mut expected = mass(None);
                for contact in &contacts {
                    expected += contact.2 * mass(Some(contact));
                }
                let observed = observation_multiplicity(&mask, None);
                assert_eq!(workspace.mixture_mass(&plans, &contacts, &observed, minimum, &prior).to_bits(),
                           expected.to_bits());
            }
        }
    }

    /// Build a two-tip tree: root -> {tipA, tipB}, each branch length t.
    fn two_tip(t: f64) -> Tree {
        let mut tree = Tree::new();
        let root = tree.add_node(None, 0.0, None);
        tree.add_node(Some(root), t, Some(0));
        tree.add_node(Some(root), t, Some(1));
        tree.root = root;
        tree
    }

    #[test]
    fn two_tip_matches_hand_computation() {
        // With both tips present (state 1), the likelihood is
        // sum_r pi_r P_r1(t) P_r1(t).
        let model = BinaryGainLoss::new(0.4, 0.8);
        let t = 0.6;
        let tree = two_tip(t);
        let pi = model.stationary();
        let p = model.transition(t);
        let expected = pi[0] * p.at(0, 1) * p.at(0, 1) + pi[1] * p.at(1, 1) * p.at(1, 1);
        let states = vec![Some(1usize), Some(1usize)];
        let tips = TipStates { states: &states };
        let got = log_likelihood(&tree, &model, &tips, &pi);
        assert!(
            (got - expected.ln()).abs() < 1e-10,
            "{got} vs {}",
            expected.ln()
        );
    }

    #[test]
    fn missing_differs_from_absent() {
        let model = BinaryGainLoss::new(0.4, 0.8);
        let t = 0.6;
        let tree = two_tip(t);
        let pi = model.stationary();

        let absent = vec![Some(1usize), Some(0usize)];
        let missing = vec![Some(1usize), None];
        let ll_absent = log_likelihood(&tree, &model, &TipStates { states: &absent }, &pi);
        let ll_missing = log_likelihood(&tree, &model, &TipStates { states: &missing }, &pi);
        // Missing marginalizes over tip B, so its likelihood is strictly larger
        // than pinning B to a specific state.
        assert!(ll_missing > ll_absent, "{ll_missing} !> {ll_absent}");
    }

    #[test]
    fn scaling_stable_on_wide_polytomy() {
        let mut tree = Tree::new();
        let root = tree.add_node(None, 0.0, None);
        let n = 2000;
        for taxon in 0..n {
            tree.add_node(Some(root), 100.0, Some(taxon));
        }
        let model = BinaryGainLoss::new(0.5, 0.5);
        let columns = vec![vec![Some(1); n], vec![Some(0); n]];
        let expected = -((2 * n) as f64) * 2.0f64.ln();
        let score = log_likelihood_block(&tree, &model, &columns, &[0.5, 0.5]);
        assert!((score - expected).abs() < 1e-8, "{score} vs {expected}");
    }

    #[test]
    fn scaling_stable_on_deep_tree() {
        // A caterpillar tree deep enough to underflow without scaling.
        let model = BinaryGainLoss::new(0.5, 0.5);
        let mut tree = Tree::new();
        let root = tree.add_node(None, 0.0, None);
        let mut current = root;
        let mut taxon = 0;
        for _ in 0..200 {
            let internal = tree.add_node(Some(current), 0.1, None);
            tree.add_node(Some(current), 0.1, Some(taxon));
            taxon += 1;
            current = internal;
        }
        tree.nodes[current].taxon = Some(taxon);
        tree.root = root;
        let pi = model.stationary();
        let states: Vec<Option<usize>> = (0..=taxon).map(|_| Some(1usize)).collect();
        let ll = log_likelihood(&tree, &model, &TipStates { states: &states }, &pi);
        assert!(ll.is_finite(), "deep-tree likelihood underflowed: {ll}");
    }
}
