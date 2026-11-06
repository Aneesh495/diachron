//! Reverse tree messages for derivatives with respect to branch durations.

use std::collections::HashMap;

use super::{
    child_counts, combine_bucket, combine_counts, fill_count_messages, root_retention, Counts,
    PruningPlan, StateCounts, COUNT_BUCKETS, ZERO_FLAG,
};
use crate::subst::{BinaryGainLoss, SubstModel};
use crate::tree::Tree;

type Pair = [f64; 2];

fn multiply(a: Pair, b: Pair) -> Pair {
    let product = [a[0] * b[0], a[1] * b[1]];
    let scale = product[0].max(product[1]);
    if scale > 0.0 {
        [product[0] / scale, product[1] / scale]
    } else {
        product
    }
}

struct BinaryWorkspace {
    inside: Vec<Pair>,
    outside: Vec<Pair>,
    suffix: Vec<Pair>,
}

impl BinaryWorkspace {
    fn new(n: usize) -> Self {
        Self {
            inside: vec![[0.0; 2]; n],
            outside: vec![[0.0; 2]; n],
            suffix: Vec::new(),
        }
    }

    fn projected(&self, plan: &PruningPlan<'_>, child: usize) -> Pair {
        let p = &plan.transitions[child];
        let [a, b] = self.inside[child];
        [
            p.at(0, 0) * a + p.at(0, 1) * b,
            p.at(1, 0) * a + p.at(1, 1) * b,
        ]
    }

    fn score(
        &mut self,
        plan: &PruningPlan<'_>,
        column: &[Option<usize>],
        prior: Pair,
        change: &[Pair],
        gradient: &mut [f64],
    ) -> f64 {
        let value = plan.score_binary(column, &prior, &mut self.inside);
        gradient.fill(0.0);
        if !value.is_finite() {
            return value;
        }
        self.outside[plan.tree.root] = prior;
        for &node in plan.order.iter().rev() {
            let children = &plan.tree.nodes[node].children;
            self.suffix.resize(children.len() + 1, [1.0, 1.0]);
            self.suffix[children.len()] = [1.0, 1.0];
            for (i, &child) in children.iter().enumerate().rev() {
                self.suffix[i] = multiply(self.suffix[i + 1], self.projected(plan, child));
            }
            let mut prefix = self.outside[node];
            for (i, &child) in children.iter().enumerate() {
                let siblings = multiply(prefix, self.suffix[i + 1]);
                let p = &plan.transitions[child];
                let outside = [
                    siblings[0] * p.at(0, 0) + siblings[1] * p.at(1, 0),
                    siblings[0] * p.at(0, 1) + siblings[1] * p.at(1, 1),
                ];
                let [a, b] = self.inside[child];
                let denominator = outside[0] * a + outside[1] * b;
                let [gain, loss] = change[child];
                gradient[child] = (b - a) * (siblings[0] * gain - siblings[1] * loss) / denominator;
                self.outside[child] = multiply(outside, [1.0, 1.0]);
                prefix = multiply(prefix, self.projected(plan, child));
            }
        }
        value
    }
}

struct SelectionWorkspace {
    inside: Vec<StateCounts>,
    outside: Vec<StateCounts>,
    suffix: Vec<StateCounts>,
}

impl SelectionWorkspace {
    fn new(n: usize) -> Self {
        let zero = [[0.0; COUNT_BUCKETS]; 2];
        Self {
            inside: vec![zero; n],
            outside: vec![zero; n],
            suffix: Vec::new(),
        }
    }

    fn score(
        &mut self,
        plan: &PruningPlan<'_>,
        multiplicity: &[usize],
        minimum: usize,
        prior: Pair,
        change: &[Pair],
        gradient: &mut [f64],
    ) -> f64 {
        fill_count_messages(plan, multiplicity, &mut self.inside);
        let root = plan.tree.root;
        let mass = root_retention(&self.inside[root], minimum, &prior);
        gradient.fill(0.0);
        if mass <= 0.0 {
            return f64::NEG_INFINITY;
        }
        let mut unit = [0.0; COUNT_BUCKETS];
        unit[0] = 1.0;
        for (state, p) in prior.into_iter().enumerate() {
            self.outside[root][state].fill(0.0);
            self.outside[root][state][ZERO_FLAG + minimum..].fill(p);
        }
        for &node in plan.order.iter().rev() {
            let children = &plan.tree.nodes[node].children;
            self.suffix.resize(children.len() + 1, [unit; 2]);
            self.suffix[children.len()] = [unit; 2];
            let mut prefix = [unit; 2];
            for (i, &child) in children.iter().enumerate().rev() {
                for state in 0..2 {
                    let projected = child_counts(plan, child, state, &self.inside[child]);
                    self.suffix[i][state] = combine_counts(&projected, &self.suffix[i + 1][state]);
                }
            }
            for (i, &child) in children.iter().enumerate() {
                let p = &plan.transitions[child];
                let [gain, loss] = change[child];
                let mut outside = [[0.0; COUNT_BUCKETS]; 2];
                let mut derivative = 0.0;
                for (state, prefix_counts) in prefix.iter_mut().enumerate() {
                    let siblings = combine_counts(prefix_counts, &self.suffix[i + 1][state]);
                    let context = self.outside[node][state];
                    for bucket in 0..COUNT_BUCKETS {
                        let weight: f64 = siblings
                            .iter()
                            .enumerate()
                            .map(|(other, value)| value * context[combine_bucket(bucket, other)])
                            .sum();
                        outside[0][bucket] += p.at(state, 0) * weight;
                        outside[1][bucket] += p.at(state, 1) * weight;
                        let difference =
                            self.inside[child][1][bucket] - self.inside[child][0][bucket];
                        derivative += weight * difference * if state == 0 { gain } else { -loss };
                    }
                    let projected = child_counts(plan, child, state, &self.inside[child]);
                    *prefix_counts = combine_counts(prefix_counts, &projected);
                }
                self.outside[child] = outside;
                gradient[child] = derivative / mass;
            }
        }
        mass.ln()
    }
}

/// Stable sum of component masses and their probability-weighted derivatives.
struct Mixture {
    log_mass: f64,
    gradient: Vec<f64>,
}

impl Mixture {
    fn new(n: usize) -> Self {
        Self {
            log_mass: f64::NEG_INFINITY,
            gradient: vec![0.0; n],
        }
    }

    fn add(&mut self, log_mass: f64, gradient: &[f64]) {
        if log_mass == f64::NEG_INFINITY {
            return;
        }
        let peak = self.log_mass.max(log_mass);
        let old = (self.log_mass - peak).exp();
        let new = (log_mass - peak).exp();
        for (mean, derivative) in self.gradient.iter_mut().zip(gradient) {
            *mean = (*mean * old + derivative * new) / (old + new);
        }
        self.log_mass = peak + (old + new).ln();
    }
}

/// Return the total conditioned log likelihood and its branch-duration score.
pub fn binary_mixture_branch_gradient(
    tree: &Tree,
    model: &BinaryGainLoss,
    columns: &[Vec<Option<usize>>],
    rates: &[f64],
    minimum: usize,
    contacts: &[(usize, usize, f64)],
) -> Result<(f64, Vec<f64>), &'static str> {
    let n = tree.nodes.len();
    let stationary = model.stationary();
    let prior = [stationary[0], stationary[1]];
    let plans: Vec<_> = rates
        .iter()
        .map(|&rate| {
            let change: Vec<_> = tree
                .nodes
                .iter()
                .map(|node| {
                    let exposure = node.branch_len * rate;
                    let retain = (-(model.gain + model.loss) * exposure).exp();
                    [retain * model.gain * rate, retain * model.loss * rate]
                })
                .collect();
            (PruningPlan::new(tree, model, rate), change)
        })
        .collect();
    let mut patterns = BinaryWorkspace::new(n);
    let mut selection = SelectionWorkspace::new(n);
    let mut score = 0.0;
    let mut gradient = vec![0.0; n];
    let mut component_gradient = vec![0.0; n];
    let mut normalizers = HashMap::new();
    let unselected_mass =
        (rates.len() as f64).ln() + (1.0 + contacts.iter().map(|edge| edge.2).sum::<f64>()).ln();
    let mut working = vec![None; tree.n_tips()];

    for column in columns {
        if minimum > 0
            && (column.iter().filter(|s| **s == Some(1)).count() < minimum
                || !column.contains(&Some(0)))
        {
            return Err("column fails the retention condition");
        }
        let mut numerator = Mixture::new(n);
        let base_weight = 1.0
            + contacts
                .iter()
                .filter(|edge| column[edge.1].is_none())
                .map(|edge| edge.2)
                .sum::<f64>();
        for (plan, change) in &plans {
            let value = patterns.score(plan, column, prior, change, &mut component_gradient);
            numerator.add(value + base_weight.ln(), &component_gradient);
            for &(donor, recipient, intensity) in contacts {
                let Some(observed) = column[recipient] else {
                    continue;
                };
                if intensity == 0.0 || column[donor].is_some_and(|state| state != observed) {
                    continue;
                }
                working.copy_from_slice(column);
                working[donor] = Some(observed);
                working[recipient] = None;
                let value = patterns.score(plan, &working, prior, change, &mut component_gradient);
                numerator.add(value + intensity.ln(), &component_gradient);
            }
        }
        let denominator = if minimum == 0 {
            unselected_mass
        } else {
            let mask: Vec<_> = column
                .iter()
                .map(|state| usize::from(state.is_some()))
                .collect();
            let normalizer = normalizers.entry(mask).or_insert_with_key(|mask| {
                let mut mixture = Mixture::new(n);
                let base_weight = 1.0
                    + contacts
                        .iter()
                        .filter(|edge| mask[edge.1] == 0)
                        .map(|edge| edge.2)
                        .sum::<f64>();
                let mut changed = mask.clone();
                for (plan, change) in &plans {
                    let value = selection.score(
                        plan,
                        mask,
                        minimum,
                        prior,
                        change,
                        &mut component_gradient,
                    );
                    mixture.add(value + base_weight.ln(), &component_gradient);
                    for &(donor, recipient, intensity) in contacts {
                        if mask[recipient] == 0 || intensity == 0.0 {
                            continue;
                        }
                        changed[donor] += 1;
                        changed[recipient] = 0;
                        let value = selection.score(
                            plan,
                            &changed,
                            minimum,
                            prior,
                            change,
                            &mut component_gradient,
                        );
                        mixture.add(value + intensity.ln(), &component_gradient);
                        changed[donor] -= 1;
                        changed[recipient] = 1;
                    }
                }
                mixture
            });
            for (derivative, correction) in gradient.iter_mut().zip(&normalizer.gradient) {
                *derivative -= correction;
            }
            normalizer.log_mass
        };
        score += numerator.log_mass - denominator;
        for (derivative, increment) in gradient.iter_mut().zip(&numerator.gradient) {
            *derivative += increment;
        }
    }
    if !score.is_finite() || gradient.iter().any(|g| !g.is_finite()) {
        return Err("nonfinite binary likelihood or branch score");
    }
    Ok((score, gradient))
}
