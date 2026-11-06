# Dated sampler pilots

Each pilot used the same pinned dated IE-CoR target, four independent starts,
1,000 warmup iterations, 3,000 sampling iterations and thinning by ten.
The replica-exchange arm used eight replicas per start. Only the cold draws
entered these summaries.

| Arm | Root split R-hat | Root bulk ESS | Distinct topologies by chain |
| --- | ---: | ---: | --- |
| Control | 2.115 | 5.42 | 25, 26, 20, 23 |
| Native attachment-age SPR | 2.155 | 5.68 | 28, 22, 34, 26 |
| Lexical replica exchange | 1.990 | 6.37 | 21, 26, 24, 30 |
| Informed fixed-root attachment | 2.917 | 5.67 | 63, 58, 92, 77 |
| Informed root-stem attachment | 2.529 | 5.17 | 65, 57, 91, 68 |

Every arm had a maximum between-chain clade frequency difference of 1.00.
No replica made a round trip. Both informed arms increased the number of
distinct topologies visited, while root-age agreement worsened relative to
the control. These results leave the default sampler and declared full
thresholds in place.
[Machine-readable pilot results](sampler-pilots.json)
[Source-bound trajectories and checkpoints](https://github.com/Aneesh495/diachron/releases/tag/pilots-47f25ca)

## Calendar-time mixing with fixed topology

At `c6b81929cccb2f6e3c1591090627f9f47b6226a3`, two prescribed topologies
each ran four chains on the observed dated target. The roots used the first
retained trees of chains one and four in the full observed run. Each
conditional chain stayed within its prescribed topology and replayed.

| Anchor | Root split R-hat | Root bulk ESS | Root tail ESS |
| --- | ---: | ---: | ---: |
| Chain one | 1.458 | 8.81 | 21.61 |
| Chain four | 1.182 | 15.77 | 42.34 |

Calendar mixing remained poor on both fixed topologies at the pilot budget.
The results motivate checking coupled clock and dated-height moves separately
from topology transitions. [Conditional diagnostic](conditional-time-c6b8192.json)

## Coupled calendar scale

Joint scaling changed within-topology root bulk ESS from 8.81 to 33.77 on
the first anchor and from 15.77 to 6.50 on the second. All calendar-scale
proposals were accepted on the fixed trees. The free-tree pilot had root
split R-hat 2.079, bulk ESS 5.46, and maximum between-chain clade difference
1.00. The fixed and free trajectories passed replay.
[Calendar-scale pilot measurements](calendar-scale-3b6836f.json)
[Raw pilot archive](https://github.com/Aneesh495/diachron/releases/tag/pilots-47f25ca)

## Gradient-guided joint time and clock move

The fixed-tree root bulk ESS values were 12.69 and 6.61, with split R-hat
1.251 and 1.646. The move accepted 708 of 1,386 and 646 of 1,473 attempts,
respectively. The free-tree pilot improved root split R-hat from 2.115 in
the control to 1.580, with bulk ESS 7.44 and maximum between-chain clade
difference 1.00. The four free chains took 75.2 seconds. The trajectories
passed replay and the convergence thresholds remain unmet.
[Joint-move pilot measurements](time-block-c8d6788.json)
[Raw pilot archive](https://github.com/Aneesh495/diachron/releases/tag/pilots-47f25ca)

## Whole-backbone parsimony guide

At `c2154c34b2676f9dcc8206ab082077852345156e`, the four-chain guide-only
pilot visited 48 to 81 distinct topologies per chain. It reached root split
R-hat 2.779 and bulk ESS 5.16 in 107.7 seconds. The arm combining the guide
with joint time proposals reached R-hat 1.751 and bulk ESS 6.33 in 138.6
seconds. Both arms reached maximum between-chain clade difference 1.00 and
passed saved-draw replay.
[Whole-backbone pilot measurements](whole-backbone-c2154c3.json)
[Raw pilot archive](https://github.com/Aneesh495/diachron/releases/tag/pilots-47f25ca)

## Four-step Hamiltonian time move

At `d1f6444ddb5918cd5738f76bfe1db71c4d03a7d1`, the two fixed-topology
anchors reached root split R-hat 1.036 and 1.103 and root bulk ESS 79.39 and
27.47. Both improved on the corresponding conditional controls. The free-tree
fit took 176.0 seconds and reached root split R-hat 1.336, bulk ESS 9.65,
tail ESS 95.69, and maximum between-chain clade difference 1.00. Its four
chains visited 31, 27, 25, and 37 topologies. All three fits passed replay.
[Hamiltonian pilot measurements](hamiltonian-d1f6444.json)
[Raw Hamiltonian archive](https://github.com/Aneesh495/diachron/releases/download/pilots-47f25ca/hamiltonian_d1f6444.tar.gz)

## Lexical likelihood guide with Hamiltonian time moves

At `70467e9a790c3ba1aae73eb2b3434421f1cf7c48`, likelihood-guided placement
with four-step Hamiltonian time moves took 203.8 seconds. The chains visited
46, 51, 56, and 37 topologies. Root split R-hat reached 1.324, with bulk ESS
9.57 and tail ESS 11.86. The maximum clade difference remained 1.00. The
trajectories passed replay. The guide increased the number of topologies
visited in three chains but did not resolve their disagreement. These pilots
did not trigger a new full campaign or change its declared thresholds.
[Lexical placement measurements](likelihood-placement-70467e9.json)
[Raw lexical placement archive](https://github.com/Aneesh495/diachron/releases/download/pilots-47f25ca/lexical_placement_70467e9.tar.gz)
