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
