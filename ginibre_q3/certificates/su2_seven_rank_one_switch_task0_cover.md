# SU(2) seven-factor rank-one switch task 0

Date: 2026-07-27

## Primary reduced-cover certificate

In the reduced 7,812-cell source snapshot, indices `[0,43)` are
exactly the endpoint task

```text
orbit=1 parity=0 selected=1 position=0 kind=2 rank=1.
```

Command:

```text
Q3_MAX_THREADS=1
verify_su2_seven_shallow_z3 --small-switch-cover-range 0 43
```

Source, executable, and transcript SHA-256:

```text
cdeba6c1beffc977f9d7df13e25ab734dec6b4186b5291eab008f9e518c8bbbd
  verify_su2_seven_shallow_z3.cpp
42ead492a4fbd3257c0b0af5b0d49ee85a940b15508851b9b2ebf00a6c01a3b8
  verify_su2_seven_residual_z3.cpp
5367dab66cc311756dd17602bc9bb3d5c2535f862c1752c8859f0324f4d8c64f
  verify_su2_seven_shallow_z3
c436ef0a89201e0f7bddb5643ba86ba944a9f1f2e11251999ff51fbc38c8c54a
  su2_rank_one_switch_cover_7812_task0.log
```

The terminal transcript is

```text
SU2_SEVEN_SHALLOW_SMALL_SWITCH_COVER_Z3
  tasks=43 counterexamples=UNSAT result=PASS
```

This complete shard proves the task.  The older overlapping shards
below are retained as an independent overcomplete check.

This certificate binds two exact one-worker shards of the original
94,212-cell rank-one switch cover.  At that source snapshot, indices
`[0,523)` are exactly the endpoint task

```text
orbit=1 parity=0 selected=1 position=0 kind=2 rank=1
```

and its local/channel subdivisions.

## Prefix shard

Command:

```text
Q3_MAX_THREADS=1
verify_su2_seven_shallow_z3 --small-switch-cover-range 0 523
```

Source and executable SHA-256:

```text
7f855ef8d6dc9c5326555c44952f71042c450f1394df217226acb1f29fe409b9
  verify_su2_seven_shallow_z3.cpp
e54a025fd3154159e2951df492bfb8a1001b5bad5c9896adc8af09708aef8e08
  verify_su2_seven_residual_z3.cpp
d3cf960ca8d1bfd08d3f42633657561fdabfcae235ba92b0082841123f1be11d
  verify_su2_seven_shallow_z3
```

The live source-bound transcript reached

```text
SU2_SEVEN_SHALLOW_SMALL_SWITCH_COVER_Z3
  progress=38/523 assigned=39 failed=0
```

with one worker.  The worker obtains indices in increasing order and
increments `completed` only after an `UNSAT` result.  This certifies
every cell in `[0,38)`.  The still-running continuation is not used.

## Suffix shard

Command:

```text
Q3_MAX_THREADS=1
verify_su2_seven_shallow_z3 --small-switch-cover-range 35 523
```

Source and executable SHA-256:

```text
7f855ef8d6dc9c5326555c44952f71042c450f1394df217226acb1f29fe409b9
  verify_su2_seven_shallow_z3.cpp
6a2f17f7b43300b7314f6a85790e39d2dd6e4080551a7945ef2940271de21485
  verify_su2_seven_residual_z3.cpp
f791143573d92931f1e98785a0c3fcbe90792fb464ab4881ea1fb531beac3a7a
  verify_su2_seven_shallow_z3
```

The second residual source differs only by exact endpoint consequences
and the selected-simple-current local formula of Lemma 5A7B30.  Its
terminal transcript is

```text
SU2_SEVEN_SHALLOW_SMALL_SWITCH_COVER_Z3
  tasks=488 counterexamples=UNSAT result=PASS
```

This certifies every cell in `[35,523)`.

## Independent overlapping conclusion

The two certified intervals overlap and their union is `[0,523)`.
Therefore all cells of the displayed rank-one endpoint task are
`UNSAT`.  Lemmas 5A7B28--31 imply

```text
|B_C|<=L(C)-1
```

throughout that task, so its seven-factor direct-payment inequality
holds uniformly for every level.
