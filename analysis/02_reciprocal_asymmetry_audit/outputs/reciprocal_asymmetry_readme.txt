Ovule project reciprocal-cross family audit
============================================

Definitions used in this audit:
- Retained row: Table1 remove == no.
- Reciprocal interspecific pair: both interspecific_value_sp1_sp2 and interspecific_value_sp2_sp1 are numeric.
- Directional RI pair: reciprocal interspecific pair plus numeric positive intraspecific values for both species.
- Asymmetry model-ready pair: directional RI pair plus positive ovule_number for both species in Table2.
- High-confidence pair: asymmetry model-ready pair with both confidence values <= 2.
- Strict PI/K2P subset: asymmetry model-ready pair with nonmissing K2P distance. In these files this also maps to the PI subset.

Retained pair rows: 560 pairs across 50 families
Both reciprocal interspecific directions: 396 pairs across 46 families
Directional RI possible: reciprocal + both intraspecific controls: 383 pairs across 45 families
Asymmetry model ready: directional RI + both ovule numbers: 305 pairs across 41 families
Asymmetry model ready with nonzero ovule difference: 217 pairs across 27 families
High-confidence ovules only: confidence <= 2: 213 pairs across 31 families
Strict PI/K2P asymmetry subset: 63 pairs across 35 families
Strict PI/K2P + high-confidence ovules: 42 pairs across 27 families

Family-level modelability counts:

reciprocal interspecific pairs:
  families with >= 1: 46
  families with >= 2: 31
  families with >= 3: 26
  families with >= 5: 16
  families with >= 10: 10
  families with >= 20: 4

directional RI pairs:
  families with >= 1: 45
  families with >= 2: 29
  families with >= 3: 24
  families with >= 5: 16
  families with >= 10: 10
  families with >= 20: 4

asymmetry pairs with ovules:
  families with >= 1: 41
  families with >= 2: 26
  families with >= 3: 21
  families with >= 5: 11
  families with >= 10: 4
  families with >= 20: 3

asymmetry pairs with nonzero ovule difference:
  families with >= 1: 27
  families with >= 2: 18
  families with >= 3: 13
  families with >= 5: 6
  families with >= 10: 2
  families with >= 20: 2

high-confidence asymmetry pairs:
  families with >= 1: 31
  families with >= 2: 19
  families with >= 3: 15
  families with >= 5: 8
  families with >= 10: 3
  families with >= 20: 3

strict PI/K2P asymmetry pairs:
  families with >= 1: 35
  families with >= 2: 14
  families with >= 3: 5
  families with >= 5: 2
  families with >= 10: 0
  families with >= 20: 0
