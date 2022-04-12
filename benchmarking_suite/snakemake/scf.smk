
''' SCF and SCF-derived quantities' workflow.'''

locstr = "raw/{struct}_ch{charge}_mult{mult}_{mftype}_{xc}_zeta{basis}_df{dfint}_ecp{ecp}"
rule runscf:
  output:
    out = locstr + "/scfcalc.chk",
  params:
    loc = locstr,
    json = locstr + "/scfcalc.in.json",
  run: 

    jsondump(dict(
        struct_source = relpath(find_structure(wildcards.struct, int(wildcards.charge), int(wildcards.mult)), params['loc']),
        structargs = dict(
            charge        =  int(wildcards.charge),
            multiplicity  =  int(wildcards.mult),
            basis         =  wildcards.basis,
            ecp           =  wildcards.ecp,
        ),
        mfargs=dict(
          mftype  = wildcards.mftype,
          xc      = wildcards.xc,
          usedf   = bool(int(wildcards.dfint)),
        ),
      ),params['json'])

    subscript("scfcalc.py", params['loc'], time='1-0', wait=True, ptype='rome', **PYSUB)

locstr = "raw/{struct}_latc{latc}_ch{charge}_mult{mult}_{mftype}_{xc}_zeta{basis}_df{dfint}_ecp{ecp}"
rule runscf_pbc:
  output:
    out = locstr + "/scfcalc.chk",
  params:
    loc = locstr,
    json = locstr + "/scfcalc.in.json",
  run: 

    jsondump(dict(
        struct_source = relpath(find_structure(wildcards.struct, int(wildcards.charge), int(wildcards.mult)), params['loc']),
        structargs = dict(
            charge        =  int(wildcards.charge),
            multiplicity  =  int(wildcards.mult),
            basis         =  wildcards.basis,
            ecp           =  wildcards.ecp,
            latc          =  float(wildcards.latc),
        ),
        mfargs=dict(
          mftype  = wildcards.mftype,
          xc      = wildcards.xc,
          usedf   = bool(int(wildcards.dfint)),
        ),
      ),params['json'])

    subscript("scfcalc.py", params['loc'], time='1-0', wait=True, **PYSUB)

ruleorder: runscf_pbc > runscf

default_structures = {
    ("VCp2", 0): "../geometries/xyz/VCp2_0_4.xyz",
    ("n2",0): "../geometries/xyz/n2.xyz",
    ("v",0): "../geometries/xyz/v.xyz",
  }

def find_structure(struct, charge, mult):
  best_match = f"../geometries/xyz/{struct}_{charge}_{mult}.xyz"
  if exists(best_match):
    return best_match
  else:
    return default_structures[struct, charge]
