''' MCSCF and MCSCF-derived quantites' workflow.'''

locstr = "{scfdir}/cisd"
rule runcisd:
  input:
    scfout = "{scfdir}/scfcalc.out.json",
  output:
    out = locstr + "/mcdata.h5",
  params:
    json = locstr + "/mccalc.in.json",
    loc = locstr,
  run:
    with open(input['scfout'],'r') as inpf:
      scfout = json.load(inpf)

    jsondump(dict(
        chkfile = relpath(f"{wildcards.scfdir}/{scfout['chkfile']}", params['loc']),
        mctype = 'cisd',
        mcargs = {},
        nroots  = 13,
      ), params['json'])

    subscript("mccalc.py", params['loc'], time='1-0', wait=True, **PYSUB)

locstr = "{scfdir}/ccsd"
rule runccsd:
  input:
    scfout = "{scfdir}/scfcalc.out.json",
  output:
    out = locstr + "/mcdata.h5",
  params:
    json = locstr + "/mccalc.in.json",
    loc = locstr,
  run:
    with open(input['scfout'],'r') as inpf:
      scfout = json.load(inpf)

    jsondump(dict(
        chkfile = relpath(f"{wildcards.scfdir}/{scfout['chkfile']}", params['loc']),
        mctype = 'ccsd',
        mcargs = {},
        nroots  = 13,
      ), params['json'])

    subscript("mccalc.py", params['loc'], time='1-0', wait=True, **PYSUB)

locstr = "{scfdir}/casci_{nelecas}_{ncas}"
rule runcasci:
  input:
    scfout = "{scfdir}/scfcalc.out.json",
  output:
    out = locstr + "/mcdata.h5",
  params:
    json = locstr + "/mccalc.in.json",
    loc = locstr,
  run:
    with open(input['scfout'],'r') as inpf:
      scfout = json.load(inpf)

    jsondump(dict(
        chkfile = relpath(f"{wildcards.scfdir}/{scfout['chkfile']}", params['loc']),
        mctype = 'casci',
        mcargs = dict(
            nelecas = int(wildcards.nelecas), 
            ncas = int(wildcards.ncas),
          ),
        nroots  = 13,
      ), params['json'])

    subscript("mccalc.py", params['loc'], time='0-1', wait=True, **PYSUB)

locstr = "{scfdir}/casci_avas"
rule runcasci_avas:
  input:
    scfout = "{scfdir}/scfcalc.out.json",
  output:
    out = locstr + "/mcdata.h5",
  params:
    json = locstr + "/mccalc.in.json",
    loc = locstr,
  run:
    with open(input['scfout'],'r') as inpf:
      scfout = json.load(inpf)

    jsondump(dict(
        chkfile = relpath(f"{wildcards.scfdir}/{scfout['chkfile']}", params['loc']),
        mctype = 'casci',
        mcargs = dict(use_avas = True),
        nroots  = 13,
      ), params['json'])

    subscript("mccalc.py", params['loc'], time='0-1', wait=True, **PYSUB)

locstr = "{scfdir}/casscf_{nelecas}_{ncas}_state{state}_nr{nroots}"
rule runcasscf:
  input:
    scfout = "{scfdir}/scfcalc.out.json",
  output:
    out = locstr + "/mcdata.h5",
  params:
    json = locstr + "/mccalc.in.json",
    loc = locstr,
  run:
    with open(input['scfout'],'r') as inpf:
      scfout = json.load(inpf)

    jsondump(dict(
        chkfile = relpath(f"{wildcards.scfdir}/{scfout['chkfile']}", params['loc']),
        mctype = 'casscf',
        mcargs = dict(
            nelecas = int(wildcards.nelecas), 
            ncas = int(wildcards.ncas),
            state = wildcards.state, # average or specific.
          ),
        nroots = int(wildcards.nroots),
      ), params['json'])

    subscript("mccalc.py", params['loc'], time='1-0', wait=True, **PYSUB)

locstr = "{scfdir}/casscf_avas{thresh}_state{state}_nr{nroots}"
rule runcasscf_avas:
  input:
    scfout = "{scfdir}/scfcalc.out.json",
  output:
    out   = locstr + "/mcdata.h5",
  params:
    json = locstr + "/mccalc.in.json",
    loc = locstr,
  run:
    with open(input['scfout'],'r') as inpf:
      scfout = json.load(inpf)

    jsondump(dict(
        chkfile = relpath(f"{wildcards.scfdir}/{scfout['chkfile']}", params['loc']),
        mctype = 'casscf',
        mcargs = dict(
            use_avas = True,
            avas_thresh = float(wildcards.thresh),
            state = wildcards.state, # average or specific.
          ),
        nroots  = int(wildcards.nroots),
      ), params['json'])

    subscript("mccalc.py", params['loc'], time='1-0', wait=True, **PYSUB)

locstr = "{scfdir}"
rule meausure_scf_observables:
  input:
    scfout = locstr + "/scfcalc.out.json",
    iaohdfn = locstr + "/iao_pzbasis.h5",
  output:
    out = locstr + "/scf_observables.out.json"
  params:
    json = locstr + "/scf_observables.in.json",
    loc = locstr,
  run:
    with open(input['scfout'],'r') as inpf:
      scfout = json.load(inpf)

    jsondump(dict(
        chkfile = relpath(f"{wildcards.scfdir}/{scfout['chkfile']}", params['loc']),
        basishdfn = relpath(input['iaohdfn'], params['loc']),
      ), params['json'])

    subscript("measure_observables.py", params['loc'], time='0-1', wait=True, **PYSUB)

locstr = "{scfdir}/{mcdir}"
rule measure_mc_observables:
  input:
    mcdata = locstr + "/mcdata.h5",
    basis  = "{scfdir}/iao_pzbasis.h5",
    scfout  = "{scfdir}/scfcalc.out.json",
  output:
    out = locstr + "/observables.out.json"
  params:
    json = locstr + "/observables.in.json",
    loc = locstr,
  run: 
    with open(input['scfout'],'r') as inpf:
      scfout = json.load(inpf)

    jsondump(dict(
        chkfile = relpath(f"{wildcards.scfdir}/{scfout['chkfile']}", params['loc']),
        basishdfn = relpath(f"{wildcards.scfdir}/iao_pzbasis.h5", params['loc']),
        mcdata = "mcdata.h5",
      ) ,params['json'])

    subscript("measure_observables.py", params['loc'], time='0-1', wait=True, **PYSUB)
