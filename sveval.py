#!/usr/bin/python
import argparse
import cyvcf2
import sys
import shutil
import subprocess
import tempfile
import os
import json


def extract_truvari(f):

    name = "truvari"
    _, sizemin, sizemax, svtype = f.split("/")[-2].rsplit(".", 3)

    with open(f) as fh:
        d = json.load(fh)
        d["TP"] = d.pop("TP-call")
        d["name"] = name
        d["sizemin"] = sizemin
        d["sizemax"] = sizemax
        d["svtype"] = svtype
        return d

def extract_wittyer(path):

    with open(path) as fh:
        d = json.load(fh)
        for row in d["PerSampleStats"][0]["DetailedStats"]:
            if row["VariantType"] not in {"Insertion", "Deletion", "XXXDuplication"}: continue
            for b in row["PerBinStats"]:
                ev = [x for x in b["Stats"] if x["StatsType"] == "Event"][0]
                if b['Bin'].endswith("+"):
                    sz = b['Bin'][:-1]
                    b['Bin'] = f'[{sz},10000000)'
                d = dict(sizemin=b['Bin'][1:].split(',')[0],
                        sizemax=b['Bin'][:-1].split(',')[1],
                        svtype=row["VariantType"][:3].upper(),
                        precision=ev["Precision"],
                        recall=ev["Recall"],
                        name="wittyer",
                        TP=ev["QueryTpCount"],
                        FN=ev["TruthFnCount"],
                        FP=ev["QueryFpCount"],
                        f1=ev["Fscore"],
                        )
                if d['sizemin'] == '1': continue

                yield d

import svbench as svb

def svbench_evaluate(true_v, test_v, include_regions, caller, sizes=(50, 500, 5000), slop=150, name="HG002"):
    true = svb.CallSet(dataset=name).load_vcf(true_v, other_cols=["FILTER"])

    true = true.filter_include_bed(include_regions, inplace=False)
    true = true.add_intervals(slop=slop)

    test = getattr(svb, "load_" + caller)(test_v,
            dataset=name).filter_include_bed(include_regions)

    callers = [test]

    #score = svb.score(true, callers)
    #print(len(score))
    #score = score[0]
    #print(score.scores)

    inss = [i.filter_by_svtype("INS", inplace=False) for i in callers]
    dels = [i.filter_by_svtype("DEL", inplace=False) for i in callers]

    true_ins = true.filter_by_svtype("INS", inplace=False)
    true_ins = true.add_intervals(slop=slop)

    true_del = true.filter_by_svtype("DEL", inplace=False)
    true_del = true.add_intervals(slop=slop)

    trues = {"INS": true_ins, "DEL": true_del}
    tests = {"INS": inss, "DEL": dels}
    sizes = list(sizes) + [10000000]

    for ev in ("INS", "DEL"):
        svb.score(trues[ev], tests[ev], ref_size_bins=sizes,
                min_ref_size=sizes[0] - 20, stratify=True)
        sc = tests[ev][0]

        for k, row in sc.size_scores.iterrows():
            r = dict(row)
            r["precision"] = r.pop("Precision")
            r["recall"] = r.pop("Recall")
            if r["Ref size ranges"] == "All ranges": continue
            rsizes = r.pop("Ref size ranges")[1:-1].split(", ")
            r["sizemin"] = rsizes[0].strip()
            r["sizemax"] = rsizes[1].strip()
            r["svtype"] = ev
            r["f1"] = r.pop("F1")
            r["name"] = "svbench"
            yield r

def filter_filters(input_vcf, filters):
    ivcf = cyvcf2.VCF(input_vcf)
    (_, ofile) = tempfile.mkstemp(suffix=".tmp.vcf.gz")
    ovcf = cyvcf2.Writer(ofile, ivcf)
    ovcf.write_header()

    n = 0
    tot = 0

    for v in ivcf:
        tot += 1
        vf = ["PASS"] if v.FILTER is None else v.FILTER.split(";")
        if len(filters & set(vf)) == 0:
            continue
        if v.ALT[0] == "<TRA>": continue
        if v.CHROM == "hs37d5": continue

        if v.INFO.get("SVTYPE") == "BND": continue

        ovcf.write_record(v)
        n += 1

    ovcf.close()
    ivcf.close()
    print(f"[sveval] ({input_vcf}) wrote {n} out of {tot} variants matching filters", file=sys.stderr)
    subprocess.check_call(["tabix", "-f", ofile]) 
    return ofile

def wittyer(true_v, test_v, bin_sizes):

    odir = "sveval.wittyer.test"

    shutil.rmtree(odir, ignore_errors=True)
    if bin_sizes.endswith(",inf"):
        bin_sizes = bin_sizes[:-4]

    cmd = [
    "/opt/Wittyer/Wittyer",
    "-b", "/opt/data/hg002/HG002_SVs_Tier1_v0.6.bed",
    "--binSizes", f"{bin_sizes}",
    "-t", true_v,
    "-i", test_v,
    "--em", "SimpleCounting",
    "-o", odir]
    print(" ".join(cmd), file=sys.stderr)

    return (odir, subprocess.Popen(cmd,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE))

def truvari(true_v, test_v, sizemin, sizemax, fasta, svtype):

    out = f"sveval.truvari.{sizemin}.{sizemax}.{svtype}"
    shutil.rmtree(out, ignore_errors=True)

    cmd = [
            "python", "/truvari/truvari.py",
            "--passonly", "--pctsim=0", "--giabreport", "--no-ref",
            "--includebed", REGIONS,
            "--svtype", svtype,
            "-b", true_v,
            "-c", test_v,
            "-f", fasta,
            "-s", str(sizemin),
            "--sizemax", str(sizemax),
            "-o", out,
            ]
    if svtype == "INS":
        cmd.extend(["--pctsize", "0.1", "-O", "0.0", "-r", "100"])
    else:
        cmd.extend(["-O", "0.4", "-r", "20"])

    return (out, subprocess.Popen(cmd,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE))


TRUTH_VCF = "/opt/data/hg002/HG002_SVs_Tier1_v0.6.vcf.gz"
REGIONS = "/opt/data/hg002/HG002_SVs_Tier1_v0.6.bed"

def main(args=sys.argv[1:]):

    p = argparse.ArgumentParser()
    p.add_argument("--bin-sizes", default="50,500,5000,inf")
    p.add_argument("-f", "--allowed-filters", default="PASS,.,LongReadHomRef")
    p.add_argument("-c", "--caller", choices=("dysgu", "manta"), help="caller name for svbench", required=True)
    p.add_argument("--fasta", help="fasta reference", required=True)
    p.add_argument("vcf", help="path to SV vcf to evaluate")

    a = p.parse_args(args)

    filters = set(a.allowed_filters.strip().split(","))

    test_v = filter_filters(a.vcf, filters)
    true_v = filter_filters(TRUTH_VCF, filters)

    ps = []

    dirs = []
    d, p = wittyer(true_v, test_v, a.bin_sizes)
    dirs.append(d)
    ps.append(p)

    bins = [1000000000 if x == "inf" else int(x) for x in a.bin_sizes.split(",")]
    for i, b in enumerate(bins):
        if i == len(bins) - 1: break

        for svt in ("DEL", "INS"):
          d, p = truvari(true_v, test_v, b, bins[i + 1], a.fasta, svt)
          ps.append(p)
          dirs.append(d)

    for p in ps:
        if p.stderr is not None:
          p.stderr.read()
        p.wait()

    #os.unlink(test_v)
    #os.unlink(true_v)

    for p in ps:
        if p.returncode != 0:
            if p.stderr is not None:
                print(p.stderr.read(), file=sys.stderr)
            #sys.exit(p.returncode)

    fmt = "{svtype}\t{name}\t{sizemin}\t{sizemax}\t{precision}\t{recall}\t{TP}\t{FP}\t{FN}"
    print("#" + fmt.replace("{", "").replace("}", ""))
    for d in dirs:
       if "sveval.truvari" in d and os.path.exists(os.path.join(d, "summary.txt")):
           print(fmt.format(**extract_truvari(os.path.join(d, "summary.txt"))))
       elif "sveval.wittyer" in d and os.path.exists(os.path.join(d, "Wittyer.Stats.json")):
           for d in extract_wittyer(os.path.join(d, "Wittyer.Stats.json")):
              print(fmt.format(**d))

    for d in svbench_evaluate(true_v, test_v, REGIONS, a.caller, bins[:-1]):
       print(fmt.format(**d))

if __name__ == "__main__":
    main()
