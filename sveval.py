#!/usr/bin/python
import argparse
import cyvcf2
import sys
import shutil
import subprocess
import tempfile

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

        ovcf.write_record(v)
        n += 1

    ovcf.close()
    ivcf.close()
    print(f"[sveval] ({input_vcf}) wrote {n} out of {tot} variants matching filters", file=sys.stderr)
    subprocess.check_call(["tabix", "-f", ofile]) 
    return ofile

def wittyer(true_v, test_v):

    odir = "_witty.test"

    shutil.rmtree(odir, ignore_errors=True)

    cmd = [
    "/opt/Wittyer/Wittyer",
    "-b", "/opt/data/hg002/HG002_SVs_Tier1_v0.6.bed",
    "--binSizes", '!1,30,50,500,5000',
    "-t", true_v,
    "-i", test_v,
    "--em", "SimpleCounting",
    "-o", odir]

    return subprocess.Popen(cmd, stderr=subprocess.PIPE, stdout=sys.stdout)

def truvari(true_v, test_v, sizemin, sizemax, fasta, svtype):

    out = f"sveval.truvari.{sizemin}.{sizemax}.{svtype}"
    shutil.rmtree(out, ignore_errors=True)

    cmd = [
            "python", "/truvari/truvari.py",
            "--passonly", "--pctsim=0", "--giabreport", "--no-ref",
            "--includebed", REGIONS,
            "-r", "20",
            "--svtype", svtype,
            "-O", "0.4",
            "-b", true_v,
            "-c", test_v,
            "-f", fasta,
            "-s", str(sizemin),
            "--sizemax", str(sizemax),
            "-o", out,
            ]

    return subprocess.Popen(cmd, stderr=subprocess.PIPE)


TRUTH_VCF = "/opt/data/hg002/HG002_SVs_Tier1_v0.6.vcf.gz"
REGIONS = "/opt/data/hg002/HG002_SVs_Tier1_v0.6.bed"

def main(args=sys.argv[1:]):

    p = argparse.ArgumentParser()
    p.add_argument("--bin-sizes", default="50,500,5000,inf")
    p.add_argument("-f", "--allowed-filters", default="PASS,.,LongReadHomRef")
    p.add_argument("--fasta", help="fasta reference", required=True)
    p.add_argument("vcf", help="path to SV vcf to evaluate")

    a = p.parse_args(args)

    filters = set(a.allowed_filters.strip().split(","))

    test_v = filter_filters(a.vcf, filters)
    true_v = filter_filters(TRUTH_VCF, filters)

    ps = []
    ps.append(wittyer(true_v, test_v))

    bins = [1000000000 if x == "inf" else int(x) for x in a.bin_sizes.split(",")]
    for i, b in enumerate(bins):
        if i == len(bins) - 1: break

        ps.append(truvari(true_v, test_v, b, bins[i + 1], a.fasta, "DEL"))
        ps.append(truvari(true_v, test_v, b, bins[i + 1], a.fasta, "INS"))

    for p in ps:
        p.wait()

    #os.unlink(test_v)
    #os.unlink(true_v)

    for p in ps:
        if p.returncode != 0:
            print(p.stderr.read())
            sys.exit(p.returncode)





if __name__ == "__main__":
    main()
