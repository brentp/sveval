FROM microsoft/dotnet:2.0-sdk as wittyer
WORKDIR /src
COPY ./witty.er /src
RUN cd Wittyer \
      && dotnet publish -f netcoreapp2.0 -r linux-x64 -c Release -o /output \
      && chmod +x /output/Wittyer
      
FROM microsoft/dotnet:2.0.9-runtime as dotnet
LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
#RUN apt-get -y update && apt-get -y install tabix
#ENTRYPOINT ["/opt/Wittyer/Wittyer"]


FROM ubuntu:20.04
COPY --from=wittyer /output /opt/Wittyer
COPY --from=dotnet /opt/Wittyer /opt/Wittyer
RUN apt update && DEBIAN_FRONTEND=noninteractive apt install -qy python3 python3-pip git wget tabix libunwind8 icu-devtools bcftools && \
    ln -sf /usr/bin/python3 /usr/bin/python && ln -sf /usr/bin/pip3 /usr/bin/pip && \
    pip install cython progressbar2 python-levenshtein "intervaltree<2.1.0" \
pyvcf pyfaidx pysam python-dateutil tabulate numpy cyvcf2 matplotlib joblib && \
    git clone -b no-filter-as-pass https://github.com/brentp/truvari.git && python truvari/truvari.py -h &&   \
    git clone https://github.com/kcleal/svbench && cd svbench && pip install -r requirements.txt && pip install .


RUN mkdir -p /opt/data/hg002/ && cd /opt/data/hg002/ && \
    wget -q ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz && \
    wget -q ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed && \
    wget -q ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1plusTier2_v0.6.1.bed

ADD ./sveval.py /usr/bin/sveval.py
RUN chmod +x /usr/bin/sveval.py

ENTRYPOINT ["/usr/bin/sveval.py"]
