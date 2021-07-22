FROM microsoft/dotnet:2.0-sdk as wittyer
WORKDIR /src
COPY ./witty.er /src
RUN cd Wittyer \
      && dotnet publish -f netcoreapp2.0 -r linux-x64 -c Release -o /output \
      && chmod +x /output/Wittyer
      
FROM microsoft/dotnet:2.0.9-runtime
LABEL git_repository=https://git.illumina.com/DASTE/Ilmn.Das.App.Wittyer.git
WORKDIR /opt/Wittyer
RUN apt-get -y update && apt-get -y install tabix
COPY --from=wittyer /output /opt/Wittyer
#ENTRYPOINT ["/opt/Wittyer/Wittyer"]


FROM ubuntu:20.04 as truvari
RUN apt update && apt install -qy python3 python3-pip git && \
    ln -sf /usr/bin/python3 /usr/bin/python && ln -sf /usr/bin/pip3 /usr/bin/pip && \
    pip install cython progressbar2 python-levenshtein "intervaltree<2.1.0" pyvcf pyfaidx pysam python-dateutil tabulate numpy && \
    git clone -b no-filter-as-pass https://github.com/brentp/truvari.git && python truvari/truvari.py -h && \
    git clone https://github.com/kcleal/svbench && cd svbench && pip install -r requirements.txt && pip install .
