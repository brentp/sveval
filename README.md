# sveval: evaulate variant callers on hg002 truth set.

This uses:
+ [Wittyer](https://github.com/Illumina/witty.er)
+ [svbench](https://github.com/kcleal/svbench)
+ [truvari](https://github.com/spiralgenetics/truvari)

to evaluate structural variant calls. 

Run via docker as:

```
docker run -it -v $path:$mapped_path brentp/sveval:v0.0.1 --fasta $fasta $dysgu_vcf --caller dysgu

docker run -it -v $path:$mapped_path brentp/sveval:v0.0.1 --fasta $fasta $manta_vcf --caller manta
```

where `$vcf` **must be in GRCh37 coordinates**.

This will report a table like:
```
#svtype	name	sizemin	sizemax	precision	recall	TP	FP	FN
DEL	wittyer	50	 500	0.8861863240018357	0.566733929086192	1931	248	1503
DEL	wittyer	500	 5000	0.9617486338797814	0.8758389261744967	528	21	74
DEL	wittyer	5000	10000000	0.984	0.8840579710144928	123	2	16
INS	wittyer	50	 500	0.9257575757575758	0.32461355529131986	1222	98	2840
INS	wittyer	500	 5000	0.9626168224299065	0.15391539153915393	206	8	940
INS	wittyer	5000	10000000	NaN	0.078125	0	0	118
DEL	truvari	50	500	0.8925318761384335	0.5650043240126837	1960	236	1509
INS	truvari	50	500	0.8233240223463687	0.28024720703589256	1179	253	3028
DEL	truvari	500	5000	0.9435336976320583	0.8691275167785235	518	31	78
INS	truvari	500	5000	0.4507042253521127	0.11531531531531532	128	156	982
DEL	truvari	5000	1000000000	0.976	0.8840579710144928	122	3	16
INS	truvari	5000	1000000000	0	0	0	0	128
INS	svbench	50	500	0.910547396528705	0.1269900381714924	1364	134	9377
INS	svbench	500	5000	0.9056603773584906	0.06304728546409807	144	15	2140
INS	svbench	5000	10000000	nan	0.0	0	0	386
DEL	svbench	50	500	0.9668874172185431	0.2038916302020296	2190	75	8551
DEL	svbench	500	5000	0.8501742160278746	0.2136602451838879	488	86	1796
DEL	svbench	5000	10000000	0.9117647058823529	0.32124352331606215	124	12	262
```
