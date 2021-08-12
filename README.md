# CMPP
**Annotate C(CAZy)M(MEROPS)P(PHI)P(P450) for fungi.**  
Dependency: ncbi-blast (v2.9), hmmscan, venn and [functional](wulab.cf/CMPP/funcational.tar.gz)
## Install
Clone repo to local and download functional database.  
```bash
git clone https://github.com/JinyuanSun/CMPP.git
cd CMPP
wget http://wulab.cf/CMPP/funcational.tar.gz
tar vczf functional.tar.gz
```
## Run CMPP pipeline
```bash
./run_CMPP_anno.sh -i $protein.fas -d functional
```
Check the output in `CMPP_out`.
