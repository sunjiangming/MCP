# MCCP
This refer to 18 datasets that were extracted from the ExCAPE-DB database [1], a repository storing public available chemogenomics data. Signature descriptors [2] of heights 0-3 were generated for all compounds in ExCAPE-DB. You may join the split files prior to any analysis.

## Activity output
BioActs.tar.part1</br>
BioActs.tar.part2</br>
</br>
The tar file contains 18 files. File name is actually target's name (Entrez ID).</br>
Example:</br>
Ambit_InChIKey  A(active)/N(inactive)</br>
AADCDMQTJNYOSS-YPVXKFLCNA-N     A</br>
AADWBQLQKGQDIM-UHFFFAOYNA-N     A</br>

## Compound's structure
compound.smi.part1</br>
compound.smi.part2</br>
compound.smi.part3</br>
compound.smi.part4</br></br>
Example:</br>
AMBIT_InChIKey  AMBIT_InChI     AMBIT_SMILES</br>
HAFSOLHTAKKOHZ-CMDGGOBGNA-N     InChI=1/C17H16N2O/c1-12-4-6-14(7-5-12)17(20)9-8-15-10-16(11-18)19(3)13(15)2/h4-10H,1-3H3/b9-8+  O=C(/C=C/C1=C(N(C(=C1)C#N)C)C)C2=CC=C(C=C2)C
</br>GAWSSNLFMXRZOG-AEMIGZHMNA-N     InChI=1/C14H11N5O3/c20-13-6-5-10(19(21)22)7-9(13)8-15-18-14-16-11-3-1-2-4-12(11)17-14/h1-8,15H,(H2,16,17,18)/b9-8+/f/h16,18H    O=C1/C(=C/NNC=2NC=3C(N2)=CC=CC3)/C=C(N(=O)=O)C=C1

## Compound's signatures in the libsvm format
compound2signatures.libsvm.txt.part1</br>
compound2signatures.libsvm.txt.part2</br>
compound2signatures.libsvm.txt.part3</br>
compound2signatures.libsvm.txt.part4</br>
compound2signatures.libsvm.txt.part5</br></br>
Example:</br>
HAFSOLHTAKKOHZ-CMDGGOBGNA-N     3:1.0 4:5.0 9:2.0 14:1.0 19:2.0 45:2.0 57:2.0 136:2.0 207:2.0 212:1.0 213:1.0 223:1.0 232:1.0 233:1.0 234:1.0 235:1.0 236:1.0 237:1.0 238:1.0 239:1.0 240:1.0 241:1.0 242:1.0 243:1.0 244:1.0 245:1.0 246:1.0 247:1.0 248:1.0 249:1.0 250:1.0 251:1.0 252:1.0 253:1.0 254:1.0 255:1.0 256:1.0 257:1.0 258:1.0 259:1.0 260:1.0 261:1.0 262:1.0 263:1.0 264:1.0 265:2.0 266:2.0 267:1.0
</br>GAWSSNLFMXRZOG-AEMIGZHMNA-N     4:4.0 14:1.0 97:2.0 110:2.0 136:3.0 232:1.0 237:1.0 239:1.0 268:1.0 269:1.0 270:1.0 271:2.0 272:2.0 273:1.0 274:1.0 275:1.0 276:1.0 277:1.0 278:1.0 279:1.0 280:1.0 281:1.0 282:1.0 283:2.0 284:2.0 285:2.0 286:1.0 287:1.0 288:1.0 289:1.0 290:1.0 291:1.0 292:1.0 293:1.0 294:1.0 295:1.0 296:1.0 297:1.0 298:1.0 299:2.0 300:2.0 301:2.0 302:2.0 303:1.0 304:1.0 305:1.0 306:1.0 307:1.0 308:1.0 309:1.0

## Contact</br>
Jiangming Sun, PhD</br>
AstraZeneca R&D, Gothenburg</br>
Email_1: Jiangming.Sun at astrazeneca.com</br>
Email_2: sunjiangming at gmail.com</br>

reference
[1] Sun, J.; Jeliazkova, N.; Chupakin, V.; Golib-Dzib, J.-F.; Engkvist, O.; Carlsson, L.; Wegner, J.; Ceulemans, H.; Georgiev, I.; Jeliazkov, V.; Kochev, N.; Ashby, T. J.; Chen, H., ExCAPE-DB: an integrated large scale dataset facilitating Big Data analysis in chemogenomics. Journal of Cheminformatics 2017, 9 (1), 17
[2] Carbonell, P.; Carlsson, L.; Faulon, J.-L., Stereo Signature Molecular Descriptor. Journal of Chemical Information and Modeling 2013, 53 (4), 887-897.
